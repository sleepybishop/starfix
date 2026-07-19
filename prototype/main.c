#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>

#include "bno055.h"
#include "camera.h"

#include "../src/starfix_arena.h"
#include "../src/starfix_centroid.h"
#include "../src/starfix_identify.h"
#include "../src/starfix_attitude.h"
#include "../src/starfix_solver.h"
#include "../src/starfix_ekf.h"
#include "../src/starfix_status.h"

#define CAMERA_DEV "/dev/video0"
#define I2C_DEV "/dev/i2c-1"

#define CAM_WIDTH 1920
#define CAM_HEIGHT 1080
#define FOV_DEG 12.0

#define MAX_CENTROIDS 100
#define MAX_MATCHES 20

/* Pre-allocated memory for starfix algorithms */
static uint8_t mempool[10 * 1024 * 1024];
static starfix_arena_t arena;
#define RESET_ARENA() (arena.beg = mempool, arena.end = mempool + sizeof(mempool))

int main(void) {
    printf("StarFix Prototype: IMX462 + BNO055 + Pi Zero W\n");

    /* Initialize hardware */
    if (camera_init(CAMERA_DEV, CAM_WIDTH, CAM_HEIGHT) < 0) {
        printf("Failed to initialize camera on %s\n", CAMERA_DEV);
        return 1;
    }
    printf("Camera initialized successfully.\n");

    if (bno055_init(I2C_DEV, BNO055_I2C_ADDR1) < 0) {
        printf("Failed to initialize BNO055 on %s (trying alt addr)\n", I2C_DEV);
        if (bno055_init(I2C_DEV, BNO055_I2C_ADDR2) < 0) {
            printf("Failed to initialize BNO055 entirely.\n");
            return 1;
        }
    }
    printf("BNO055 IMU initialized successfully.\n");

    /* Load StarFix Database */
    starfix_catalog_star_t* catalog = NULL;
    starfix_hash_entry_t* hash_table = NULL;
    uint32_t num_stars = 0, num_entries = 0, bin_factor = 0;
    
    printf("Loading StarFix database...\n");
    if (starfix_load_db("../data/starfix_db.bin", &catalog, &num_stars, &hash_table, &num_entries, &bin_factor) != 0) {
        printf("Failed to load star database.\n");
        return 1;
    }
    
    /* Propagate precession */
    time_t t_now = time(NULL);
    struct tm* tm_now = gmtime(&t_now);
    double current_year = 1900.0 + tm_now->tm_year + (tm_now->tm_yday / 365.25);
    starfix_propagate_precession(catalog, num_stars, current_year);

    /* Initialize EKF */
    starfix_ekf_t ekf;
    starfix_ekf_init(&ekf, 0.0, 0.0); /* Init at 0,0 for testing */
    
    double focal_length = CAM_WIDTH / (2.0 * tan((FOV_DEG * M_PI / 180.0) / 2.0));

    /* Main run loop */
    printf("Entering main tracking loop...\n");
    int frame_count = 0;
    while (frame_count < 10) { /* Run 10 frames for prototype test */
        uint8_t* frame_buffer = NULL;
        if (camera_capture_frame(&frame_buffer) < 0) {
            printf("Frame drop.\n");
            continue;
        }

        /* 1. Read IMU Zenith Vector */
        double gx, gy, gz;
        if (bno055_read_gravity(&gx, &gy, &gz) < 0) {
            printf("IMU read failed.\n");
            continue;
        }
        starfix_vector3_t imu_zenith = {gx, gy, gz};

        /* 2. Find Centroids */
        starfix_telemetry_t telem = {0};
        starfix_centroid_t centroids[MAX_CENTROIDS];
        starfix_status_t c_stat = starfix_find_centroids(frame_buffer, CAM_WIDTH, CAM_HEIGHT, 
                                            45, MAX_CENTROIDS, centroids, &telem);
        
        if (c_stat != STARFIX_SUCCESS || telem.num_stars_detected < 4) {
            printf("Frame %d: Found only %d stars, skipping.\n", frame_count, telem.num_stars_detected);
            frame_count++;
            continue;
        }

        /* 3. Identify Stars */
        RESET_ARENA();
        starfix_match_t matches[MAX_MATCHES];
        starfix_status_t id_stat = starfix_identify_stars(telem.num_stars_detected, centroids, CAM_WIDTH, CAM_HEIGHT, 
                                            FOV_DEG, num_stars, catalog, num_entries, hash_table, 
                                            bin_factor, MAX_MATCHES, matches, &arena, &telem);
        
        if (id_stat != STARFIX_SUCCESS || telem.identify_matches < 3) {
            printf("Frame %d: Star ID failed. Matches: %d\n", frame_count, telem.identify_matches);
            frame_count++;
            continue;
        }

        /* 4. Estimate Attitude */
        starfix_vector3_t cam_vecs[MAX_MATCHES];
        starfix_vector3_t cat_vecs[MAX_MATCHES];
        for (uint32_t i = 0; i < telem.identify_matches; i++) {
            int c_idx = matches[i].centroid_idx;
            int cat_idx = matches[i].catalog_idx;

            double cx = (centroids[c_idx].u - CAM_WIDTH / 2.0) / focal_length;
            double cy = -(centroids[c_idx].v - CAM_HEIGHT / 2.0) / focal_length;
            double norm = sqrt(cx*cx + cy*cy + 1.0);
            cam_vecs[i].x = cx/norm; cam_vecs[i].y = cy/norm; cam_vecs[i].z = 1.0/norm;

            double dec_s = starfix_catalog_dec(&catalog[cat_idx]);
            double ra_s = starfix_catalog_ra(&catalog[cat_idx]);
            cat_vecs[i].x = cos(dec_s)*cos(ra_s);
            cat_vecs[i].y = cos(dec_s)*sin(ra_s);
            cat_vecs[i].z = sin(dec_s);
        }

        starfix_quaternion_t q_est;
        double R_est[3][3];
        RESET_ARENA();
        starfix_status_t att_stat = starfix_solve_attitude_ransac(telem.identify_matches, cam_vecs, cat_vecs, 
                                                    0.1, &q_est, R_est, &arena, &telem);
        
        if (att_stat != STARFIX_SUCCESS) {
            printf("Frame %d: Attitude solve failed.\n", frame_count);
            frame_count++;
            continue;
        }

        /* 5. Solve Position using IMU Zenith */
        double solved_lat, solved_lon;
        double gha_aries = 100.0; /* Mock value for prototype */
        RESET_ARENA();
        starfix_status_t pos_stat = starfix_solve_position(1, &imu_zenith, R_est, &gha_aries, ekf.lat, ekf.lon, 
                                            &solved_lat, &solved_lon, &arena, &telem);
        
        if (pos_stat == STARFIX_SUCCESS) {
            starfix_ekf_correct(&ekf, solved_lat, solved_lon, &telem);
            printf("Frame %d: Pos FIX [Lat: %.4f, Lon: %.4f] (Inliers: %d)\n", 
                   frame_count, solved_lat, solved_lon, telem.ransac_inliers);
        } else {
            printf("Frame %d: Position fix failed.\n", frame_count);
        }

        frame_count++;
    }

    /* Cleanup */
    camera_close();
    free(catalog);
    free(hash_table);

    printf("Prototype execution completed.\n");
    return 0;
}
