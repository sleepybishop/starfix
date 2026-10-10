#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "../src/starfix_arena.h"
#include "../src/starfix_attitude.h"
#include "../src/starfix_centroid.h"
#include "../src/starfix_ekf.h"
#include "../src/starfix_identify.h"
#include "../src/starfix_solver.h"
#include "../src/starfix_status.h"
#include "bno055.h"
#include "camera.h"

#define CAMERA_DEV "/dev/video0"
#define I2C_DEV "/dev/i2c-1"

#define CAM_WIDTH 1920
#define CAM_HEIGHT 1080
#define FOV_DEG 12.0

#define MAX_CENTROIDS 100
#define MAX_MATCHES 20

/* Pre-allocated memory for starfix algorithms */
static uint8_t mempool[10 * 1024 * 1024];
static starfix_arena_t arena = {mempool, mempool + sizeof(mempool)};
static uint32_t visited_mask[CAM_WIDTH * CAM_HEIGHT];

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
    starfix_db_t db = {0};

    printf("Loading StarFix database...\n");
    if (starfix_load_db("../data/starfix_db.bin", &db) != 0) {
        printf("Failed to load star database.\n");
        return 1;
    }

    /* Propagate precession */
    time_t t_now = time(NULL);
    struct tm* tm_now = gmtime(&t_now);
    double current_year = 1900.0 + tm_now->tm_year + (tm_now->tm_yday / 365.25);
    starfix_propagate_astrometry(db.stars, db.num_stars, current_year);

    /* Initialize EKF */
    starfix_ekf_t ekf;
    starfix_ekf_init(&ekf, 0.0, 0.0); /* Init at 0,0 for testing */

    double focal_length = starfix_focal_length(CAM_WIDTH, FOV_DEG);

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
        starfix_centroid_t centroids[MAX_CENTROIDS];
        int num_centroids =
            starfix_find_centroids(frame_buffer, CAM_WIDTH, CAM_HEIGHT, 45, MAX_CENTROIDS,
                                   centroids, visited_mask, (uint32_t)frame_count + 1U);

        if (num_centroids < 4) {
            printf("Frame %d: Centroid detection returned %d, skipping.\n", frame_count,
                   num_centroids);
            frame_count++;
            continue;
        }

        /* 3. Identify Stars */
        starfix_match_t matches[MAX_MATCHES];
        int num_matches = starfix_identify_stars(num_centroids, centroids, CAM_WIDTH, CAM_HEIGHT,
                                                 FOV_DEG, &db, MAX_MATCHES, matches, arena, NULL);

        if (num_matches < 3) {
            printf("Frame %d: Star ID failed (%d)\n", frame_count, num_matches);
            frame_count++;
            continue;
        }

        /* 4. Estimate Attitude */
        starfix_vector3_t cam_vecs[MAX_MATCHES];
        starfix_vector3_t cat_vecs[MAX_MATCHES];
        for (int i = 0; i < num_matches; i++) {
            int c_idx = matches[i].centroid_idx;
            int cat_idx = matches[i].catalog_idx;

            cam_vecs[i] = starfix_camera_ray(centroids[c_idx].u - CAM_WIDTH / 2.0,
                                             centroids[c_idx].v - CAM_HEIGHT / 2.0, focal_length);
            cat_vecs[i] = starfix_celestial_ray(db.stars[cat_idx].ra, db.stars[cat_idx].dec);
        }

        double R_est[3][3];

        int inliers =
            starfix_solve_attitude_ransac(num_matches, cam_vecs, cat_vecs, 0.1, R_est, arena);

        if (inliers < 0) {
            printf("Frame %d: Attitude solve failed.\n", frame_count);
            frame_count++;
            continue;
        }

        /* 5. Solve Position using IMU Zenith */
        double solved_lat, solved_lon;
        double gha_aries = 100.0; /* Mock value for prototype */

        starfix_status_t pos_stat =
            starfix_solve_position(1, &imu_zenith, R_est, &gha_aries, &solved_lat, &solved_lon);

        if (pos_stat == STARFIX_SUCCESS) {
            starfix_ekf_correct(&ekf, solved_lat, solved_lon);
            printf("Frame %d: Pos FIX [Lat: %.4f, Lon: %.4f] (Inliers: %d)\n", frame_count,
                   solved_lat, solved_lon, inliers);
        } else {
            printf("Frame %d: Position fix failed.\n", frame_count);
        }

        frame_count++;
    }

    /* Cleanup */
    camera_close();
    starfix_free_db(&db);

    printf("Prototype execution completed.\n");
    return 0;
}
