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
#include "../src/starfix_nmea.h"
#include "../src/starfix_photo.h"
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

/* pre-allocated memory for starfix algorithms */
static uint8_t mempool[10 * 1024 * 1024];
static starfix_arena_t arena = {mempool, mempool + sizeof(mempool)};
static uint32_t visited_mask[CAM_WIDTH * CAM_HEIGHT];

/* compute greenwich hour angle of aries in degrees from unix UTC timestamp */
static double current_gha_aries_utc(double t_now) {
    double jd = t_now / 86400.0 + 2440587.5;
    double d_ut1 = jd - 2451545.0;
    double t = d_ut1 / 36525.0;
    double gmst =
        280.46061837 + 360.98564736629 * d_ut1 + 0.000387933 * t * t - (t * t * t) / 38710000.0;
    gmst = fmod(gmst, 360.0);
    if (gmst < 0.0) gmst += 360.0;
    return gmst;
}

int main(void) {
    printf("StarFix Prototype: IMX462 + BNO055 + Pi Zero W\n");

    /* initialize hardware */
    if (camera_init(CAMERA_DEV, CAM_WIDTH, CAM_HEIGHT) < 0) {
        printf("Failed to initialize camera on %s\n", CAMERA_DEV);
        return 1;
    }
    printf("Camera initialized successfully.\n");

    if (bno055_init(I2C_DEV, BNO055_I2C_ADDR1) < 0) {
        printf("Failed to initialize BNO055 on %s (trying alt addr)\n", I2C_DEV);
        if (bno055_init(I2C_DEV, BNO055_I2C_ADDR2) < 0) {
            printf("Failed to initialize BNO055 entirely.\n");
            camera_close();
            return 1;
        }
    }
    printf("BNO055 IMU initialized successfully.\n");

    /* load starfix database */
    starfix_db_t db = {0};

    printf("Loading StarFix database...\n");
    if (starfix_load_db("../data/starfix_db.bin", &db) != 0 &&
        starfix_load_db("data/starfix_db.bin", &db) != 0) {
        printf("Failed to load star database.\n");
        bno055_close();
        camera_close();
        return 1;
    }

    /* propagate astrometry to current epoch */
    time_t t_now = time(NULL);
    struct tm* tm_now = gmtime(&t_now);
    double current_year = tm_now ? 1900.0 + tm_now->tm_year + (tm_now->tm_yday / 365.25) : 2026.0;
    starfix_propagate_astrometry(db.stars, db.num_stars, current_year);

    /* EKF state initialized on first valid celestial fix */
    starfix_ekf_t ekf;
    int ekf_initialized = 0;

    /* main run loop */
    printf("Entering main tracking loop...\n");
    int frame_count = 0;
    while (frame_count < 10) {
        uint8_t* frame_buffer = NULL;
        double captured_utc;
        if (camera_capture_frame(&frame_buffer, &captured_utc) < 0) {
            printf("Frame drop.\n");
            continue;
        }

        /* 1. read IMU zenith vector and compass heading */
        double gx = 0.0, gy = 0.0, gz = 1.0;
        if (bno055_read_gravity(&gx, &gy, &gz) < 0) {
            printf("IMU read failed.\n");
            continue;
        }
        starfix_vector3_t imu_zenith = {gx, gy, gz};

        double heading_deg = 0.0;
        if (ekf_initialized && bno055_read_heading(&heading_deg) == 0) {
            starfix_ekf_predict(&ekf, 0.0, heading_deg, 1.0);
        }

        /* 2. find centroids */
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

        /* 3. identify stars */
        starfix_match_t matches[MAX_MATCHES];
        int num_matches = starfix_identify_stars(num_centroids, centroids, CAM_WIDTH, CAM_HEIGHT,
                                                 FOV_DEG, &db, MAX_MATCHES, matches, arena, NULL);

        if (num_matches < 3) {
            printf("Frame %d: Star ID failed (%d)\n", frame_count, num_matches);
            frame_count++;
            continue;
        }

        /* 4. solve attitude and position using IMU zenith */
        starfix_photo_options_t options = {
            .width = CAM_WIDTH,
            .height = CAM_HEIGHT,
            .fov_deg = FOV_DEG,
            .gha_aries = current_gha_aries_utc(captured_utc),
            .zenith_camera = imu_zenith,
        };
        starfix_photo_result_t result;
        starfix_status_t pos_stat =
            starfix_solve_photo(num_centroids, centroids, db.num_stars, db.stars, num_matches,
                                matches, &options, &result, arena);

        if (pos_stat == STARFIX_SUCCESS) {
            if (!ekf_initialized) {
                starfix_ekf_init(&ekf, result.latitude, result.longitude);
                ekf_initialized = 1;
            } else {
                starfix_ekf_correct(&ekf, result.latitude, result.longitude);
            }
            printf(
                "Frame %d: Pos FIX [Lat: %.4f, Lon: %.4f] EKF [Lat: %.4f, Lon: %.4f] "
                "(Inliers: %u)\n",
                frame_count, result.latitude, result.longitude, ekf.lat, ekf.lon, result.inliers);
            char nmea_buf[256];
            double heading_deg = NAN;
            if (bno055_read_heading(&heading_deg) < 0) heading_deg = NAN;
            if (starfix_format_nmea_fix(ekf.lat, ekf.lon, (time_t)captured_utc, 1,
                                        (int)result.inliers, NAN, heading_deg, nmea_buf,
                                        sizeof(nmea_buf)) == STARFIX_SUCCESS) {
                fputs(nmea_buf, stdout);
            }
        } else {
            printf("Frame %d: Position fix failed.\n", frame_count);
        }

        frame_count++;
    }

    /* cleanup */
    bno055_close();
    camera_close();
    starfix_free_db(&db);

    printf("Prototype execution completed.\n");
    return 0;
}
