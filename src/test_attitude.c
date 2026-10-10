#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "starfix_arena.h"
#include "starfix_astrometry.h"
#include "starfix_attitude.h"
#include "starfix_status.h"

static uint8_t mempool[10 * 1024 * 1024];
static starfix_arena_t arena = {mempool, mempool + sizeof(mempool)};

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void) {
    int failures = 0;
    printf("--- Running Attitude Solver Comprehensive Tests ---\n");

    starfix_vector3_t w_dummy[2] = {{1, 0, 0}, {0, 1, 0}};
    starfix_vector3_t v_dummy[2] = {{1, 0, 0}, {0, 1, 0}};

    double R_est[3][3];

    /* 1. NULL Pointer & Invalid Size Checks */
    printf("Test 1: Boundary check... ");
    int ret1 = starfix_solve_attitude(1, w_dummy, v_dummy, R_est);
    int ret2 = starfix_solve_attitude(2, NULL, v_dummy, R_est);
    int ret3 = starfix_solve_attitude_ransac(2, NULL, v_dummy, 0.1, R_est, arena);

    if (ret1 < 0 && ret2 < 0 && ret3 < 0) {
        printf("Passed\n");
    } else {
        failures++;
        printf("FAILED (got %d, %d, %d)\n", ret1, ret2, ret3);
    }

    /* 2. Happy Path Attitude Estimation */
    printf("Test 2: Happy path attitude estimation... ");
    starfix_vector3_t v[5] = {{1.0, 0.0, 0.0},
                              {0.0, 1.0, 0.0},
                              {0.0, 0.0, 1.0},
                              {0.577350269, 0.577350269, 0.577350269},
                              {0.707106781, 0.707106781, 0.0}};

    double yaw = 30.0 * M_PI / 180.0;
    double pitch = 10.0 * M_PI / 180.0;
    double R_true[3][3];
    R_true[0][0] = cos(pitch) * cos(yaw);
    R_true[0][1] = cos(pitch) * sin(yaw);
    R_true[0][2] = -sin(pitch);
    R_true[1][0] = -sin(yaw);
    R_true[1][1] = cos(yaw);
    R_true[1][2] = 0.0;
    R_true[2][0] = sin(pitch) * cos(yaw);
    R_true[2][1] = sin(pitch) * sin(yaw);
    R_true[2][2] = cos(pitch);

    starfix_vector3_t w[5];
    int i;
    for (i = 0; i < 5; i++) {
        w[i].x = R_true[0][0] * v[i].x + R_true[0][1] * v[i].y + R_true[0][2] * v[i].z;
        w[i].y = R_true[1][0] * v[i].x + R_true[1][1] * v[i].y + R_true[1][2] * v[i].z;
        w[i].z = R_true[2][0] * v[i].x + R_true[2][1] * v[i].y + R_true[2][2] * v[i].z;

        /* add small noise */
        w[i].x += 0.001;
        w[i].y -= 0.001;
        w[i].z += 0.0005;

        double norm = sqrt(w[i].x * w[i].x + w[i].y * w[i].y + w[i].z * w[i].z);
        w[i].x /= norm;
        w[i].y /= norm;
        w[i].z /= norm;
    }

    int status = starfix_solve_attitude(5, w, v, R_est);
    if (status == 0) {
        double max_err = 0.0;
        int r, c;
        for (r = 0; r < 3; r++) {
            for (c = 0; c < 3; c++) {
                double err = fabs(R_est[r][c] - R_true[r][c]);
                if (err > max_err) {
                    max_err = err;
                }
            }
        }
        double max_err_deg = max_err * 180.0 / M_PI;
        if (max_err_deg < 0.1) {
            printf("Passed\n");
            printf("  Max rotation error: %f degrees\n", max_err_deg);
        } else {
            failures++;
            printf("FAILED (large error: %f deg)\n", max_err_deg);
        }
    } else {
        failures++;
        printf("FAILED (solver returned %d)\n", status);
    }

    /* 3. Collinear/Degenerate Vectors Check */
    printf("Test 3: Collinear degenerate vectors... ");
    /* define two identical/collinear vectors (attitude not fully observable) */
    starfix_vector3_t v_degen[2] = {{1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}};
    starfix_vector3_t w_degen[2] = {{0.0, 1.0, 0.0}, {0.0, 1.0, 0.0}};

    /* davenport solver should execute and return a valid eigenvector (status 0) */
    /* but attitude might be arbitrary or under-constrained */
    int status_degen = starfix_solve_attitude(2, w_degen, v_degen, R_est);
    if (status_degen == 0) {
        printf("Passed (handled gracefully)\n");
    } else {
        failures++;
        printf("FAILED (returned error code: %d)\n", status_degen);
    }

    /* 4. RANSAC Outlier Rejection Check */
    printf("Test 4: RANSAC outlier rejection... ");
    srand(42); /* seed random number generator for RANSAC reproducibility */

    starfix_vector3_t w_ransac[6];
    starfix_vector3_t v_ransac[6];

    /* 4 inliers */
    for (i = 0; i < 4; i++) {
        v_ransac[i] = v[i];
        w_ransac[i].x = R_true[0][0] * v[i].x + R_true[0][1] * v[i].y + R_true[0][2] * v[i].z;
        w_ransac[i].y = R_true[1][0] * v[i].x + R_true[1][1] * v[i].y + R_true[1][2] * v[i].z;
        w_ransac[i].z = R_true[2][0] * v[i].x + R_true[2][1] * v[i].y + R_true[2][2] * v[i].z;
    }

    /* 2 outliers (corrupted values) */
    v_ransac[4] = v[4];
    w_ransac[4].x = 1.0;
    w_ransac[4].y = 0.0;
    w_ransac[4].z = 0.0; /* pointing in a completely different direction */

    v_ransac[5].x = 0.0;
    v_ransac[5].y = 0.0;
    v_ransac[5].z = 1.0;
    w_ransac[5].x = 0.0;
    w_ransac[5].y = 1.0;
    w_ransac[5].z = 0.0;

    /* run RANSAC attitude solver */
    int ransac_inliers = starfix_solve_attitude_ransac(6, w_ransac, v_ransac, 0.2, R_est, arena);

    if (ransac_inliers == 4) {
        double max_err = 0.0;
        int r, c;
        for (r = 0; r < 3; r++) {
            for (c = 0; c < 3; c++) {
                double err = fabs(R_est[r][c] - R_true[r][c]);
                if (err > max_err) {
                    max_err = err;
                }
            }
        }
        double max_err_deg = max_err * 180.0 / M_PI;
        if (max_err_deg < 0.1) {
            printf("Passed (found %d inliers)\n", ransac_inliers);
        } else {
            failures++;
            printf("FAILED (rotation error too large: %f deg)\n", max_err_deg);
        }
    } else {
        failures++;
        printf("FAILED (got %d inliers, expected 4)\n", ransac_inliers);
    }

    /* 5. Atmospheric Refraction Correction Check */
    printf("Test 5: Atmospheric Refraction Correction... ");
    starfix_vector3_t z_cam = {0.0, 0.0, 1.0};
    double theta_a = 10.0 * M_PI / 180.0;
    starfix_vector3_t w_refraction;
    w_refraction.x = sin(theta_a);
    w_refraction.y = 0.0;
    w_refraction.z = cos(theta_a);

    int ref_status = starfix_correct_refraction(1, &w_refraction, &z_cam, 10.0, 1010.0);

    double apparent_alt = M_PI / 2.0 - theta_a;
    double true_alt;
    starfix_unrefract(1010.0, 10.0, apparent_alt, &true_alt);
    double theta_t = M_PI / 2.0 - true_alt;

    double expected_x = sin(theta_t);
    double expected_z = cos(theta_t);

    if (ref_status == 0 && fabs(w_refraction.x - expected_x) < 1e-12 &&
        fabs(w_refraction.z - expected_z) < 1e-12) {
        printf("Passed\n");
    } else {
        failures++;
        printf("FAILED (status=%d, x=%f (exp=%f), z=%f (exp=%f))\n", ref_status, w_refraction.x,
               expected_x, w_refraction.z, expected_z);
    }

    /* 6. Stellar Aberration Correction Check */
    printf("Test 6: Stellar Aberration Correction... ");
    starfix_vector3_t v_aberration[1] = {{1.0, 0.0, 0.0}};
    double v_obs_km_s[3] = {0.0, 30.0, 0.0}; /* ~30 km/s Earth velocity */

    int ab_status = starfix_correct_aberration(1, v_aberration, v_obs_km_s);
    double c_kms = 299792.458;
    double beta_y = 30.0 / c_kms;
    double norm_ab = sqrt(1.0 + beta_y * beta_y);
    double expected_ab_x = 1.0 / norm_ab;
    double expected_ab_y = beta_y / norm_ab;

    if (ab_status == 0 && fabs(v_aberration[0].x - expected_ab_x) < 1e-7 &&
        fabs(v_aberration[0].y - expected_ab_y) < 1e-7) {
        printf("Passed\n");
    } else {
        failures++;
        printf("FAILED (status=%d, x=%f (exp=%f), y=%f (exp=%f))\n", ab_status, v_aberration[0].x,
               expected_ab_x, v_aberration[0].y, expected_ab_y);
    }

    printf("Test 7: Shared geometry conventions... ");
    double focal = starfix_focal_length(200, 90);
    starfix_vector3_t center = starfix_camera_ray(0, 0, focal);
    starfix_vector3_t corner = starfix_camera_ray(100, -100, focal);
    starfix_vector3_t extreme = starfix_camera_ray(DBL_MAX / 4, 0, 1);
    starfix_vector3_t east = starfix_celestial_ray(M_PI / 2, 0);
    starfix_vector3_t pole = starfix_celestial_ray(0, M_PI / 2);
    const double cyclic[3][3] = {{0, 0, 1}, {1, 0, 0}, {0, 1, 0}};
    starfix_vector3_t rotated = starfix_rotate(cyclic, (starfix_vector3_t){1, 2, 3});
    if (fabs(focal - 100) < 1e-12 && fabs(center.x) < 1e-12 && fabs(center.y) < 1e-12 &&
        fabs(center.z - 1) < 1e-12 && fabs(corner.x - 1 / sqrt(3.0)) < 1e-12 &&
        fabs(corner.y + 1 / sqrt(3.0)) < 1e-12 && fabs(corner.z - 1 / sqrt(3.0)) < 1e-12 &&
        fabs(extreme.x - 1) < 1e-12 && isfinite(extreme.z) && extreme.z > 0 &&
        fabs(east.x) < 1e-12 && fabs(east.y - 1) < 1e-12 && fabs(east.z) < 1e-12 &&
        fabs(pole.x) < 1e-12 && fabs(pole.z - 1) < 1e-12 && fabs(rotated.x - 3) < 1e-12 &&
        fabs(rotated.y - 1) < 1e-12 && fabs(rotated.z - 2) < 1e-12 &&
        fabs(starfix_dot(rotated, rotated) - 14) < 1e-12) {
        printf("Passed\n");
    } else {
        failures++;
        printf("FAILED (projection, rotation, or normalization)\n");
    }

    return failures ? 1 : 0;
}
