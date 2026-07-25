#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "starfix_arena.h"
#include "starfix_status.h"

static uint8_t mempool[10 * 1024 * 1024];
static starfix_arena_t arena;
#define RESET_ARENA() (arena.beg = mempool, arena.end = mempool + sizeof(mempool))

#include "starfix_attitude.h"

/* test runner to verify the C attitude solver under multiple conditions */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main() {
    starfix_telemetry_t telem = {0};
    (void)telem;
    printf("--- Running Attitude Solver Comprehensive Tests ---\n");

    starfix_vector3_t w_dummy[2] = {{1, 0, 0}, {0, 1, 0}};
    starfix_vector3_t v_dummy[2] = {{1, 0, 0}, {0, 1, 0}};
    starfix_quaternion_t q_est;
    double R_est[3][3];

    /* 1. NULL Pointer & Invalid Size Checks */
    printf("Test 1: Boundary check... ");
    int ret1 = starfix_solve_attitude(1, w_dummy, v_dummy, &q_est, R_est);
    int ret2 = starfix_solve_attitude(2, NULL, v_dummy, &q_est, R_est);

    if (ret1 < 0 && ret2 < 0) {
        printf("Passed\n");
    } else {
        printf("FAILED (got %d, %d)\n", ret1, ret2);
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

    int status = starfix_solve_attitude(5, w, v, &q_est, R_est);
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
            printf("FAILED (large error: %f deg)\n", max_err_deg);
        }
    } else {
        printf("FAILED (solver returned %d)\n", status);
    }

    /* 3. Collinear/Degenerate Vectors Check */
    printf("Test 3: Collinear degenerate vectors... ");
    /* define two identical/collinear vectors (attitude not fully observable) */
    starfix_vector3_t v_degen[2] = {{1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}};
    starfix_vector3_t w_degen[2] = {{0.0, 1.0, 0.0}, {0.0, 1.0, 0.0}};

    /* davenport solver should execute and return a valid eigenvector (status 0) */
    /* but attitude might be arbitrary or under-constrained */
    int status_degen = starfix_solve_attitude(2, w_degen, v_degen, &q_est, R_est);
    if (status_degen == 0) {
        printf("Passed (handled gracefully)\n");
    } else {
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
    starfix_status_t ransac_status =
        (RESET_ARENA(),
         starfix_solve_attitude_ransac(6, w_ransac, v_ransac, 0.2, &q_est, R_est, &arena, &telem));
    int ransac_inliers = (int)telem.ransac_inliers;

    if (ransac_status == STARFIX_SUCCESS && ransac_inliers == 4) {
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
            printf("FAILED (rotation error too large: %f deg)\n", max_err_deg);
        }
    } else {
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

    /* Expected:
       A = 60.2 arcseconds * PI / (180.0 * 3600.0) ~ 2.9186e-4
       R = A * tan(10 degrees) ~ 5.146e-5
       theta_t = theta_a + R */
    double A_val = 60.2 * M_PI / (3600.0 * 180.0);
    double R_val = A_val * tan(theta_a);
    double theta_t = theta_a + R_val;

    double expected_x = sin(theta_t);
    double expected_z = cos(theta_t);

    if (ref_status == 0 && fabs(w_refraction.x - expected_x) < 1e-7 &&
        fabs(w_refraction.z - expected_z) < 1e-7) {
        printf("Passed\n");
    } else {
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
        printf("FAILED (status=%d, x=%f (exp=%f), y=%f (exp=%f))\n", ab_status, v_aberration[0].x,
               expected_ab_x, v_aberration[0].y, expected_ab_y);
    }

    return 0;
}
