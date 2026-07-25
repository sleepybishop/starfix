#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "starfix_arena.h"
#include "starfix_status.h"

static uint8_t mempool[10 * 1024 * 1024];
static starfix_arena_t arena;
#define RESET_ARENA() (arena.beg = mempool, arena.end = mempool + sizeof(mempool))

#include "starfix_solver.h"

/* test runner to verify the C position solver under multiple conditions */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main() {
    printf("--- Running Position Solver Comprehensive Tests ---\n");

    const double R_dummy[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    double gha_dummy = 120.0;
    double solved_lat = 0.0, solved_lon = 0.0;

    /* 1. Invalid Arguments & NULL Pointer Robustness */
    printf("Test 1: Argument validation... ");
    int ret1 = (RESET_ARENA(), starfix_solve_position(0, NULL, R_dummy, &gha_dummy, 34.0, -118.0,
                                                      &solved_lat, &solved_lon, &arena, NULL));
    int ret2 = (RESET_ARENA(), starfix_solve_position(1, NULL, R_dummy, &gha_dummy, 34.0, -118.0,
                                                      &solved_lat, &solved_lon, &arena, NULL));

    if (ret1 < 0 && ret2 < 0) {
        printf("Passed\n");
    } else {
        printf("FAILED (got %d, %d)\n", ret1, ret2);
    }

    /* 2. Happy Path Fix Solver */
    printf("Test 2: Happy path position fix solver... ");
    double true_lat = 34.0522;
    double true_lon = -118.2437;
    double ap_lat = 34.0;
    double ap_lon = -118.0;
    double gha_aries = 120.0;

    /* simulated camera attitude: RA = 83.82, Dec = -5.39, Roll = 15.0 deg */
    double cam_ra = 83.82 * M_PI / 180.0;
    double cam_dec = -5.39 * M_PI / 180.0;
    double cam_roll = 15.0 * M_PI / 180.0;

    double z_cam[3] = {cos(cam_dec) * cos(cam_ra), cos(cam_dec) * sin(cam_ra), sin(cam_dec)};
    double x_zero[3] = {-sin(cam_ra), cos(cam_ra), 0.0};
    double y_zero[3] = {z_cam[1] * x_zero[2] - z_cam[2] * x_zero[1],
                        z_cam[2] * x_zero[0] - z_cam[0] * x_zero[2],
                        z_cam[0] * x_zero[1] - z_cam[1] * x_zero[0]};

    double R[3][3];
    int i;
    for (i = 0; i < 3; i++) {
        R[0][i] = cos(cam_roll) * x_zero[i] + sin(cam_roll) * y_zero[i];
        R[1][i] = -sin(cam_roll) * x_zero[i] + cos(cam_roll) * y_zero[i];
        R[2][i] = z_cam[i];
    }

    double true_lat_rad = true_lat * M_PI / 180.0;
    double true_lon_rad = true_lon * M_PI / 180.0;
    double true_lst_rad = (gha_aries * M_PI / 180.0) + true_lon_rad;

    double true_zenith_cel[3] = {cos(true_lat_rad) * cos(true_lst_rad),
                                 cos(true_lat_rad) * sin(true_lst_rad), sin(true_lat_rad)};

    starfix_vector3_t zenith_cam;
    zenith_cam.x =
        R[0][0] * true_zenith_cel[0] + R[0][1] * true_zenith_cel[1] + R[0][2] * true_zenith_cel[2];
    zenith_cam.y =
        R[1][0] * true_zenith_cel[0] + R[1][1] * true_zenith_cel[1] + R[1][2] * true_zenith_cel[2];
    zenith_cam.z =
        R[2][0] * true_zenith_cel[0] + R[2][1] * true_zenith_cel[1] + R[2][2] * true_zenith_cel[2];

    int status = (RESET_ARENA(),
                  starfix_solve_position(1, &zenith_cam, (const double (*)[3])R, &gha_aries, ap_lat,
                                         ap_lon, &solved_lat, &solved_lon, &arena, NULL));

    if (status >= 0) {
        double error_lat = solved_lat - true_lat;
        double error_lon = (solved_lon - true_lon) * cos(true_lat_rad);
        double dist_error_nm = sqrt(error_lat * error_lat + error_lon * error_lon) * 60.0;
        if (dist_error_nm < 0.1) {
            printf("Passed\n");
            printf("  Solved Lat = %f, Lon = %f (error = %f nm)\n", solved_lat, solved_lon,
                   dist_error_nm);
        } else {
            printf("FAILED (large error: %f nm)\n", dist_error_nm);
        }
    } else {
        printf("FAILED (solver returned %d)\n", status);
    }

    return 0;
}
