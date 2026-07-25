#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "starfix_arena.h"
#include "starfix_status.h"

static uint8_t mempool[10 * 1024 * 1024];
static starfix_arena_t arena;
#define RESET_ARENA() (arena.beg = mempool, arena.end = mempool + sizeof(mempool))

#include "starfix_graph.h"

/* test runner to verify the C implementation of the Factor Graph solver */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main() {
    printf("--- Running Factor Graph Comprehensive Tests ---\n");

    int num_nodes = 50;
    int steps_per_node = 50;
    double dt_sec = 10.0;
    double init_lat = 34.0522;
    double init_lon = -118.2437;

    /* 1. NULL Pointer & Sizing Validations */
    printf("Test 1: Sizing and NULL validations... ");
    double scale = 1.0, bias = 0.0;
    int ret1 =
        (RESET_ARENA(), starfix_solve_graph(0, steps_per_node, dt_sec, init_lat, init_lon, NULL, 0,
                                            NULL, NULL, &scale, &bias, &arena, NULL));
    int ret2 =
        (RESET_ARENA(), starfix_solve_graph(num_nodes, steps_per_node, dt_sec, init_lat, init_lon,
                                            NULL, 0, NULL, NULL, &scale, &bias, &arena, NULL));

    if (ret1 < 0 && ret2 < 0) {
        printf("Passed\n");
    } else {
        printf("FAILED (got %d, %d)\n", ret1, ret2);
    }

    /* 2. Happy Path Factor Graph Optimization */
    printf("Test 2: Happy path factor graph optimization... ");

    int total_steps = num_nodes * steps_per_node;
    starfix_odometry_t* odom =
        (starfix_odometry_t*)malloc((size_t)total_steps * sizeof(starfix_odometry_t));
    starfix_pose_t* poses =
        (starfix_pose_t*)malloc(((size_t)num_nodes + 1) * sizeof(starfix_pose_t));

    if (odom == NULL || poses == NULL) {
        printf("Memory allocation failed.\n");
        return -1;
    }

    /* true parameters */
    double true_speed = 12.0; /* knots */
    double true_scale = 0.93;
    double true_bias = 4.0 * M_PI / 180.0; /* 4.0 degrees in radians */

    /* generate raw measurements with biases and tacking every 500 steps (1.4 hours) */
    int i;
    for (i = 0; i < total_steps; i++) {
        double current_heading = (i < total_steps / 2) ? 45.0 : 135.0;
        odom[i].speed = true_speed * true_scale;
        odom[i].heading = current_heading + 4.0; /* 4 degrees compass bias */
    }

    /* simulate ground truth trajectory to generate odometry and celestial fixes */
    double* true_lats = (double*)malloc(((size_t)num_nodes + 1) * sizeof(double));
    double* true_lons = (double*)malloc(((size_t)num_nodes + 1) * sizeof(double));
    true_lats[0] = init_lat;
    true_lons[0] = init_lon;

    double cur_lat = init_lat;
    double cur_lon = init_lon;
    for (i = 1; i <= num_nodes; i++) {
        int j;
        for (j = 0; j < steps_per_node; j++) {
            int s_idx = (i - 1) * steps_per_node + j;
            double current_heading = (s_idx < total_steps / 2) ? 45.0 : 135.0;
            double dist_nm = true_speed * dt_sec / 3600.0;
            double lat_rad = cur_lat * M_PI / 180.0;
            cur_lat += dist_nm * cos(current_heading * M_PI / 180.0) / 60.0;
            cur_lon += dist_nm * sin(current_heading * M_PI / 180.0) / (60.0 * cos(lat_rad));
        }
        true_lats[i] = cur_lat;
        true_lons[i] = cur_lon;
    }

    int num_fixes = 25;
    starfix_fix_t* fixes = (starfix_fix_t*)malloc((size_t)num_fixes * sizeof(starfix_fix_t));
    for (i = 0; i < num_fixes; i++) {
        fixes[i].node_idx = (i + 1) * 2;
        fixes[i].lat = true_lats[(i + 1) * 2];
        fixes[i].lon = true_lons[(i + 1) * 2];
    }

    /* solve graph */
    double est_scale = 1.0;
    double est_bias = 0.0;
    int status = (RESET_ARENA(), starfix_solve_graph(num_nodes, steps_per_node, dt_sec, init_lat,
                                                     init_lon, odom, num_fixes, fixes, poses,
                                                     &est_scale, &est_bias, &arena, NULL));

    if (status >= 0) {
        /* compute final optimized pose error vs true pose at end node */
        double err_lat = poses[num_nodes].lat - true_lats[num_nodes];
        double err_lon = (poses[num_nodes].lon - true_lons[num_nodes]) *
                         cos(true_lats[num_nodes] * M_PI / 180.0);
        double dist_err_nm = sqrt(err_lat * err_lat + err_lon * err_lon) * 60.0;

        /* verify scale and bias calibration */
        double err_scale = fabs(est_scale - true_scale);
        double err_bias = fabs(est_bias - true_bias) * 180.0 / M_PI;

        if (dist_err_nm < 1.5 && err_scale < 0.05 && err_bias < 0.5) {
            printf("Passed\n");
            printf("  Optimized Final Pose Error: %f nm\n", dist_err_nm);
            printf("  Estimated Speed Scale:     %f (True: %f)\n", est_scale, true_scale);
            printf("  Estimated Heading Bias:    %f deg (True: 4.0 deg)\n",
                   est_bias * 180.0 / M_PI);
        } else {
            printf("FAILED (large error: dist_err=%f nm, scale_err=%f, bias_err=%f deg)\n",
                   dist_err_nm, err_scale, err_bias);
        }
    } else {
        printf("FAILED (solver returned %d)\n", status);
    }

    free(odom);
    free(poses);
    free(fixes);
    free(true_lats);
    free(true_lons);
    return 0;
}
