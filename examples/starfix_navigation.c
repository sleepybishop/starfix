#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "starfix_arena.h"

static uint8_t mempool[10 * 1024 * 1024];
static starfix_arena_t arena = {mempool, mempool + sizeof(mempool)};

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "starfix_ekf.h"
#include "starfix_graph.h"

static int do_fuse_ekf() {
    printf("Starting EKF Sensor Fusion Simulation...\n");
    starfix_ekf_t ekf;
    starfix_ekf_init(&ekf, 34.0522, -118.2437);

    double true_speed = 12.0;
    double true_heading = 45.0;
    double compass_bias = 4.0;
    double speed_scale = 0.93;
    double dt = 10.0;

    double cur_true_lat = 34.0522;
    double cur_true_lon = -118.2437;

    int step;
    srand(42);
    for (step = 0; step < 4320; step++) {
        double current_heading = ((step * dt) < 10800) ? true_heading : 135.0;

        /* update true location */
        double dist_nm = true_speed * dt / 3600.0;
        cur_true_lat += dist_nm * cos(current_heading * M_PI / 180.0) / 60.0;
        cur_true_lon += dist_nm * sin(current_heading * M_PI / 180.0) /
                        (60.0 * cos(cur_true_lat * M_PI / 180.0));

        /* simulate noisy measurements */
        double m_speed = true_speed * speed_scale + (rand() / (double)RAND_MAX - 0.5) * 0.6;
        double m_heading = current_heading + compass_bias + (rand() / (double)RAND_MAX - 0.5) * 2.0;

        starfix_ekf_predict(&ekf, m_speed, m_heading, dt);

        /* celestial fixes every 20 minutes (120 steps) */
        if ((step + 1) % 120 == 0) {
            double n_lat = (rand() / (double)RAND_MAX - 0.5) * 0.08; /* ~3.5 nm error */
            double n_lon = (rand() / (double)RAND_MAX - 0.5) * 0.08;
            starfix_telemetry_t telem = {0};
            starfix_ekf_correct(&ekf, cur_true_lat + n_lat, cur_true_lon + n_lon, &telem);
        }
    }

    double err_lat = ekf.lat - cur_true_lat;
    double err_lon = (ekf.lon - cur_true_lon) * cos(cur_true_lat * M_PI / 180.0);
    double dist_err_nm = sqrt(err_lat * err_lat + err_lon * err_lon) * 60.0;

    printf("Simulation Finished!\n");
    printf("=================================================================\n");
    printf("Total distance traveled: 144.0 nautical miles\n");
    printf("-----------------------------------------------------------------\n");
    printf("EKF Sensor Fusion Error: %.4f nm\n", dist_err_nm);
    printf("=================================================================\n\n");
    printf("Online Calibration of Sensor Biases:\n");
    printf("-----------------------------------------------------------------\n");
    printf("Parameter            | Estimated Value      | True Value     \n");
    printf("-----------------------------------------------------------------\n");
    printf("Compass Bias (deg)   |           %10.4f |          4.0000\n", ekf.bias * 180.0 / M_PI);
    printf("Speed Scale Factor   |           %10.4f |          0.9300\n", ekf.scale);
    printf("-----------------------------------------------------------------\n");

    return 0;
}

static int do_fuse_graph() {
    printf("Generating Sailboat Path & Sensor Readings...\n");
    int num_nodes = 50;
    int steps_per_node = 50;
    double dt_sec = 10.0;
    double init_lat = 34.0522;
    double init_lon = -118.2437;
    int total_steps = num_nodes * steps_per_node;

    starfix_odometry_t* odom =
        (starfix_odometry_t*)malloc((size_t)total_steps * sizeof(starfix_odometry_t));
    starfix_pose_t* poses =
        (starfix_pose_t*)malloc((size_t)(num_nodes + 1) * sizeof(starfix_pose_t));

    if (odom == NULL || poses == NULL) {
        printf("Memory allocation failed.\n");
        return -1;
    }

    double true_speed = 12.0;
    double true_scale = 0.93;

    srand(42);
    int i;
    for (i = 0; i < total_steps; i++) {
        double current_heading = (i < total_steps / 2) ? 45.0 : 135.0;
        odom[i].speed = true_speed * true_scale + (rand() / (double)RAND_MAX - 0.5) * 0.1;
        odom[i].heading = current_heading + 4.0 + (rand() / (double)RAND_MAX - 0.5) * 0.2;
    }

    double* true_lats = (double*)malloc((size_t)(num_nodes + 1) * sizeof(double));
    double* true_lons = (double*)malloc((size_t)(num_nodes + 1) * sizeof(double));
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
        fixes[i].lat = true_lats[(i + 1) * 2] + (rand() / (double)RAND_MAX - 0.5) * 0.001;
        fixes[i].lon = true_lons[(i + 1) * 2] + (rand() / (double)RAND_MAX - 0.5) * 0.001;
    }

    double est_scale = 1.0;
    double est_bias = 0.0;
    starfix_telemetry_t telem = {0};

    starfix_status_t status =
        starfix_solve_graph(num_nodes, steps_per_node, dt_sec, init_lat, init_lon, odom, num_fixes,
                            fixes, poses, &est_scale, &est_bias, arena, &telem);

    if (status == STARFIX_SUCCESS) {
        double err_lat = poses[num_nodes].lat - true_lats[num_nodes];
        double err_lon = (poses[num_nodes].lon - true_lons[num_nodes]) *
                         cos(true_lats[num_nodes] * M_PI / 180.0);
        double dist_err_nm = sqrt(err_lat * err_lat + err_lon * err_lon) * 60.0;

        printf("Factor Graph Optimization Finished!\n");
        printf("=================================================================\n");
        printf("Total distance traveled: 166.7 nautical miles\n");
        printf("-----------------------------------------------------------------\n");
        printf("Factor Graph Smoother Error: %.4f nm\n", dist_err_nm);
        printf("=================================================================\n\n");
        printf("Calibrated Sensor Biases:\n");
        printf("-----------------------------------------------------------------\n");
        printf("Parameter                 | Estimated Value      | True Value     \n");
        printf("-----------------------------------------------------------------\n");
        printf("Compass Bias (deg)        |           %10.4f |          4.0000\n",
               est_bias * 180.0 / M_PI);
        printf("Speed Scale Factor        |           %10.4f |          0.9300\n", est_scale);
        printf("-----------------------------------------------------------------\n");
    } else {
        printf("Factor Graph Optimization failed with code %d\n", status);
    }

    free(odom);
    free(poses);
    free(fixes);
    free(true_lats);
    free(true_lons);
    return 0;
}

int main(int argc, char** argv) {
    if (argc == 2 && strcmp(argv[1], "--fuse-ekf") == 0) return do_fuse_ekf();
    if (argc == 2 && strcmp(argv[1], "--fuse-graph") == 0) return do_fuse_graph();
    fprintf(stderr, "Usage: %s --fuse-ekf | --fuse-graph\n", argv[0]);
    return 1;
}
