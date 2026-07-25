#include "starfix_graph.h"

#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "../deps/nanoqsp/nanoqsp.h"
#include "starfix_status.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* helper to convert degrees to radians */
static double deg_to_rad(double deg) { return deg * M_PI / 180.0; }

starfix_status_t starfix_solve_graph(int num_nodes, int steps_per_node, double dt_sec,
                                     double initial_lat, double initial_lon,
                                     const starfix_odometry_t* odom, int num_fixes,
                                     const starfix_fix_t* fixes, starfix_pose_t* poses,
                                     double* est_scale, double* est_bias, starfix_arena_t* arena,
                                     starfix_telemetry_t* telem) {
    int i, j, k, iter;

    if (num_nodes <= 0 || steps_per_node <= 0 || dt_sec <= 0.0 || odom == NULL || poses == NULL ||
        est_scale == NULL || est_bias == NULL) {
        return STARFIX_ERR_INVALID_PARAM;
    }

    int N = 2 * (num_nodes + 1) + 2;
    int M = 4 + 2 * num_nodes + 2 * num_fixes;

    /* allocate state vector X */
    double* X = starfix_alloc_array(arena, double, (size_t)N);
    /* allocate Jacobian matrix J (M x N) and residuals vector r (M) */
    double* J = starfix_alloc_array(arena, double, (size_t)M*(size_t)N);
    double* r = starfix_alloc_array(arena, double, (size_t)M);
    double* b = starfix_alloc_array(arena, double, (size_t)M);
    double* W = starfix_alloc_array(arena, double, (size_t)M);
    double* delta_X = starfix_alloc_array(arena, double, (size_t)N);

    if (X == NULL || J == NULL || r == NULL || b == NULL || W == NULL || delta_X == NULL) {
        return STARFIX_ERR_NULL_POINTER;
    }

    /* 1. initialize state vector X using raw dead reckoning */
    X[0] = initial_lat;
    X[1] = initial_lon;

    double dr_lat = initial_lat;
    double dr_lon = initial_lon;

    for (i = 1; i <= num_nodes; i++) {
        int start_step = (i - 1) * steps_per_node;
        for (j = 0; j < steps_per_node; j++) {
            int s_idx = start_step + j;
            double dist_nm = odom[s_idx].speed * dt_sec / 3600.0;
            double dr_lat_rad = deg_to_rad(dr_lat);
            double d_lat = dist_nm * cos(deg_to_rad(odom[s_idx].heading)) / 60.0;
            double d_lon =
                dist_nm * sin(deg_to_rad(odom[s_idx].heading)) / (60.0 * cos(dr_lat_rad));
            dr_lat += d_lat;
            dr_lon += d_lon;
        }
        X[2 * i] = dr_lat;
        X[2 * i + 1] = dr_lon;
    }

    X[N - 2] = 1.0; /* speed log scale initial guess */
    X[N - 1] = 0.0; /* compass heading bias initial guess */

    /* standard deviation uncertainties */
    double sigma_prior_pos = 0.5; /* deg */
    double sigma_prior_scale = 0.15;
    double sigma_prior_bias = deg_to_rad(10.0);

    double sigma_odom_lat = 0.2 / 60.0; /* deg per keyframe */
    double sigma_odom_lon = 0.2 / 60.0;

    double sigma_celestial = 3.5 / 60.0; /* deg */

    /* 2. Gauss-Newton optimization loop */
    int max_iters = 10;
    int converged = 0;

    for (iter = 0; iter < max_iters; iter++) {
        /* reset Jacobian and residuals */
        memset(J, 0, (size_t)M * (size_t)N * sizeof(double));

        int r_idx = 0;

        /* 2.1 prior factors */
        r[r_idx] = X[0] - initial_lat;
        W[r_idx] = 1.0 / sigma_prior_pos;
        J[r_idx * N + 0] = 1.0;
        r_idx++;

        r[r_idx] = X[1] - initial_lon;
        W[r_idx] = 1.0 / sigma_prior_pos;
        J[r_idx * N + 1] = 1.0;
        r_idx++;

        r[r_idx] = X[N - 2] - 1.0;
        W[r_idx] = 1.0 / sigma_prior_scale;
        J[r_idx * N + (N - 2)] = 1.0;
        r_idx++;

        r[r_idx] = X[N - 1] - 0.0;
        W[r_idx] = 1.0 / sigma_prior_bias;
        J[r_idx * N + (N - 1)] = 1.0;
        r_idx++;

        /* 2.2 odometry factors */
        double s_val = X[N - 2];
        double b_val = X[N - 1];

        if (fabs(s_val) < 1e-6) s_val = 1e-6;

        for (i = 1; i <= num_nodes; i++) {
            double lat_prev = X[2 * (i - 1)];
            double lon_prev = X[2 * (i - 1) + 1];
            double lat_curr = X[2 * i];
            double lon_curr = X[2 * i + 1];

            double sum_cos = 0.0;
            double sum_sin = 0.0;
            int start_step = (i - 1) * steps_per_node;

            for (j = 0; j < steps_per_node; j++) {
                int s_idx = start_step + j;
                double h_rad = deg_to_rad(odom[s_idx].heading) - b_val;
                double dist_nm = odom[s_idx].speed * dt_sec / 3600.0;
                sum_cos += dist_nm * cos(h_rad);
                sum_sin += dist_nm * sin(h_rad);
            }

            double lat_prev_rad = deg_to_rad(lat_prev);
            double cos_lat_prev = cos(lat_prev_rad);
            if (fabs(cos_lat_prev) < 1e-6) cos_lat_prev = 1e-6;

            double exp_d_lat = sum_cos / (60.0 * s_val);
            double exp_d_lon = sum_sin / (60.0 * s_val * cos_lat_prev);

            /* lat odometry residual */
            r[r_idx] = lat_curr - lat_prev - exp_d_lat;
            W[r_idx] = 1.0 / sigma_odom_lat;
            J[r_idx * N + 2 * (i - 1)] = -1.0;
            J[r_idx * N + 2 * i] = 1.0;
            J[r_idx * N + (N - 2)] = exp_d_lat / s_val;
            J[r_idx * N + (N - 1)] = -sum_sin / (60.0 * s_val);
            r_idx++;

            /* lon odometry residual */
            r[r_idx] = lon_curr - lon_prev - exp_d_lon;
            W[r_idx] = 1.0 / sigma_odom_lon;
            J[r_idx * N + 2 * (i - 1)] = -exp_d_lon * tan(lat_prev_rad);
            J[r_idx * N + 2 * (i - 1) + 1] = -1.0;
            J[r_idx * N + 2 * i + 1] = 1.0;
            J[r_idx * N + (N - 2)] = exp_d_lon / s_val;
            J[r_idx * N + (N - 1)] = sum_cos / (60.0 * s_val * cos_lat_prev);
            r_idx++;
        }

        /* 2.3 celestial fix factors */
        for (i = 0; i < num_fixes; i++) {
            int n_idx = fixes[i].node_idx;
            if (n_idx < 0 || n_idx > num_nodes) continue;

            r[r_idx] = X[2 * n_idx] - fixes[i].lat;
            W[r_idx] = 1.0 / sigma_celestial;
            J[r_idx * N + 2 * n_idx] = 1.0;
            r_idx++;

            r[r_idx] = X[2 * n_idx + 1] - fixes[i].lon;
            W[r_idx] = 1.0 / sigma_celestial;
            J[r_idx * N + 2 * n_idx + 1] = 1.0;
            r_idx++;
        }

        /* 2.4 apply weights to Jacobian and residuals */
        for (i = 0; i < M; i++) {
            r[i] *= W[i];
            b[i] = -r[i];
            for (j = 0; j < N; j++) {
                J[i * N + j] *= W[i];
            }
        }

        /* 2.5 solve update step using nanoqsp sparse CSR solver */
        for (i = 0; i < N; i++) {
            delta_X[i] = 0.0;
        }

        double* D_dense = starfix_alloc_array(arena, double, (size_t)N*(size_t)N);
        double* d_vec = starfix_alloc_array(arena, double, (size_t)N);
        if (D_dense == NULL || d_vec == NULL) {
            return STARFIX_ERR_MEMORY;
        }

        /* Compute D_dense = J^T * J and d_vec = J^T * b */
        for (i = 0; i < N; i++) {
            double d_sum = 0.0;
            for (k = 0; k < M; k++) {
                d_sum += J[k * N + i] * b[k];
            }
            d_vec[i] = d_sum;

            for (j = 0; j < N; j++) {
                double h_sum = 0.0;
                for (k = 0; k < M; k++) {
                    h_sum += J[k * N + i] * J[k * N + j];
                }
                D_dense[i * N + j] = h_sum;
            }
        }

        /* Calculate infinity norm for Tikhonov regularization */
        double norm_inf = 0.0;
        for (i = 0; i < N; i++) {
            double row_sum = 0.0;
            for (j = 0; j < N; j++) {
                row_sum += fabs(D_dense[i * N + j]);
            }
            if (row_sum > norm_inf) {
                norm_inf = row_sum;
            }
        }
        double epsilon = 1e-8;
        if (norm_inf > 1e-12) {
            epsilon = norm_inf * 1e-8;
        }
        for (i = 0; i < N; i++) {
            D_dense[i * N + i] += epsilon;
        }

        /* Count non-zeros */
        int nnz = 0;
        for (i = 0; i < N * N; i++) {
            if (fabs(D_dense[i]) > 1e-15) {
                nnz++;
            }
        }

        /* Allocate CSR buffers */
        double* val = starfix_alloc_array(arena, double, (size_t)nnz);
        int* col_ind = starfix_alloc_array(arena, int, (size_t)nnz);
        int* row_ptr = starfix_alloc_array(arena, int, (size_t)N + 1);
        if (val == NULL || col_ind == NULL || row_ptr == NULL) {
            return STARFIX_ERR_MEMORY;
        }

        /* Populate CSR arrays */
        int idx = 0;
        for (i = 0; i < N; i++) {
            row_ptr[i] = idx;
            for (j = 0; j < N; j++) {
                if (fabs(D_dense[i * N + j]) > 1e-15) {
                    val[idx] = D_dense[i * N + j];
                    col_ind[idx] = j;
                    idx++;
                }
            }
        }
        row_ptr[N] = idx;

        NanoqspCSR D_sparse;
        D_sparse.n = N;
        D_sparse.nnz = nnz;
        D_sparse.values = val;
        D_sparse.col_indices = col_ind;
        D_sparse.row_ptr = row_ptr;

        double step_bound = 0.2; /* limit updates to +/- 0.2 per iteration */
        double* lb = starfix_alloc_array(arena, double, (size_t)N);
        double* ub = starfix_alloc_array(arena, double, (size_t)N);
        if (lb == NULL || ub == NULL) {
            return STARFIX_ERR_MEMORY;
        }
        for (i = 0; i < N; i++) {
            lb[i] = -step_bound;
            ub[i] = step_bound;
        }

        int nanoqsp_ws_size = 65536;
        double* nanoqsp_ws = starfix_alloc_array(arena, double, nanoqsp_ws_size);
        if (nanoqsp_ws == NULL) return STARFIX_ERR_INVALID_PARAM;
        NanoqspConfig config = {.strategy = NANOQSP_STRATEGY_COORDINATE_DESCENT,
                                .max_iterations = 1000,
                                .tolerance = 1e-6,
                                .workspace = nanoqsp_ws,
                                .workspace_size = nanoqsp_ws_size};
        config.strategy = NANOQSP_STRATEGY_SPECTRAL_GRADIENT;
        config.max_iterations = 2000;
        config.tolerance = 1e-9;

        int result = nanoqsp_solve_box_sparse(N, &D_sparse, d_vec, lb, ub, delta_X, &config);

        if (result < 0) {
            /* solver failed */
            return STARFIX_ERR_MEMORY;
        }

        /* 2.6 update state vector */
        double step_norm = 0.0;
        for (i = 0; i < N; i++) {
            X[i] += delta_X[i];
            step_norm += delta_X[i] * delta_X[i];
        }
        step_norm = sqrt(step_norm);

        if (step_norm < 1e-5) {
            converged = 1;
            break;
        }
    }

    /* 3. copy optimized states back to output arrays */
    for (i = 0; i <= num_nodes; i++) {
        poses[i].lat = X[2 * i];
        poses[i].lon = X[2 * i + 1];
    }
    *est_scale = X[N - 2];
    *est_bias = X[N - 1];

    if (telem) {
        telem->solver_iterations = (uint32_t)(converged ? (iter + 1) : max_iters);
    }
    return STARFIX_SUCCESS;
}
