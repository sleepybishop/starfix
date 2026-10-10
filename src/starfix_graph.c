#include "starfix_graph.h"

#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "../deps/nanoqsp/nanoqsp.h"
#include "starfix_status.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static double deg_to_rad(double deg) { return deg * M_PI / 180.0; }

typedef struct {
    int count;
    int cols[6];
    double values[6];
    double residual, weight;
} factor_t;

typedef struct {
    int row, col;
    double value;
} normal_term_t;

static int compare_normal_terms(const void* lhs, const void* rhs) {
    const normal_term_t* a = lhs;
    const normal_term_t* b = rhs;
    if (a->row < b->row) return -1;
    if (a->row > b->row) return 1;
    return (a->col > b->col) - (a->col < b->col);
}

static factor_t scalar_factor(int col, double residual, double sigma) {
    return (factor_t){1, {col}, {1.0}, residual, 1.0 / sigma};
}

static void initialize_trajectory(int nodes, int steps, double dt, const starfix_odometry_t* odom,
                                  double lat, double lon, double scale, double bias, double* x) {
    x[0] = lat;
    x[1] = lon;
    for (int i = 1; i <= nodes; i++) {
        for (int j = 0; j < steps; j++) {
            const starfix_odometry_t* o = &odom[(i - 1) * steps + j];
            double distance = o->speed * dt / (3600.0 * scale);
            double heading = deg_to_rad(o->heading) - bias;
            double dlat = distance * cos(heading) / 60.0;
            double dlon = distance * sin(heading) / (60.0 * cos(deg_to_rad(lat)));
            lat += dlat;
            lon += dlon;
        }
        x[2 * i] = lat;
        x[2 * i + 1] = lon;
    }
    x[2 * nodes + 2] = scale;
    x[2 * nodes + 3] = bias;
}

static void linearize(int nodes, int steps, double dt, const starfix_odometry_t* odom,
                      int num_fixes, const starfix_fix_t* fixes, double lat, double lon,
                      const double* x, factor_t* factors) {
    int s = 2 * nodes + 2, b = s + 1;
    *factors++ = scalar_factor(0, x[0] - lat, 0.5);
    *factors++ = scalar_factor(1, x[1] - lon, 0.5);
    *factors++ = scalar_factor(s, x[s] - 1.0, 0.15);
    *factors++ = scalar_factor(b, x[b], deg_to_rad(10.0));

    double scale = x[s];
    if (scale < 1e-6) scale = 1e-6;
    for (int i = 1; i <= nodes; i++) {
        int prev = 2 * (i - 1), curr = 2 * i;
        double north = 0.0, east = 0.0;
        for (int j = 0; j < steps; j++) {
            const starfix_odometry_t* o = &odom[(i - 1) * steps + j];
            double heading = deg_to_rad(o->heading) - x[b];
            double distance = o->speed * dt / 3600.0;
            north += distance * cos(heading);
            east += distance * sin(heading);
        }

        double mid = deg_to_rad((x[prev] + x[curr]) / 2.0);
        double cos_mid = cos(mid);
        if (fabs(cos_mid) < 1e-6) cos_mid = 1e-6;
        /* scale is measured speed / true speed, matching the EKF */
        double dlat = north / (60.0 * scale);
        double dlon = east / (60.0 * cos_mid * scale);
        double dlon_dlat = -dlon * tan(mid) * 0.5 * (M_PI / 180.0);
        double weight = 1.0 / (0.2 / 60.0); /* 0.2 arcmin odometry sigma */
        *factors++ = (factor_t){4,
                                {prev, curr, s, b},
                                {-1.0, 1.0, dlat / scale, -east / (60.0 * scale)},
                                x[curr] - x[prev] - dlat,
                                weight};
        *factors++ = (factor_t){
            6,
            {prev, prev + 1, curr, curr + 1, s, b},
            {dlon_dlat, -1.0, dlon_dlat, 1.0, dlon / scale, north / (60.0 * cos_mid * scale)},
            x[curr + 1] - x[prev + 1] - dlon,
            weight};
    }

    for (int i = 0; i < num_fixes; i++) {
        int col = 2 * fixes[i].node_idx;
        *factors++ = scalar_factor(col, x[col] - fixes[i].lat, 3.5 / 60.0);
        *factors++ = scalar_factor(col + 1, x[col + 1] - fixes[i].lon, 3.5 / 60.0);
    }
}

/* accumulate J^T J and -J^T r in factor order; columns in each factor are sorted */
static ptrdiff_t normal_equations(int count, const factor_t* factors, normal_term_t* terms,
                                  double* rhs) {
    ptrdiff_t used = 0;
    for (int k = 0; k < count; k++) {
        const factor_t* f = &factors[k];
        double values[6];
        double residual = -(f->residual * f->weight);
        for (int i = 0; i < f->count; i++) values[i] = f->values[i] * f->weight;
        for (int i = 0; i < f->count; i++) {
            if (!(fabs(values[i]) > 0.0)) continue;
            rhs[f->cols[i]] += values[i] * residual;
            for (int j = 0; j < f->count; j++) {
                if (!(fabs(values[j]) > 0.0)) continue;
                terms[used++] = (normal_term_t){f->cols[i], f->cols[j], values[i] * values[j]};
            }
        }
    }
    return used;
}

static int pack_normal_matrix(int n, ptrdiff_t count, normal_term_t* terms, double* values,
                              int* cols, int* rows) {
    qsort(terms, (size_t)count, sizeof(*terms), compare_normal_terms);
    double norm = 0.0;
    ptrdiff_t at = 0;
    int nnz = 0;
    for (int row = 0; row < n; row++) {
        rows[row] = nnz;
        double row_sum = 0.0;
        while (at < count && terms[at].row == row) {
            int col = terms[at].col;
            double sum = 0.0;
            do {
                sum += terms[at++].value;
            } while (at < count && terms[at].row == row && terms[at].col == col);
            if (fabs(sum) > 1e-15) {
                values[nnz] = sum;
                cols[nnz++] = col;
                row_sum += fabs(sum);
            }
        }
        if (row_sum > norm) norm = row_sum;
    }
    rows[n] = nnz;

    /* tikhonov regularization requires a diagonal entry for every state */
    double epsilon = norm > 1e-12 ? norm * 1e-8 : 1e-8;
    for (int row = 0; row < n; row++) {
        int j = rows[row], end = rows[row + 1];
        while (j < end && cols[j] != row) j++;
        if (j == end) return STARFIX_ERR_SOLVER_SINGULAR;
        values[j] += epsilon;
    }
    return nnz;
}

starfix_status_t starfix_solve_graph(int num_nodes, int steps_per_node, double dt_sec,
                                     double initial_lat, double initial_lon,
                                     const starfix_odometry_t* odom, int num_fixes,
                                     const starfix_fix_t* fixes, starfix_pose_t* poses,
                                     double* est_scale, double* est_bias, starfix_arena_t arena) {
    if (!odom || !poses || !est_scale || !est_bias || (num_fixes > 0 && !fixes))
        return STARFIX_ERR_NULL_POINTER;
    if (num_nodes <= 0 || steps_per_node <= 0 || num_fixes < 0 || !isfinite(dt_sec) ||
        dt_sec <= 0.0 || !isfinite(initial_lat) || !isfinite(initial_lon)) {
        return STARFIX_ERR_INVALID_PARAM;
    }
    if (initial_lat <= -90.0 || initial_lat >= 90.0 || initial_lon < -180.0 ||
        initial_lon > 180.0 || !isfinite(*est_scale) || *est_scale < 1e-6 || !isfinite(*est_bias) ||
        num_nodes > (INT_MAX - 4) / 2 || num_fixes > (INT_MAX - 4 - 2 * num_nodes) / 2 ||
        num_nodes > INT_MAX / steps_per_node) {
        return STARFIX_ERR_INVALID_PARAM;
    }

    int total_steps = num_nodes * steps_per_node;
    for (int i = 0; i < total_steps; i++) {
        if (!isfinite(odom[i].speed) || !isfinite(odom[i].heading)) {
            return STARFIX_ERR_DATA_CORRUPTION;
        }
    }
    for (int i = 0; i < num_fixes; i++) {
        if (fixes[i].node_idx < 0 || fixes[i].node_idx > num_nodes || !isfinite(fixes[i].lat) ||
            !isfinite(fixes[i].lon)) {
            return STARFIX_ERR_INVALID_PARAM;
        }
    }

    int n = 2 * (num_nodes + 1) + 2;
    int m = 4 + 2 * num_nodes + 2 * num_fixes;
    /* four priors, 4x4 + 6x6 terms per interval, and two terms per fix */
    uint64_t capacity = 4ULL + 52U * (uint64_t)num_nodes + 2U * (uint64_t)num_fixes;
    if (capacity > (uint64_t)PTRDIFF_MAX) return STARFIX_ERR_INVALID_PARAM;
    ptrdiff_t max_terms = (ptrdiff_t)capacity;

    double* x = starfix_alloc_array(&arena, double, n);
    factor_t* factors = starfix_alloc_array(&arena, factor_t, m);
    double* delta = starfix_alloc_array(&arena, double, n);
    double* rhs = starfix_alloc_array(&arena, double, n);
    normal_term_t* terms = starfix_alloc_array(&arena, normal_term_t, max_terms);
    double* values = starfix_alloc_array(&arena, double, max_terms);
    int* cols = starfix_alloc_array(&arena, int, max_terms);
    int* rows = starfix_alloc_array(&arena, int, n + 1);
    NanoqspCSR matrix = {.n = n, .values = values, .col_indices = cols, .row_ptr = rows};
    double* lower = starfix_alloc_array(&arena, double, n);
    double* upper = starfix_alloc_array(&arena, double, n);
    NanoqspConfig config = {.strategy = NANOQSP_STRATEGY_SPECTRAL_GRADIENT,
                            .max_iterations = 2000,
                            .tolerance = 1e-9,
                            .workspace = starfix_alloc_array(&arena, double, 65536),
                            .workspace_size = 65536};
    if (!x || !factors || !delta || !rhs || !terms || !values || !cols || !rows || !lower ||
        !upper || !config.workspace)
        return STARFIX_ERR_MEMORY;
    for (int i = 0; i < n; i++) {
        lower[i] = -0.2;
        upper[i] = 0.2;
    }

    initialize_trajectory(num_nodes, steps_per_node, dt_sec, odom, initial_lat, initial_lon,
                          *est_scale, *est_bias, x);
    int iterations = 0;
    while (iterations < 10) {
        linearize(num_nodes, steps_per_node, dt_sec, odom, num_fixes, fixes, initial_lat,
                  initial_lon, x, factors);
        memset(rhs, 0, (size_t)n * sizeof(*rhs));
        ptrdiff_t count = normal_equations(m, factors, terms, rhs);
        matrix.nnz = pack_normal_matrix(n, count, terms, values, cols, rows);
        if (matrix.nnz < 0) return matrix.nnz;

        memset(delta, 0, (size_t)n * sizeof(*delta));
        int status = nanoqsp_solve_box_sparse(n, &matrix, rhs, lower, upper, delta, &config);
        if (status < 0 && status != NANOQSP_ERR_MAX_ITERATIONS)
            return STARFIX_ERR_SOLVER_DIVERGENCE;
        double step = 0.0;
        for (int i = 0; i < n; i++) {
            x[i] += delta[i];
            step += delta[i] * delta[i];
        }
        iterations++;
        if (sqrt(step) < 1e-5) break;
    }

    for (int i = 0; i <= num_nodes; i++) poses[i] = (starfix_pose_t){x[2 * i], x[2 * i + 1]};
    *est_scale = x[n - 2];
    *est_bias = x[n - 1];
    return STARFIX_SUCCESS;
}
