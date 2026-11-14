#include "starfix_graph.h"

#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include "starfix_geometry.h"

typedef struct {
    double m[2][2];
} mat2_t;

typedef struct {
    double v[2];
} vec2_t;

/* integrate one interval in nautical miles: speed = measured / scale, heading = measured - bias */
static void interval_displacement(const starfix_odometry_t* odom, int steps, double dt_sec,
                                  double scale, double bias, double* north, double* east) {
    *north = 0.0;
    *east = 0.0;
    for (int j = 0; j < steps; j++) {
        double dist = (odom[j].speed / scale) * dt_sec / 3600.0;
        double heading = starfix_deg_to_rad(odom[j].heading) - bias;
        *north += dist * cos(heading);
        *east += dist * sin(heading);
    }
}

static void initialize_trajectory(int num_nodes, int steps_per_node, double dt_sec,
                                  const starfix_odometry_t* odom, double lat, double lon,
                                  double scale, double bias, double* x) {
    int n = 2 * (num_nodes + 1) + 2;
    x[0] = lat;
    x[1] = lon;
    for (int i = 0; i < num_nodes; i++) {
        double north, east;
        interval_displacement(odom + (size_t)i * (size_t)steps_per_node, steps_per_node, dt_sec,
                              scale, bias, &north, &east);
        double dlat = north / 60.0;
        double mid = starfix_deg_to_rad(lat + 0.5 * dlat);
        double dlon = east / (60.0 * fmax(1e-6, cos(mid)));
        lat += dlat;
        lon = starfix_wrap_deg(lon + dlon);
        x[2 * (i + 1)] = lat;
        x[2 * (i + 1) + 1] = lon;
    }
    x[n - 2] = scale;
    x[n - 1] = bias;
}

static double evaluate_cost(int num_nodes, int steps_per_node, double dt_sec,
                            const starfix_odometry_t* odom, int num_fixes,
                            const starfix_fix_t* fixes, double initial_lat, double initial_lon,
                            const double* x) {
    int n = 2 * (num_nodes + 1) + 2;
    double scale = x[n - 2];
    double bias = x[n - 1];
    if (!(scale >= 1e-4) || !isfinite(scale) || !isfinite(bias)) return HUGE_VAL;

    double cos_init = fmax(1e-6, cos(starfix_deg_to_rad(initial_lat)));
    double r0 = (x[0] - initial_lat) / 10.0;
    double r1 = (starfix_wrap_deg(x[1] - initial_lon) * cos_init) / 10.0;
    double rs = (scale - 1.0) / 0.25;
    double rb = bias / starfix_deg_to_rad(15.0);
    double cost = 0.5 * (r0 * r0 + r1 * r1 + rs * rs + rb * rb);

    double odom_w = 1.0 / (0.2 / 60.0);
    for (int i = 0; i < num_nodes; i++) {
        int prev = 2 * i, curr = 2 * (i + 1);
        double north, east;
        interval_displacement(odom + (size_t)i * (size_t)steps_per_node, steps_per_node, dt_sec,
                              scale, bias, &north, &east);
        double mid = starfix_deg_to_rad(0.5 * (x[prev] + x[curr]));
        double cos_mid = fmax(1e-6, cos(mid));
        double dlat = north / 60.0;
        double dlon = east / (60.0 * cos_mid);
        double res_lat = (x[curr] - x[prev] - dlat) * odom_w;
        double res_lon = (starfix_wrap_deg(x[curr + 1] - x[prev + 1]) - dlon) * cos_mid * odom_w;
        cost += 0.5 * (res_lat * res_lat + res_lon * res_lon);
    }

    double fix_w = 1.0 / (3.5 / 60.0);
    for (int i = 0; i < num_fixes; i++) {
        int col = 2 * fixes[i].node_idx;
        double cos_fix = fmax(1e-6, cos(starfix_deg_to_rad(fixes[i].lat)));
        double rl = (x[col] - fixes[i].lat) * fix_w;
        double rlon = starfix_wrap_deg(x[col + 1] - fixes[i].lon) * cos_fix * fix_w;
        cost += 0.5 * (rl * rl + rlon * rlon);
    }
    return isfinite(cost) ? cost : HUGE_VAL;
}

static void add_interval_row(int k, const double jp0[2], const double jp1[2], const double jg[2],
                             double wr, mat2_t* D, mat2_t* E, mat2_t* B, mat2_t* C, vec2_t* r_p,
                             vec2_t* r_g) {
    for (int a = 0; a < 2; a++) {
        r_p[k].v[a] -= jp0[a] * wr;
        r_p[k + 1].v[a] -= jp1[a] * wr;
        r_g->v[a] -= jg[a] * wr;
        for (int b = 0; b < 2; b++) {
            D[k].m[a][b] += jp0[a] * jp0[b];
            D[k + 1].m[a][b] += jp1[a] * jp1[b];
            E[k].m[a][b] += jp0[a] * jp1[b];
            B[k].m[a][b] += jp0[a] * jg[b];
            B[k + 1].m[a][b] += jp1[a] * jg[b];
            C->m[a][b] += jg[a] * jg[b];
        }
    }
}

static void build_bordered_system(int num_nodes, int steps_per_node, double dt_sec,
                                  const starfix_odometry_t* odom, int num_fixes,
                                  const starfix_fix_t* fixes, double initial_lat,
                                  double initial_lon, const double* x, mat2_t* D, mat2_t* E,
                                  mat2_t* B, mat2_t* C, vec2_t* r_p, vec2_t* r_g) {
    int m_poses = num_nodes + 1;
    int n = 2 * m_poses + 2;
    double scale = x[n - 2];
    double bias = x[n - 1];

    memset(D, 0, (size_t)m_poses * sizeof(*D));
    memset(E, 0, (size_t)num_nodes * sizeof(*E));
    memset(B, 0, (size_t)m_poses * sizeof(*B));
    *C = (mat2_t){0};
    memset(r_p, 0, (size_t)m_poses * sizeof(*r_p));
    *r_g = (vec2_t){0};

    /* initial priors */
    double cos_init = fmax(1e-6, cos(starfix_deg_to_rad(initial_lat)));
    double w_lat0 = 1.0 / 10.0;
    double w_lon0 = cos_init / 10.0;
    double w_scale = 1.0 / 0.25;
    double w_bias = 1.0 / starfix_deg_to_rad(15.0);

    D[0].m[0][0] += w_lat0 * w_lat0;
    r_p[0].v[0] -= (w_lat0 * w_lat0) * (x[0] - initial_lat);
    D[0].m[1][1] += w_lon0 * w_lon0;
    r_p[0].v[1] -= (w_lon0 * w_lon0) * starfix_wrap_deg(x[1] - initial_lon);

    C->m[0][0] += w_scale * w_scale;
    r_g->v[0] -= (w_scale * w_scale) * (scale - 1.0);
    C->m[1][1] += w_bias * w_bias;
    r_g->v[1] -= (w_bias * w_bias) * bias;

    /* odometry chain factors */
    double odom_w = 1.0 / (0.2 / 60.0);
    for (int k = 0; k < num_nodes; k++) {
        int prev = 2 * k, curr = 2 * (k + 1);
        double north, east;
        interval_displacement(odom + (size_t)k * (size_t)steps_per_node, steps_per_node, dt_sec,
                              scale, bias, &north, &east);
        double mid = starfix_deg_to_rad(0.5 * (x[prev] + x[curr]));
        double cos_mid = fmax(1e-6, cos(mid));
        double dlat = north / 60.0;
        double dlon = east / (60.0 * cos_mid);

        /* north-south arc factor */
        double jp0_n[2] = {-odom_w, 0.0};
        double jp1_n[2] = {odom_w, 0.0};
        double jg_n[2] = {(dlat / scale) * odom_w, (-east / (60.0 * scale)) * odom_w};
        double wr_n = (x[curr] - x[prev] - dlat) * odom_w;
        add_interval_row(k, jp0_n, jp1_n, jg_n, wr_n, D, E, B, C, r_p, r_g);

        /* east-west metric arc factor (scaled by cos_mid) */
        double w_lon = odom_w * cos_mid;
        double dlon_dlat = -dlon * tan(mid) * 0.5 * (M_PI / 180.0);
        double jp0_e[2] = {dlon_dlat * w_lon, -w_lon};
        double jp1_e[2] = {dlon_dlat * w_lon, w_lon};
        double jg_e[2] = {(dlon / scale) * w_lon, (north / (60.0 * cos_mid * scale)) * w_lon};
        double wr_e = (starfix_wrap_deg(x[curr + 1] - x[prev + 1]) - dlon) * w_lon;
        add_interval_row(k, jp0_e, jp1_e, jg_e, wr_e, D, E, B, C, r_p, r_g);
    }

    /* celestial fix priors */
    double fix_w = 1.0 / (3.5 / 60.0);
    for (int i = 0; i < num_fixes; i++) {
        int node = fixes[i].node_idx;
        int col = 2 * node;
        double cos_fix = fmax(1e-6, cos(starfix_deg_to_rad(fixes[i].lat)));
        double w_lat = fix_w;
        double w_lon = fix_w * cos_fix;
        D[node].m[0][0] += w_lat * w_lat;
        r_p[node].v[0] -= (w_lat * w_lat) * (x[col] - fixes[i].lat);
        D[node].m[1][1] += w_lon * w_lon;
        r_p[node].v[1] -= (w_lon * w_lon) * starfix_wrap_deg(x[col + 1] - fixes[i].lon);
    }
}

/* 2x2 lower-triangular cholesky factorization: A = L * L^T */
static int chol2(const mat2_t* A, mat2_t* L) {
    double a00 = A->m[0][0];
    if (!(a00 > 1e-18) || !isfinite(a00)) return 0;
    double l00 = sqrt(a00);
    double l10 = A->m[1][0] / l00;
    double d11 = A->m[1][1] - l10 * l10;
    if (!(d11 > 1e-18) || !isfinite(d11)) return 0;
    *L = (mat2_t){{{l00, 0.0}, {l10, sqrt(d11)}}};
    return 1;
}

/* solve L * y = b for 2-vector */
static vec2_t solve_l2_vec(const mat2_t* L, vec2_t b) {
    double y0 = b.v[0] / L->m[0][0];
    double y1 = (b.v[1] - L->m[1][0] * y0) / L->m[1][1];
    return (vec2_t){{y0, y1}};
}

/* solve L^T * x = y for 2-vector */
static vec2_t solve_lt2_vec(const mat2_t* L, vec2_t y) {
    double x1 = y.v[1] / L->m[1][1];
    double x0 = (y.v[0] - L->m[1][0] * x1) / L->m[0][0];
    return (vec2_t){{x0, x1}};
}

/* solve L * X = M for 2x2 matrix */
static mat2_t solve_l2_mat(const mat2_t* L, const mat2_t* M) {
    mat2_t X;
    for (int col = 0; col < 2; col++) {
        vec2_t c = {{M->m[0][col], M->m[1][col]}};
        vec2_t s = solve_l2_vec(L, c);
        X.m[0][col] = s.v[0];
        X.m[1][col] = s.v[1];
    }
    return X;
}

/* direct O(N) block-tridiagonal + 2x2 schur complement solve in place */
static int solve_bordered_system(int num_nodes, double lambda, mat2_t* D, mat2_t* E, mat2_t* B,
                                 mat2_t C, vec2_t* r_p, vec2_t r_g, double* delta) {
    int m_poses = num_nodes + 1;
    int n = 2 * m_poses + 2;

    for (int k = 0; k < m_poses; k++) {
        D[k].m[0][0] += lambda;
        D[k].m[1][1] += lambda;
    }
    C.m[0][0] += lambda;
    C.m[1][1] += lambda;

    /* forward sweep: overwrite D with L, E with U, B with V, r_p with y */
    for (int k = 0; k < m_poses; k++) {
        mat2_t Lk;
        if (!chol2(&D[k], &Lk)) return 0;
        D[k] = Lk;

        if (k > 0) {
            const mat2_t* Uprev = &E[k - 1];
            const mat2_t* Vprev = &B[k - 1];
            vec2_t yprev = r_p[k - 1];
            for (int row = 0; row < 2; row++) {
                r_p[k].v[row] -= Uprev->m[0][row] * yprev.v[0] + Uprev->m[1][row] * yprev.v[1];
                for (int col = 0; col < 2; col++) {
                    B[k].m[row][col] -=
                        Uprev->m[0][row] * Vprev->m[0][col] + Uprev->m[1][row] * Vprev->m[1][col];
                }
            }
        }
        r_p[k] = solve_l2_vec(&Lk, r_p[k]);
        B[k] = solve_l2_mat(&Lk, &B[k]);

        if (k < num_nodes) {
            mat2_t Uk = solve_l2_mat(&Lk, &E[k]);
            E[k] = Uk;
            for (int a = 0; a < 2; a++) {
                for (int b = 0; b < 2; b++) {
                    D[k + 1].m[a][b] -= Uk.m[0][a] * Uk.m[0][b] + Uk.m[1][a] * Uk.m[1][b];
                }
            }
        }
    }

    /* 2x2 schur complement for global (scale, bias) */
    mat2_t S = C;
    vec2_t r_s = r_g;
    for (int k = 0; k < m_poses; k++) {
        const mat2_t* Vk = &B[k];
        vec2_t yk = r_p[k];
        for (int a = 0; a < 2; a++) {
            r_s.v[a] -= Vk->m[0][a] * yk.v[0] + Vk->m[1][a] * yk.v[1];
            for (int b = 0; b < 2; b++) {
                S.m[a][b] -= Vk->m[0][a] * Vk->m[0][b] + Vk->m[1][a] * Vk->m[1][b];
            }
        }
    }
    mat2_t L_s;
    if (!chol2(&S, &L_s)) return 0;
    vec2_t dg = solve_lt2_vec(&L_s, solve_l2_vec(&L_s, r_s));
    delta[n - 2] = dg.v[0];
    delta[n - 1] = dg.v[1];

    /* backward substitution for poses */
    for (int k = m_poses - 1; k >= 0; k--) {
        const mat2_t* Vk = &B[k];
        vec2_t rhs = {{r_p[k].v[0] - (Vk->m[0][0] * dg.v[0] + Vk->m[0][1] * dg.v[1]),
                       r_p[k].v[1] - (Vk->m[1][0] * dg.v[0] + Vk->m[1][1] * dg.v[1])}};
        if (k < num_nodes) {
            const mat2_t* Uk = &E[k];
            double dp0 = delta[2 * (k + 1)];
            double dp1 = delta[2 * (k + 1) + 1];
            rhs.v[0] -= Uk->m[0][0] * dp0 + Uk->m[0][1] * dp1;
            rhs.v[1] -= Uk->m[1][0] * dp0 + Uk->m[1][1] * dp1;
        }
        vec2_t dp = solve_lt2_vec(&D[k], rhs);
        delta[2 * k] = dp.v[0];
        delta[2 * k + 1] = dp.v[1];
    }
    return 1;
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

    int m_poses = num_nodes + 1;
    int n = 2 * m_poses + 2;

    double* x = starfix_alloc_array(&arena, double, n);
    double* x_cand = starfix_alloc_array(&arena, double, n);
    double* delta = starfix_alloc_array(&arena, double, n);
    mat2_t* D = starfix_alloc_array(&arena, mat2_t, m_poses);
    mat2_t* E = starfix_alloc_array(&arena, mat2_t, num_nodes);
    mat2_t* B = starfix_alloc_array(&arena, mat2_t, m_poses);
    vec2_t* r_p = starfix_alloc_array(&arena, vec2_t, m_poses);
    if (!x || !x_cand || !delta || !D || !E || !B || !r_p) return STARFIX_ERR_MEMORY;

    initialize_trajectory(num_nodes, steps_per_node, dt_sec, odom, initial_lat, initial_lon,
                          *est_scale, *est_bias, x);
    double cost = evaluate_cost(num_nodes, steps_per_node, dt_sec, odom, num_fixes, fixes,
                                initial_lat, initial_lon, x);
    if (!isfinite(cost)) return STARFIX_ERR_SOLVER_DIVERGENCE;

    double lambda = 1e-6;
    for (int iter = 0; iter < 20; iter++) {
        int accepted = 0;
        double step_norm_sq = 0.0;
        for (int trial = 0; trial < 8; trial++) {
            mat2_t C;
            vec2_t r_g;
            build_bordered_system(num_nodes, steps_per_node, dt_sec, odom, num_fixes, fixes,
                                  initial_lat, initial_lon, x, D, E, B, &C, r_p, &r_g);
            if (!solve_bordered_system(num_nodes, lambda, D, E, B, C, r_p, r_g, delta)) {
                lambda *= 10.0;
                continue;
            }

            step_norm_sq = 0.0;
            for (int i = 0; i < n; i++) {
                double d = delta[i];
                double lim = (i < n - 2) ? 5.0 : 0.25;
                if (d > lim) d = lim;
                if (d < -lim) d = -lim;
                x_cand[i] = x[i] + d;
                step_norm_sq += d * d;
            }
            if (x_cand[n - 2] < 0.1) x_cand[n - 2] = 0.1;

            if (sqrt(step_norm_sq) < 1e-9) {
                memcpy(x, x_cand, (size_t)n * sizeof(*x));
                accepted = 1;
                break;
            }

            double next_cost = evaluate_cost(num_nodes, steps_per_node, dt_sec, odom, num_fixes,
                                             fixes, initial_lat, initial_lon, x_cand);
            if (next_cost <= cost + 1e-12 * (1.0 + cost)) {
                memcpy(x, x_cand, (size_t)n * sizeof(*x));
                cost = next_cost;
                lambda = fmax(1e-9, lambda * 0.25);
                accepted = 1;
                break;
            }
            lambda *= 8.0;
        }
        if (!accepted) break;
        if (sqrt(step_norm_sq) < 1e-6) break;
    }

    for (int i = 0; i <= num_nodes; i++) {
        poses[i] = (starfix_pose_t){x[2 * i], starfix_wrap_deg(x[2 * i + 1])};
    }
    *est_scale = x[n - 2];
    *est_bias = x[n - 1];
    return STARFIX_SUCCESS;
}
