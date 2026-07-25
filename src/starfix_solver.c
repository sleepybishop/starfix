#include "starfix_solver.h"

#include <math.h>
#include <stdlib.h>

#include "../deps/nanoqsp/nanoqsp.h"
#include "starfix_status.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* helper to convert degrees to radians */
static double deg_to_rad(double deg) { return deg * M_PI / 180.0; }

/* helper to convert radians to degrees */
static double rad_to_deg(double rad) { return rad * 180.0 / M_PI; }

starfix_status_t starfix_solve_position(int num_meas, const starfix_vector3_t* zenith_cam,
                                        const double R_est[3][3], const double* gha_aries,
                                        double ap_lat, double ap_lon, double* solved_lat,
                                        double* solved_lon, starfix_arena_t* arena,
                                        starfix_telemetry_t* telem) {
    if (num_meas <= 0) {
        return STARFIX_ERR_INVALID_PARAM;
    }
    if (zenith_cam == NULL || gha_aries == NULL || solved_lat == NULL || solved_lon == NULL) {
        return STARFIX_ERR_NULL_POINTER;
    }
    int i;
    int m = 2 * num_meas;
    for (i = 0; i < num_meas; i++) {
        STARFIX_ASSERT_VALID_DOUBLE(zenith_cam[i].x);
        STARFIX_ASSERT_VALID_DOUBLE(zenith_cam[i].y);
        STARFIX_ASSERT_VALID_DOUBLE(zenith_cam[i].z);
        STARFIX_ASSERT_VALID_DOUBLE(gha_aries[i]);
    }

    int n = 2;

    /* allocate memory for linear system matrix A (m x n) and vector b (m) */
    double* A = starfix_alloc_array(arena, double, (size_t)m*(size_t)n);
    double* b = starfix_alloc_array(arena, double, (size_t)m);

    if (A == NULL || b == NULL) {
        return STARFIX_ERR_MEMORY;
    }

    double ap_lat_rad = deg_to_rad(ap_lat);
    double cos_ap_lat = cos(ap_lat_rad);

    for (i = 0; i < num_meas; i++) {
        /* 1. transform camera zenith measurement to celestial frame: Z_cel = R^T * w_cam */
        const starfix_vector3_t* w = &zenith_cam[i];
        starfix_vector3_t Z_cel;
        Z_cel.x = R_est[0][0] * w->x + R_est[1][0] * w->y + R_est[2][0] * w->z;
        Z_cel.y = R_est[0][1] * w->x + R_est[1][1] * w->y + R_est[2][1] * w->z;
        Z_cel.z = R_est[0][2] * w->x + R_est[1][2] * w->y + R_est[2][2] * w->z;

        /* 2. calculate Local Sidereal Time (LST) at assumed position (AP) */
        double gha_rad = deg_to_rad(gha_aries[i]);
        double ap_lon_rad = deg_to_rad(ap_lon);
        double lst_rad = gha_rad + ap_lon_rad;

        /* 3. compute local East and North unit vectors in celestial coordinates at AP */
        double E_x = -sin(lst_rad);
        double E_y = cos(lst_rad);
        double E_z = 0.0;

        double N_x = -sin(ap_lat_rad) * cos(lst_rad);
        double N_y = -sin(ap_lat_rad) * sin(lst_rad);
        double N_z = cos(ap_lat_rad);

        /* 4. project celestial zenith onto North and East axes to get angular offsets */
        double x_N = Z_cel.x * N_x + Z_cel.y * N_y + Z_cel.z * N_z;
        double x_E = Z_cel.x * E_x + Z_cel.y * E_y + Z_cel.z * E_z;

        /* 5. build linear system rows */
        /* row 2i: delta_lat = x_N */
        A[(2 * i) * 2 + 0] = 1.0;
        A[(2 * i) * 2 + 1] = 0.0;
        b[2 * i] = x_N;

        /* row 2i+1: cos(ap_lat) * delta_lon = x_E */
        A[(2 * i + 1) * 2 + 0] = 0.0;
        A[(2 * i + 1) * 2 + 1] = cos_ap_lat;
        b[2 * i + 1] = x_E;
    }

    /* 6. solve the box-constrained least squares problem using nanoqsp */
    /* state offset variables: x_off = [delta_lat_rad, delta_lon_rad] */
    double x_off[2] = {0.0, 0.0};
    double lb[2] = {-deg_to_rad(5.0), -deg_to_rad(5.0)}; /* +/- 5.0 degrees search bounds */
    double ub[2] = {deg_to_rad(5.0), deg_to_rad(5.0)};

    int nanoqsp_ws_size = 1024; /* more than enough for n=2 */
    double* nanoqsp_ws = starfix_alloc_array(arena, double, nanoqsp_ws_size);
    if (nanoqsp_ws == NULL) return STARFIX_ERR_MEMORY;
    NanoqspConfig config = {.strategy = NANOQSP_STRATEGY_COORDINATE_DESCENT,
                            .max_iterations = 1000,
                            .tolerance = 1e-6,
                            .workspace = nanoqsp_ws,
                            .workspace_size = nanoqsp_ws_size};

    int result = nanoqsp_solve_least_squares(m, n, A, b, lb, ub, x_off, &config);

    if (result >= 0) {
        /* output solved coordinates in degrees */
        *solved_lat = ap_lat + rad_to_deg(x_off[0]);
        *solved_lon = ap_lon + rad_to_deg(x_off[1]);
        if (telem) {
            telem->solver_iterations = (uint32_t)result;
        }
        return STARFIX_SUCCESS;
    }
    return STARFIX_ERR_SOLVER_DIVERGENCE;
}
