#include "starfix_solver.h"

#include <math.h>
#include <stddef.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

starfix_status_t starfix_solve_position(int num_meas, const starfix_vector3_t* zenith_cam,
                                        const double R_est[3][3], const double* gha_aries,
                                        double* solved_lat, double* solved_lon,
                                        starfix_telemetry_t* telem) {
    if (num_meas <= 0) return STARFIX_ERR_INVALID_PARAM;
    if (zenith_cam == NULL || R_est == NULL || gha_aries == NULL || solved_lat == NULL ||
        solved_lon == NULL)
        return STARFIX_ERR_NULL_POINTER;
    for (int row = 0; row < 3; row++) {
        for (int col = 0; col < 3; col++) {
            STARFIX_ASSERT_VALID_DOUBLE(R_est[row][col]);
        }
    }

    double sum_x = 0.0, sum_y = 0.0, sum_z = 0.0;
    for (int i = 0; i < num_meas; i++) {
        const starfix_vector3_t* w = &zenith_cam[i];
        STARFIX_ASSERT_VALID_DOUBLE(w->x);
        STARFIX_ASSERT_VALID_DOUBLE(w->y);
        STARFIX_ASSERT_VALID_DOUBLE(w->z);
        STARFIX_ASSERT_VALID_DOUBLE(gha_aries[i]);
        double norm = hypot(hypot(w->x, w->y), w->z);
        if (!isfinite(norm) || norm < 1e-12) return STARFIX_ERR_INVALID_PARAM;
        double wx = w->x / norm, wy = w->y / norm, wz = w->z / norm;

        /* Camera -> celestial using R^T, then undo Earth's sidereal rotation.
         * Average directions in the Earth-fixed frame, so timestamps can differ
         * and longitudes straddling the antimeridian do not average to zero. */
        double x = R_est[0][0] * wx + R_est[1][0] * wy + R_est[2][0] * wz;
        double y = R_est[0][1] * wx + R_est[1][1] * wy + R_est[2][1] * wz;
        double z = R_est[0][2] * wx + R_est[1][2] * wy + R_est[2][2] * wz;
        double gha = remainder(gha_aries[i], 360.0) * M_PI / 180.0;
        sum_x += cos(gha) * x + sin(gha) * y;
        sum_y += -sin(gha) * x + cos(gha) * y;
        sum_z += z;
    }
    double horizontal = hypot(sum_x, sum_y);
    double norm = hypot(horizontal, sum_z);
    if (!isfinite(norm) || norm / (double)num_meas < 1e-12) return STARFIX_ERR_SOLVER_SINGULAR;

    double latitude = atan2(sum_z, horizontal) * 180.0 / M_PI;
    /* Longitude is undefined at a pole; return the documented convention of zero. */
    double longitude = horizontal / norm < 1e-12 ? 0.0 : atan2(sum_y, sum_x) * 180.0 / M_PI;
    if (longitude >= 180.0) longitude -= 360.0;
    *solved_lat = latitude;
    *solved_lon = longitude;
    if (telem) telem->solver_iterations = 0;
    return STARFIX_SUCCESS;
}
