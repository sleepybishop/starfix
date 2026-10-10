#ifndef STARFIX_SOLVER_H
#define STARFIX_SOLVER_H

#include "starfix_geometry.h"
#include "starfix_status.h"

#ifdef __cplusplus
extern "C" {
#endif

/* average zenith measurements at a stationary position using a shared catalog-to-camera R.
 * GHA and outputs are degrees; nonzero zenith vectors are normalized.
 * longitude is in [-180, 180), or zero at a pole; outputs are unchanged on error */
starfix_status_t starfix_solve_position(int num_meas, const starfix_vector3_t* zenith_cam,
                                        const double R_est[3][3], const double* gha_aries,
                                        double* solved_lat, double* solved_lon);

#ifdef __cplusplus
}
#endif
#endif /* STARFIX_SOLVER_H */
