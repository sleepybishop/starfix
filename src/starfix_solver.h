#ifndef STARFIX_SOLVER_H
#define STARFIX_SOLVER_H

#include "starfix_status.h"

/* starfix_solver: C API wrapper to compute latitude/longitude position fixes
   by fusing camera attitude, Greenwich Hour Angle (GHA) of Aries, and raw
   gravity (Zenith) vectors by direct spherical conversion. */

typedef struct {
    double x; /* x component */
    double y; /* y component */
    double z; /* z component */
} starfix_vector3_t;

/* solves for the absolute lat/lon position:
   - num_meas: number of zenith measurements
   - zenith_cam: array of measured zenith vectors in the camera frame
   - R_est: 3x3 camera rotation matrix (projects catalog to camera frame)
   - gha_aries: array of Greenwich Hour Angles of Aries (in degrees)
   - solved_lat/solved_lon: output degrees; longitude is in [-180, 180)
   - observations share R_est and a stationary position; directions are averaged
     after rotating each observation into the Earth-fixed frame using its timestamp
   - no initial position or workspace is needed; nonzero zenith vectors are normalized
   - at an exact pole, longitude is undefined and reported as zero
   - returns: STARFIX_SUCCESS or a negative error; outputs are unchanged on failure
   - telem: optional telemetry; solver_iterations is zero (no iterative search) */
starfix_status_t starfix_solve_position(int num_meas, const starfix_vector3_t* zenith_cam,
                                        const double R_est[3][3], const double* gha_aries,
                                        double* solved_lat, double* solved_lon,
                                        starfix_telemetry_t* telem);

#endif /* STARFIX_SOLVER_H */
