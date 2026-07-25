#ifndef STARFIX_SOLVER_H
#define STARFIX_SOLVER_H

#include "starfix_arena.h"
#include "starfix_status.h"

/* starfix_solver: C API wrapper to compute latitude/longitude position fixes
   by fusing camera attitude, Greenwich Hour Angle (GHA) of Aries, and raw
   gravity (Zenith) vectors using the nanoqsp linear least-squares solver. */

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
   - zenith_cat: array of zenith vectors in the catalog frame
   - pos: pointer to latitude/longitude position (in degrees)
   - arena: memory arena for solver allocations
   - returns: number of iterations on success, or a negative value on error */
starfix_status_t starfix_solve_position(int num_meas, const starfix_vector3_t* zenith_cam,
                                        const double R_est[3][3], const double* gha_aries,
                                        double ap_lat, double ap_lon, double* solved_lat,
                                        double* solved_lon, starfix_arena_t* arena,
                                        starfix_telemetry_t* telem);

#endif /* STARFIX_SOLVER_H */
