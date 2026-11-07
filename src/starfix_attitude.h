#ifndef STARFIX_ATTITUDE_H
#define STARFIX_ATTITUDE_H

#include "starfix_arena.h"
#include "starfix_geometry.h"
#include "starfix_status.h"

#ifdef __cplusplus
extern "C" {
#endif

/* fit >= 2 paired unit vectors using davenport's Q-method; R_est maps catalog v to camera w.
 * camera axes: +X right, +Y down, +Z forward; returns zero or a negative error.
 * collinear inputs can succeed, but do not determine a unique rotation */
starfix_status_t starfix_solve_attitude(int num_stars, const starfix_vector3_t* w,
                                        const starfix_vector3_t* v, double R_est[3][3]);

/* consensus attitude fit with threshold_deg angular tolerance; vectors and outputs as above.
 * tests every pair for <= 20 observations, otherwise 190 deterministic samples over all inputs.
 * returns inliers (>= 2) or a negative error */
int starfix_solve_attitude_ransac(int num_stars, const starfix_vector3_t* w,
                                  const starfix_vector3_t* v, double threshold_deg,
                                  double R_est[3][3], starfix_arena_t arena);

/* unrefract w in place; z_cam is zenith, temp_c celsius, press_hpa hPa */
starfix_status_t starfix_correct_refraction(int num_stars, starfix_vector3_t* w,
                                            const starfix_vector3_t* z_cam, double temp_c,
                                            double press_hpa);

#ifdef __cplusplus
}
#endif
#endif /* STARFIX_ATTITUDE_H */
