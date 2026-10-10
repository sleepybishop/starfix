#ifndef STARFIX_ATTITUDE_H
#define STARFIX_ATTITUDE_H

#include "starfix_arena.h"
#include "starfix_geometry.h"
#include "starfix_status.h"

/* fit >= 2 paired unit vectors using davenport's Q-method; R_est maps catalog v to camera w.
 * camera axes: +X right, +Y down, +Z forward; returns zero or a negative error.
 * collinear inputs can succeed, but do not determine a unique rotation */
starfix_status_t starfix_solve_attitude(int num_stars, const starfix_vector3_t* w,
                                        const starfix_vector3_t* v, double R_est[3][3]);

/* RANSAC attitude fit with threshold_deg angular tolerance; vectors and outputs as above.
 * Returns inliers or a negative error; refits all observations if no pair survives. */
int starfix_solve_attitude_ransac(int num_stars, const starfix_vector3_t* w,
                                  const starfix_vector3_t* v, double threshold_deg,
                                  double R_est[3][3], starfix_arena_t arena);

/* unrefract w in place; z_cam is zenith, temp_c celsius, press_hpa hPa */
int starfix_correct_refraction(int num_stars, starfix_vector3_t* w, const starfix_vector3_t* z_cam,
                               double temp_c, double press_hpa);

/* apply aberration to catalog unit vectors in place using ICRS observer velocity in km/s */
int starfix_correct_aberration(int num_stars, starfix_vector3_t* v, const double v_obs_km_s[3]);

#endif /* STARFIX_ATTITUDE_H */
