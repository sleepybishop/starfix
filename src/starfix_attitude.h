#ifndef STARFIX_ATTITUDE_H
#define STARFIX_ATTITUDE_H

#include "starfix_arena.h"
#include "starfix_status.h"
#include "starfix_solver.h"

/* starfix_attitude: C API for attitude determination. solves Wahba's problem
   using Davenport's Q-method (eigenvalue decomposition of 4x4 matrix K)
   to calculate the optimal attitude quaternion and 3x3 rotation matrix. */

typedef struct {
    double q[4]; /* quaternion: [q0, q1, q2, q3] where q0 is scalar */
} starfix_quaternion_t;

/* solves Wahba's problem to find optimal orientation:
   - num_stars: number of paired star vectors (must be >= 2)
   - w: camera frame unit vectors
   - v: celestial catalog frame unit vectors
   - q_est: output optimal quaternion
   - R_est: output optimal 3x3 rotation matrix (projects catalog to camera)
   - returns: 0 on success, negative value on error */
starfix_status_t starfix_solve_attitude(int num_stars, const starfix_vector3_t* w, const starfix_vector3_t* v,
                           starfix_quaternion_t* q_est, double R_est[3][3]);

/* solves Wahba's problem using RANSAC to reject outliers:
   - num_stars: number of paired star vectors
   - w: camera frame unit vectors
   - v: celestial catalog frame unit vectors
   - threshold_deg: maximum angular error for an inlier (in degrees)
   - q_est: output optimal quaternion
   - R_est: output optimal 3x3 rotation matrix
   - returns: number of inliers on success, negative value on error */
starfix_status_t starfix_solve_attitude_ransac(int num_stars, const starfix_vector3_t* w,
                                  const starfix_vector3_t* v, double threshold_deg,
                                  starfix_quaternion_t* q_est, double R_est[3][3],
                                  starfix_arena_t* arena, starfix_telemetry_t* telem);

/* applies atmospheric refraction correction to a set of observed camera vectors:
   - num_stars: number of stars
   - w: apparent camera vectors (modified in-place to true vectors)
   - z_cam: measured gravity Zenith vector in camera frame
   - temp_c: air temperature in Celsius
   - press_hpa: barometric pressure in hPa
   - returns: 0 on success, negative value on error */
int starfix_correct_refraction(int num_stars, starfix_vector3_t* w, const starfix_vector3_t* z_cam,
                               double temp_c, double press_hpa);

#endif /* STARFIX_ATTITUDE_H */
