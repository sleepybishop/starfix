#include "starfix_attitude.h"
#include "starfix_status.h"

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define MAX_JACOBI_ITERS 50

/* helper to compute eigenvalues and eigenvectors of a 4x4 symmetric matrix
   using the jacobi eigenvalue method.
   - A: input 4x4 symmetric matrix (destroyed during computation)
   - V: output 4x4 matrix of eigenvectors (columns are eigenvectors)
   - d: output 4-element array of eigenvalues */
static int jacobi_4x4(double A[4][4], double V[4][4], double d[4]) {
    int i, j, ip, iq;
    double c, s, t, theta, h, g, thresh;
    double b[4], z[4];

    /* initialize V as identity matrix */
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            V[i][j] = (i == j) ? 1.0 : 0.0;
        }
        b[i] = d[i] = A[i][i];
        z[i] = 0.0;
    }

    /* run jacobi iterations */
    for (i = 1; i <= MAX_JACOBI_ITERS; i++) {
        double sum_off = 0.0;
        for (ip = 0; ip < 3; ip++) {
            for (iq = ip + 1; iq < 4; iq++) {
                sum_off += fabs(A[ip][iq]);
            }
        }

        /* check for convergence */
        if (sum_off < 1e-15) {
            return STARFIX_SUCCESS; /* converged */
        }

        /* threshold for rotation */
        if (i < 4) {
            thresh = 0.2 * sum_off / 16.0;
        } else {
            thresh = 0.0;
        }

        for (ip = 0; ip < 3; ip++) {
            for (iq = ip + 1; iq < 4; iq++) {
                g = 100.0 * fabs(A[ip][iq]);

                /* skip rotation if off-diagonal element is extremely small */
                if (i > 4 && g < 1e-15 * fabs(d[ip]) && g < 1e-15 * fabs(d[iq])) {
                    A[ip][iq] = 0.0;
                } else if (fabs(A[ip][iq]) > thresh) {
                    h = d[iq] - d[ip];
                    if (g < 1e-15 * fabs(h)) {
                        t = A[ip][iq] / h;
                    } else {
                        theta = 0.5 * h / A[ip][iq];
                        t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
                        if (theta < 0.0) {
                            t = -t;
                        }
                    }

                    c = 1.0 / sqrt(1.0 + t * t);
                    s = t * c;
                    h = t * A[ip][iq];

                    z[ip] -= h;
                    z[iq] += h;
                    d[ip] -= h;
                    d[iq] += h;
                    A[ip][iq] = 0.0;

                    /* rotate rows ip and iq */
                    for (j = 0; j < ip; j++) {
                        g = A[j][ip];
                        h = A[j][iq];
                        A[j][ip] = c * g - s * h;
                        A[j][iq] = s * g + c * h;
                    }
                    for (j = ip + 1; j < iq; j++) {
                        g = A[ip][j];
                        h = A[j][iq];
                        A[ip][j] = c * g - s * h;
                        A[j][iq] = s * g + c * h;
                    }
                    for (j = iq + 1; j < 4; j++) {
                        g = A[ip][j];
                        h = A[iq][j];
                        A[ip][j] = c * g - s * h;
                        A[iq][j] = s * g + c * h;
                    }

                    /* accumulate eigenvectors in V */
                    for (j = 0; j < 4; j++) {
                        g = V[j][ip];
                        h = V[j][iq];
                        V[j][ip] = c * g - s * h;
                        V[j][iq] = s * g + c * h;
                    }
                }
            }
        }

        for (ip = 0; ip < 4; ip++) {
            b[ip] += z[ip];
            d[ip] = b[ip];
            z[ip] = 0.0;
        }
    }

    return STARFIX_ERR_MEMORY; /* failed to converge */
}

starfix_status_t starfix_solve_attitude(int num_stars, const starfix_vector3_t* w, const starfix_vector3_t* v,
                           starfix_quaternion_t* q_est, double R_est[3][3]) {
    int i, j, k;

    if (num_stars < 2) {
        return STARFIX_ERR_INVALID_PARAM; /* need at least 2 stars */
    }
    if (w == NULL || v == NULL || q_est == NULL || R_est == NULL) {
        return STARFIX_ERR_ATTITUDE_DEGENERATE;
    }

    /* 1. compute 3x3 covariance matrix B */
    double B[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    for (i = 0; i < num_stars; i++) {
        B[0][0] += w[i].x * v[i].x;
        B[0][1] += w[i].x * v[i].y;
        B[0][2] += w[i].x * v[i].z;

        B[1][0] += w[i].y * v[i].x;
        B[1][1] += w[i].y * v[i].y;
        B[1][2] += w[i].y * v[i].z;

        B[2][0] += w[i].z * v[i].x;
        B[2][1] += w[i].z * v[i].y;
        B[2][2] += w[i].z * v[i].z;
    }

    /* 2. compute matrix trace sigma and symmetric S */
    double sigma = B[0][0] + B[1][1] + B[2][2];
    double S[3][3];
    for (j = 0; j < 3; j++) {
        for (k = 0; k < 3; k++) {
            S[j][k] = B[j][k] + B[k][j];
        }
    }

    /* 3. compute vector Z */
    double Z[3];
    Z[0] = B[2][1] - B[1][2];
    Z[1] = B[0][2] - B[2][0];
    Z[2] = B[1][0] - B[0][1];

    /* 4. construct Davenport's K matrix (4x4 symmetric) */
    double K_mat[4][4];
    K_mat[0][0] = sigma;
    K_mat[0][1] = Z[0];
    K_mat[0][2] = Z[1];
    K_mat[0][3] = Z[2];

    K_mat[1][0] = Z[0];
    K_mat[1][1] = S[0][0] - sigma;
    K_mat[1][2] = S[0][1];
    K_mat[1][3] = S[0][2];

    K_mat[2][0] = Z[1];
    K_mat[2][1] = S[1][0];
    K_mat[2][2] = S[1][1] - sigma;
    K_mat[2][3] = S[1][2];

    K_mat[3][0] = Z[2];
    K_mat[3][1] = S[2][0];
    K_mat[3][2] = S[2][1];
    K_mat[3][3] = S[2][2] - sigma;

    /* 5. find eigenvalues and eigenvectors of K */
    double V_eig[4][4];
    double eigenvalues[4];
    int status = jacobi_4x4(K_mat, V_eig, eigenvalues);

    if (status < 0) {
        return status; /* solver error */
    }

    /* 6. find eigenvector corresponding to largest eigenvalue */
    int max_idx = 0;
    double max_val = eigenvalues[0];
    for (i = 1; i < 4; i++) {
        if (eigenvalues[i] > max_val) {
            max_val = eigenvalues[i];
            max_idx = i;
        }
    }

    /* extract quaternion (eigenvector) */
    double q0 = V_eig[0][max_idx];
    double q1 = V_eig[1][max_idx];
    double q2 = V_eig[2][max_idx];
    double q3 = V_eig[3][max_idx];

    /* normalize quaternion */
    double q_norm = sqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
    if (q_norm > 1e-12) {
        q0 /= q_norm;
        q1 /= q_norm;
        q2 /= q_norm;
        q3 /= q_norm;
    }

    /* make scalar part q0 non-negative by convention */
    if (q0 < 0) {
        q0 = -q0;
        q1 = -q1;
        q2 = -q2;
        q3 = -q3;
    }

    q_est->q[0] = q0;
    q_est->q[1] = q1;
    q_est->q[2] = q2;
    q_est->q[3] = q3;

    /* 7. convert quaternion to 3x3 rotation matrix */
    R_est[0][0] = 1.0 - 2.0 * (q2 * q2 + q3 * q3);
    R_est[0][1] = 2.0 * (q1 * q2 - q0 * q3);
    R_est[0][2] = 2.0 * (q1 * q3 + q0 * q2);

    R_est[1][0] = 2.0 * (q1 * q2 + q0 * q3);
    R_est[1][1] = 1.0 - 2.0 * (q1 * q1 + q3 * q3);
    R_est[1][2] = 2.0 * (q2 * q3 - q0 * q1);

    R_est[2][0] = 2.0 * (q1 * q3 - q0 * q2);
    R_est[2][1] = 2.0 * (q2 * q3 + q0 * q1);
    R_est[2][2] = 1.0 - 2.0 * (q1 * q1 + q2 * q2);

    return STARFIX_SUCCESS;
}

starfix_status_t starfix_solve_attitude_ransac(int num_stars, const starfix_vector3_t* w,
                                  const starfix_vector3_t* v, double threshold_deg,
                                  starfix_quaternion_t* q_est, double R_est[3][3],
                                  starfix_arena_t* arena, starfix_telemetry_t* telem) {
    if (num_stars < 2) {
        return STARFIX_ERR_INVALID_PARAM;
    }
    int _i;
    for (_i = 0; _i < num_stars; _i++) {
        STARFIX_ASSERT_VALID_DOUBLE(w[_i].x);
        STARFIX_ASSERT_VALID_DOUBLE(w[_i].y);
        STARFIX_ASSERT_VALID_DOUBLE(w[_i].z);
    }
    if (w == NULL || v == NULL || q_est == NULL || R_est == NULL) {
        return STARFIX_ERR_ATTITUDE_DEGENERATE;
    }

    /* Convert threshold to cosine value to avoid calling acos() inside the loop */
    double threshold_rad = threshold_deg * M_PI / 180.0;
    double cos_threshold = cos(threshold_rad);

    int max_iterations = 20;
    int best_inlier_count = 0;
    int* best_inliers = starfix_alloc_array(arena, int, (size_t)num_stars);
    int* current_inliers = starfix_alloc_array(arena, int, (size_t)num_stars);

    if (best_inliers == NULL || current_inliers == NULL) {
        return STARFIX_ERR_MEMORY;
    }

    starfix_vector3_t w_sample[2];
    starfix_vector3_t v_sample[2];
    starfix_quaternion_t q_tmp;
    double R_tmp[3][3];

    int iter;
    for (iter = 0; iter < max_iterations; iter++) {
        /* 1. randomly sample 2 distinct indices */
        int idx1 = rand() % num_stars;
        int idx2 = rand() % num_stars;
        while (idx2 == idx1 && num_stars > 1) {
            idx2 = rand() % num_stars;
        }

        w_sample[0] = w[idx1];
        w_sample[1] = w[idx2];
        v_sample[0] = v[idx1];
        v_sample[1] = v[idx2];

        /* 2. fit attitude model to sample */
        int status = starfix_solve_attitude(2, w_sample, v_sample, &q_tmp, R_tmp);
        if (status < 0) {
            continue;
        }

        /* 3. count inliers */
        int current_inlier_count = 0;
        int i;
        for (i = 0; i < num_stars; i++) {
            /* rotate catalog star: v_rot = R_tmp * v[i] */
            double rx = R_tmp[0][0] * v[i].x + R_tmp[0][1] * v[i].y + R_tmp[0][2] * v[i].z;
            double ry = R_tmp[1][0] * v[i].x + R_tmp[1][1] * v[i].y + R_tmp[1][2] * v[i].z;
            double rz = R_tmp[2][0] * v[i].x + R_tmp[2][1] * v[i].y + R_tmp[2][2] * v[i].z;

            /* dot product with observed vector */
            double dot = rx * w[i].x + ry * w[i].y + rz * w[i].z;
            if (dot >= cos_threshold) {
                current_inliers[current_inlier_count++] = i;
            }
        }

        /* 4. update best model if we found more inliers */
        if (current_inlier_count > best_inlier_count) {
            best_inlier_count = current_inlier_count;
            int k_in;
            for (k_in = 0; k_in < best_inlier_count; k_in++) {
                best_inliers[k_in] = current_inliers[k_in];
            }
        }
    }

    /* 5. re-fit using all inliers in the best set */
    if (best_inlier_count >= 2) {
        starfix_vector3_t* w_inliers =
            starfix_alloc_array(arena, starfix_vector3_t, (size_t)best_inlier_count);
        starfix_vector3_t* v_inliers =
            starfix_alloc_array(arena, starfix_vector3_t, (size_t)best_inlier_count);

        if (w_inliers == NULL || v_inliers == NULL) {
            return STARFIX_ERR_MEMORY;
        }
        int i;
        for (i = 0; i < best_inlier_count; i++) {
            int idx = best_inliers[i];
            w_inliers[i] = w[idx];
            v_inliers[i] = v[idx];
        }
        int status = starfix_solve_attitude(best_inlier_count, w_inliers, v_inliers, q_est, R_est);
        if (status < 0) {
            return status;
        }
    } else {
        /* fallback to standard least-squares if RANSAC failed to find inliers */
        int res = starfix_solve_attitude(num_stars, w, v, q_est, R_est);
        if (res < 0) {
            return res;
        }
        if (telem) { telem->ransac_inliers = (uint32_t)num_stars; telem->ransac_inlier_ratio = 1.0; }
    return STARFIX_SUCCESS;
    }

    if (telem) { telem->ransac_inliers = (uint32_t)best_inlier_count; telem->ransac_inlier_ratio = (double)best_inlier_count / num_stars; }
    return STARFIX_SUCCESS;
}

int starfix_correct_refraction(int num_stars, starfix_vector3_t* w, const starfix_vector3_t* z_cam,
                               double temp_c, double press_hpa) {
    if (w == NULL || z_cam == NULL || num_stars <= 0) {
        return STARFIX_ERR_INVALID_PARAM;
    }

    /* Standard refraction constant at 10 deg C and 1010 hPa:
       60.2 arcseconds = 0.00029186 radians */
    double A_std = 60.2 * M_PI / (3600.0 * 180.0);

    /* Scale based on actual temperature and pressure */
    double A = A_std * (press_hpa / 1010.0) * (283.15 / (273.15 + temp_c));

    int i;
    for (i = 0; i < num_stars; i++) {
        /* Dot product to get cosine of apparent Zenith angle */
        double cos_theta = w[i].x * z_cam->x + w[i].y * z_cam->y + w[i].z * z_cam->z;

        /* Clamp to valid range */
        if (cos_theta > 1.0) cos_theta = 1.0;
        if (cos_theta < -1.0) cos_theta = -1.0;

        double theta_a = acos(cos_theta);

        /* Only correct if angle is not completely at Zenith (where refraction is zero)
           and not close to horizon (where tan(theta) blows up) */
        if (theta_a > 1e-6 && theta_a < 1.3) {
            double R = A * tan(theta_a);
            double theta_t = theta_a + R;

            /* Orthogonal component to z_cam in the plane of w */
            double ux = w[i].x - cos_theta * z_cam->x;
            double uy = w[i].y - cos_theta * z_cam->y;
            double uz = w[i].z - cos_theta * z_cam->z;

            double u_norm = sqrt(ux * ux + uy * uy + uz * uz);
            if (u_norm > 1e-12) {
                ux /= u_norm;
                uy /= u_norm;
                uz /= u_norm;

                /* Reconstruct true vector shifted away from Zenith */
                w[i].x = cos(theta_t) * z_cam->x + sin(theta_t) * ux;
                w[i].y = cos(theta_t) * z_cam->y + sin(theta_t) * uy;
                w[i].z = cos(theta_t) * z_cam->z + sin(theta_t) * uz;

                /* Normalize */
                double w_norm = sqrt(w[i].x * w[i].x + w[i].y * w[i].y + w[i].z * w[i].z);
                if (w_norm > 1e-12) {
                    w[i].x /= w_norm;
                    w[i].y /= w_norm;
                    w[i].z /= w_norm;
                }
            }
        }
    }

    return STARFIX_SUCCESS;
}
