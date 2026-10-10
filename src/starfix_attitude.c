#include "starfix_attitude.h"

#include <math.h>
#include <stddef.h>
#include <stdlib.h>

#include "starfix_astrometry.h"
#include "starfix_status.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define MAX_JACOBI_ITERS 50

static void rotate(double* a, double* b, double c, double s) {
    double x = *a, y = *b;
    *a = c * x - s * y;
    *b = s * x + c * y;
}

/* jacobi decomposition: destroys A; V holds eigenvectors as columns */
static int jacobi_4x4(double A[4][4], double V[4][4], double d[4]) {
    double b[4], z[4];

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            V[i][j] = (i == j) ? 1.0 : 0.0;
        }
        b[i] = d[i] = A[i][i];
        z[i] = 0.0;
    }

    for (int i = 1; i <= MAX_JACOBI_ITERS; i++) {
        double sum_off = 0.0;
        for (int ip = 0; ip < 3; ip++) {
            for (int iq = ip + 1; iq < 4; iq++) {
                sum_off += fabs(A[ip][iq]);
            }
        }

        if (sum_off < 1e-15) {
            return STARFIX_SUCCESS;
        }

        double thresh = i < 4 ? 0.2 * sum_off / 16.0 : 0.0;

        for (int ip = 0; ip < 3; ip++) {
            for (int iq = ip + 1; iq < 4; iq++) {
                double g = 100.0 * fabs(A[ip][iq]);

                /* skip negligible off-diagonal elements */
                if (i > 4 && g < 1e-15 * fabs(d[ip]) && g < 1e-15 * fabs(d[iq])) {
                    A[ip][iq] = 0.0;
                    continue;
                }
                if (!(fabs(A[ip][iq]) > thresh)) continue;
                double h = d[iq] - d[ip];
                double t;
                if (g < 1e-15 * fabs(h)) {
                    t = A[ip][iq] / h;
                } else {
                    double theta = 0.5 * h / A[ip][iq];
                    t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
                    if (theta < 0.0) {
                        t = -t;
                    }
                }

                double c = 1.0 / sqrt(1.0 + t * t);
                double s = t * c;
                h = t * A[ip][iq];

                z[ip] -= h;
                z[iq] += h;
                d[ip] -= h;
                d[iq] += h;
                A[ip][iq] = 0.0;

                for (int j = 0; j < ip; j++) {
                    rotate(&A[j][ip], &A[j][iq], c, s);
                }
                for (int j = ip + 1; j < iq; j++) {
                    rotate(&A[ip][j], &A[j][iq], c, s);
                }
                for (int j = iq + 1; j < 4; j++) {
                    rotate(&A[ip][j], &A[iq][j], c, s);
                }

                for (int j = 0; j < 4; j++) {
                    rotate(&V[j][ip], &V[j][iq], c, s);
                }
            }
        }

        for (int ip = 0; ip < 4; ip++) {
            b[ip] += z[ip];
            d[ip] = b[ip];
            z[ip] = 0.0;
        }
    }

    return STARFIX_ERR_SOLVER_DIVERGENCE;
}

starfix_status_t starfix_solve_attitude(int num_stars, const starfix_vector3_t* w,
                                        const starfix_vector3_t* v, double R_est[3][3]) {
    if (num_stars < 2) {
        return STARFIX_ERR_INVALID_PARAM;
    }
    if (!w || !v || !R_est) return STARFIX_ERR_NULL_POINTER;

    double B[3][3] = {{0}};
    for (int i = 0; i < num_stars; i++) {
        double camera[3] = {w[i].x, w[i].y, w[i].z};
        double catalog[3] = {v[i].x, v[i].y, v[i].z};
        for (int axis = 0; axis < 3; axis++)
            if (!isfinite(camera[axis]) || !isfinite(catalog[axis]))
                return STARFIX_ERR_DATA_CORRUPTION;
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++) B[j][k] += camera[j] * catalog[k];
    }

    /* davenport's K for the catalog-to-camera attitude */
    double sigma = B[0][0] + B[1][1] + B[2][2];
    double Z[3] = {B[2][1] - B[1][2], B[0][2] - B[2][0], B[1][0] - B[0][1]};
    double K_mat[4][4];
    K_mat[0][0] = sigma;
    for (int j = 0; j < 3; j++) {
        K_mat[0][j + 1] = K_mat[j + 1][0] = Z[j];
        for (int k = 0; k < 3; k++) K_mat[j + 1][k + 1] = B[j][k] + B[k][j];
        K_mat[j + 1][j + 1] -= sigma;
    }

    double V_eig[4][4];
    double eigenvalues[4];
    int status = jacobi_4x4(K_mat, V_eig, eigenvalues);

    if (status < 0) {
        return status;
    }

    /* the largest eigenvalue gives the least-squares attitude */
    int max_idx = 0;
    for (int i = 1; i < 4; i++) {
        if (eigenvalues[i] > eigenvalues[max_idx]) {
            max_idx = i;
        }
    }

    double q0 = V_eig[0][max_idx];
    double q1 = V_eig[1][max_idx];
    double q2 = V_eig[2][max_idx];
    double q3 = V_eig[3][max_idx];

    double q_norm = sqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
    if (q_norm > 1e-12) {
        q0 /= q_norm;
        q1 /= q_norm;
        q2 /= q_norm;
        q3 /= q_norm;
    }

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

int starfix_solve_attitude_ransac(int num_stars, const starfix_vector3_t* w,
                                  const starfix_vector3_t* v, double threshold_deg,
                                  double R_est[3][3], starfix_arena_t arena) {
    if (num_stars < 2) {
        return STARFIX_ERR_INVALID_PARAM;
    }
    if (w == NULL || v == NULL || R_est == NULL) {
        return STARFIX_ERR_NULL_POINTER;
    }
    if (!isfinite(threshold_deg) || threshold_deg <= 0.0 || threshold_deg >= 180.0) {
        return STARFIX_ERR_INVALID_PARAM;
    }
    for (int i = 0; i < num_stars; i++) {
        if (!isfinite(w[i].x) || !isfinite(w[i].y) || !isfinite(w[i].z) || !isfinite(v[i].x) ||
            !isfinite(v[i].y) || !isfinite(v[i].z))
            return STARFIX_ERR_DATA_CORRUPTION;
    }

    /* Convert threshold to cosine value to avoid calling acos() inside the loop */
    double threshold_rad = threshold_deg * M_PI / 180.0;
    double cos_threshold = cos(threshold_rad);

    int max_iterations = 20;
    int best_inlier_count = 0;
    int* best_inliers = starfix_alloc_array(&arena, int, (size_t)num_stars);
    int* current_inliers = starfix_alloc_array(&arena, int, (size_t)num_stars);

    if (best_inliers == NULL || current_inliers == NULL) {
        return STARFIX_ERR_MEMORY;
    }

    starfix_vector3_t w_sample[2];
    starfix_vector3_t v_sample[2];
    double R_tmp[3][3];

    for (int iter = 0; iter < max_iterations; iter++) {
        int idx1 = rand() % num_stars;
        int idx2 = rand() % num_stars;
        while (idx2 == idx1 && num_stars > 1) {
            idx2 = rand() % num_stars;
        }

        w_sample[0] = w[idx1];
        w_sample[1] = w[idx2];
        v_sample[0] = v[idx1];
        v_sample[1] = v[idx2];

        int status = starfix_solve_attitude(2, w_sample, v_sample, R_tmp);
        if (status < 0) {
            continue;
        }

        int current_inlier_count = 0;
        for (int i = 0; i < num_stars; i++) {
            starfix_vector3_t rotated = starfix_rotate((const double (*)[3])R_tmp, v[i]);
            double dot = starfix_dot(rotated, w[i]);
            if (dot >= cos_threshold) {
                current_inliers[current_inlier_count++] = i;
            }
        }

        if (current_inlier_count > best_inlier_count) {
            best_inlier_count = current_inlier_count;
            int* swap = best_inliers;
            best_inliers = current_inliers;
            current_inliers = swap;
        }
    }

    /* refit the winning set; fall back to all observations if no pair survived */
    int count = num_stars;
    if (best_inlier_count >= 2) {
        starfix_vector3_t* w_inliers =
            starfix_alloc_array(&arena, starfix_vector3_t, best_inlier_count);
        starfix_vector3_t* v_inliers =
            starfix_alloc_array(&arena, starfix_vector3_t, best_inlier_count);
        if (!w_inliers || !v_inliers) return STARFIX_ERR_MEMORY;
        for (int i = 0; i < best_inlier_count; i++) {
            w_inliers[i] = w[best_inliers[i]];
            v_inliers[i] = v[best_inliers[i]];
        }
        w = w_inliers;
        v = v_inliers;
        count = best_inlier_count;
    }
    int status = starfix_solve_attitude(count, w, v, R_est);
    if (status < 0) return status;
    return count;
}

starfix_status_t starfix_correct_refraction(int num_stars, starfix_vector3_t* w,
                                            const starfix_vector3_t* z_cam, double temp_c,
                                            double press_hpa) {
    if (!w || !z_cam) return STARFIX_ERR_NULL_POINTER;
    if (num_stars <= 0) return STARFIX_ERR_INVALID_PARAM;

    for (int i = 0; i < num_stars; i++) {
        double cos_theta = starfix_dot(w[i], *z_cam);

        if (cos_theta > 1.0) cos_theta = 1.0;
        if (cos_theta < -1.0) cos_theta = -1.0;

        double theta_a = acos(cos_theta);
        double alt_a = M_PI / 2.0 - theta_a;
        double alt_t;

        /* refraction vanishes at zenith; the model requires positive altitude */
        if (theta_a > 1e-6 && alt_a > 0.0) {
            starfix_unrefract(press_hpa, temp_c, alt_a, &alt_t);
            double theta_t = M_PI / 2.0 - alt_t;

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

starfix_status_t starfix_correct_aberration(int num_stars, starfix_vector3_t* v,
                                            const double v_obs_km_s[3]) {
    if (!v || !v_obs_km_s) return STARFIX_ERR_NULL_POINTER;
    if (num_stars <= 0) return STARFIX_ERR_INVALID_PARAM;

    const double c_kms = 299792.458;
    double beta_x = v_obs_km_s[0] / c_kms;
    double beta_y = v_obs_km_s[1] / c_kms;
    double beta_z = v_obs_km_s[2] / c_kms;

    for (int i = 0; i < num_stars; i++) {
        double s_dot_b = v[i].x * beta_x + v[i].y * beta_y + v[i].z * beta_z;

        /* First order stellar aberration correction: s' = s + beta - s(s . beta) */
        double sp_x = v[i].x + beta_x - v[i].x * s_dot_b;
        double sp_y = v[i].y + beta_y - v[i].y * s_dot_b;
        double sp_z = v[i].z + beta_z - v[i].z * s_dot_b;

        double norm = sqrt(sp_x * sp_x + sp_y * sp_y + sp_z * sp_z);
        if (norm > 1e-12) {
            v[i].x = sp_x / norm;
            v[i].y = sp_y / norm;
            v[i].z = sp_z / norm;
        }
    }

    return STARFIX_SUCCESS;
}
