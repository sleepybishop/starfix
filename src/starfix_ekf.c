#include "starfix_ekf.h"

#include <math.h>
#include <stddef.h>

#include "starfix_status.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* helper to convert degrees to radians */
static double deg_to_rad(double deg) { return deg * M_PI / 180.0; }

starfix_status_t starfix_ekf_init(starfix_ekf_t* filter, double init_lat, double init_lon) {
    if (filter == NULL) {
        return STARFIX_ERR_NULL_POINTER;
    }
    int i, j;

    filter->lat = init_lat;
    filter->lon = init_lon;
    filter->scale = 1.0;
    filter->bias = 0.0;

    /* initialize covariance matrix P */
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            filter->P[i][j] = 0.0;
        }
    }

    /* initial state uncertainties: lat/lon ~0.5 deg, scale ~0.15, bias ~10 deg */
    filter->P[0][0] = 0.5 * 0.5;
    filter->P[1][1] = 0.5 * 0.5;
    filter->P[2][2] = 0.15 * 0.15;
    filter->P[3][3] = deg_to_rad(10.0) * deg_to_rad(10.0);

    /* process noise Q */
    filter->Q[0] = 1e-12; /* lat drift variance per step */
    filter->Q[1] = 1e-12; /* lon drift variance per step */
    filter->Q[2] = 1e-10; /* scale factor drift variance */
    filter->Q[3] = 3e-8;  /* heading bias drift variance */

    /* measurement noise R (3.5 nm default fix error) */
    double noise_deg = 3.5 / 60.0;
    filter->R[0] = noise_deg * noise_deg;
    filter->R[1] = noise_deg * noise_deg;
    return STARFIX_SUCCESS;
}

starfix_status_t starfix_ekf_predict(starfix_ekf_t* filter, double speed_knots, double heading_deg,
                                     double dt_sec) {
    if (filter == NULL) {
        return STARFIX_ERR_NULL_POINTER;
    }
    double corrected_speed = speed_knots / filter->scale;
    double corrected_heading = deg_to_rad(heading_deg) - filter->bias;

    double dist_nm = corrected_speed * dt_sec / 3600.0;
    double lat_rad = deg_to_rad(filter->lat);

    /* compute state prediction updates */
    double pred_d_lat = dist_nm * cos(corrected_heading) / 60.0;
    double pred_d_lon = dist_nm * sin(corrected_heading) / (60.0 * cos(lat_rad));

    filter->lat += pred_d_lat;
    filter->lon += pred_d_lon;

    /* compute Jacobian matrix F */
    double F[4][4] = {
        {1.0, 0.0, 0.0, 0.0}, {0.0, 1.0, 0.0, 0.0}, {0.0, 0.0, 1.0, 0.0}, {0.0, 0.0, 0.0, 1.0}};

    F[0][2] = -pred_d_lat / filter->scale;
    F[0][3] = dist_nm * sin(corrected_heading) / 60.0;

    F[1][0] = pred_d_lon * tan(lat_rad) * (M_PI / 180.0);
    F[1][2] = -pred_d_lon / filter->scale;
    F[1][3] = -dist_nm * cos(corrected_heading) / (60.0 * cos(lat_rad));

    /* propagate covariance: P_next = F * P * F^T + Q */
    double FP[4][4];
    int i, j, k;

    /* compute FP = F * P */
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            FP[i][j] = 0.0;
            for (k = 0; k < 4; k++) {
                FP[i][j] += F[i][k] * filter->P[k][j];
            }
        }
    }

    /* compute P = FP * F^T + Q */
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            double val = 0.0;
            for (k = 0; k < 4; k++) {
                val += FP[i][k] * F[j][k];
            }
            filter->P[i][j] = val;
        }
        /* add process noise Q to diagonal */
        filter->P[i][i] += filter->Q[i];
    }
    return STARFIX_SUCCESS;
}

starfix_status_t starfix_ekf_correct(starfix_ekf_t* filter, double celestial_lat,
                                     double celestial_lon, starfix_telemetry_t* telem) {
    if (filter == NULL) {
        return STARFIX_ERR_NULL_POINTER;
    }
    double S[2][2];
    double S_inv[2][2];
    double K[4][2];
    double KH[4][4];
    double I_KH[4][4];
    double P_next[4][4];
    int i, j, k;

    /* compute innovation covariance S = H * P * H^T + R */
    /* since H = [ [1, 0, 0, 0], [0, 1, 0, 0] ], S is just top-left 2x2 of P + R */
    S[0][0] = filter->P[0][0] + filter->R[0];
    S[0][1] = filter->P[0][1];
    S[1][0] = filter->P[1][0];
    S[1][1] = filter->P[1][1] + filter->R[1];

    /* invert 2x2 matrix S */
    double det = S[0][0] * S[1][1] - S[0][1] * S[1][0];
    if (fabs(det) < 1e-12) {
        return STARFIX_ERR_SOLVER_SINGULAR; /* singular matrix, skip update */
    }

    S_inv[0][0] = S[1][1] / det;
    S_inv[0][1] = -S[0][1] / det;
    S_inv[1][0] = -S[1][0] / det;
    S_inv[1][1] = S[0][0] / det;

    /* compute Kalman Gain K = P * H^T * S_inv */
    /* P * H^T is the 4x2 matrix representing first two columns of P */
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 2; j++) {
            K[i][j] = filter->P[i][0] * S_inv[0][j] + filter->P[i][1] * S_inv[1][j];
        }
    }

    /* compute innovation y = z - H * x */
    double y[2];
    y[0] = celestial_lat - filter->lat;
    y[1] = celestial_lon - filter->lon;
    while (y[1] > 180.0) y[1] -= 360.0;
    while (y[1] < -180.0) y[1] += 360.0;

    /* update state: x = x + K * y */
    filter->lat += K[0][0] * y[0] + K[0][1] * y[1];
    filter->lon += K[1][0] * y[0] + K[1][1] * y[1];
    filter->scale += K[2][0] * y[0] + K[2][1] * y[1];
    filter->bias += K[3][0] * y[0] + K[3][1] * y[1];

    double inno_norm = y[0] * y[0] + y[1] * y[1];
    if (telem) telem->ekf_innovation_norm = inno_norm;
    if (inno_norm > 10.0) return STARFIX_ERR_EKF_DIVERGENCE;

    /* compute covariance update: P = (I - K * H) * P */
    /* K * H is 4x4 matrix with K's columns in first two columns, rest zeros */
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            KH[i][j] = (j < 2) ? K[i][j] : 0.0;
        }
    }

    /* compute I - KH */
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            I_KH[i][j] = (i == j ? 1.0 : 0.0) - KH[i][j];
        }
    }

    /* compute P_next = (I - KH) * P */
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            P_next[i][j] = 0.0;
            for (k = 0; k < 4; k++) {
                P_next[i][j] += I_KH[i][k] * filter->P[k][j];
            }
        }
    }

    /* copy back to P and symmetrize */
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            filter->P[i][j] = (P_next[i][j] + P_next[j][i]) / 2.0;
        }
    }
    return STARFIX_SUCCESS;
}
