#include "starfix_ekf.h"

#include <math.h>
#include <stddef.h>

#include "starfix_status.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static double deg_to_rad(double deg) { return deg * M_PI / 180.0; }

/* preserve both endpoints: -180 stays -180 and +180 stays +180 */
static double wrap_longitude(double lon) {
    while (lon > 180.0) lon -= 360.0;
    while (lon < -180.0) lon += 360.0;
    return lon;
}

starfix_status_t starfix_ekf_init(starfix_ekf_t* filter, double init_lat, double init_lon) {
    if (filter == NULL) {
        return STARFIX_ERR_NULL_POINTER;
    }
    if (!isfinite(init_lat) || !isfinite(init_lon) || init_lat <= -90.0 || init_lat >= 90.0 ||
        init_lon < -180.0 || init_lon > 180.0) {
        return STARFIX_ERR_INVALID_PARAM;
    }
    double noise_deg = 3.5 / 60.0;
    /* initial uncertainty: 0.5 degrees position, 0.15 speed scale, 10 degrees heading */
    *filter = (starfix_ekf_t){
        .lat = init_lat,
        .lon = init_lon,
        .scale = 1.0,
        .P = {{0.5 * 0.5, 0, 0, 0},
              {0, 0.5 * 0.5, 0, 0},
              {0, 0, 0.15 * 0.15, 0},
              {0, 0, 0, deg_to_rad(10.0) * deg_to_rad(10.0)}},
        .Q = {1e-12, 1e-12, 1e-10, 3e-8},
        .R = {noise_deg * noise_deg, noise_deg * noise_deg},
    };
    return STARFIX_SUCCESS;
}

starfix_status_t starfix_ekf_predict(starfix_ekf_t* filter, double speed_knots, double heading_deg,
                                     double dt_sec) {
    if (filter == NULL) {
        return STARFIX_ERR_NULL_POINTER;
    }
    if (!isfinite(speed_knots) || !isfinite(heading_deg) || !isfinite(dt_sec) || dt_sec <= 0.0 ||
        !isfinite(filter->scale) || fabs(filter->scale) < 1e-9 || !isfinite(filter->lat) ||
        filter->lat <= -90.0 || filter->lat >= 90.0) {
        return STARFIX_ERR_INVALID_PARAM;
    }
    double corrected_speed = speed_knots / filter->scale;
    double corrected_heading = deg_to_rad(heading_deg) - filter->bias;

    double dist_nm = corrected_speed * dt_sec / 3600.0;
    double lat_rad = deg_to_rad(filter->lat);

    double pred_d_lat = dist_nm * cos(corrected_heading) / 60.0;
    double pred_d_lon = dist_nm * sin(corrected_heading) / (60.0 * cos(lat_rad));

    filter->lat += pred_d_lat;
    filter->lon = wrap_longitude(filter->lon + pred_d_lon);

    double F[4][4] = {
        {1.0, 0.0, 0.0, 0.0}, {0.0, 1.0, 0.0, 0.0}, {0.0, 0.0, 1.0, 0.0}, {0.0, 0.0, 0.0, 1.0}};

    F[0][2] = -pred_d_lat / filter->scale;
    F[0][3] = dist_nm * sin(corrected_heading) / 60.0;

    F[1][0] = pred_d_lon * tan(lat_rad) * (M_PI / 180.0);
    F[1][2] = -pred_d_lon / filter->scale;
    F[1][3] = -dist_nm * cos(corrected_heading) / (60.0 * cos(lat_rad));

    /* propagate covariance: P_next = F * P * F^T + Q */
    double FP[4][4];

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            FP[i][j] = 0.0;
            for (int k = 0; k < 4; k++) {
                FP[i][j] += F[i][k] * filter->P[k][j];
            }
        }
    }

    /* compute P = FP * F^T + Q */
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            double val = 0.0;
            for (int k = 0; k < 4; k++) {
                val += FP[i][k] * F[j][k];
            }
            filter->P[i][j] = val;
        }
        filter->P[i][i] += filter->Q[i];
    }
    return STARFIX_SUCCESS;
}

starfix_status_t starfix_ekf_correct(starfix_ekf_t* filter, double celestial_lat,
                                     double celestial_lon) {
    if (filter == NULL) {
        return STARFIX_ERR_NULL_POINTER;
    }
    if (!isfinite(celestial_lat) || !isfinite(celestial_lon) || celestial_lat < -90.0 ||
        celestial_lat > 90.0 || celestial_lon < -180.0 || celestial_lon > 180.0) {
        return STARFIX_ERR_INVALID_PARAM;
    }
    double S[2][2];
    double S_inv[2][2];
    double K[4][2];
    double P_next[4][4];

    /* H selects lat/lon, so S is the top-left 2x2 of P plus R */
    S[0][0] = filter->P[0][0] + filter->R[0];
    S[0][1] = filter->P[0][1];
    S[1][0] = filter->P[1][0];
    S[1][1] = filter->P[1][1] + filter->R[1];

    double det = S[0][0] * S[1][1] - S[0][1] * S[1][0];
    if (fabs(det) < 1e-12) {
        return STARFIX_ERR_SOLVER_SINGULAR; /* singular matrix, skip update */
    }

    S_inv[0][0] = S[1][1] / det;
    S_inv[0][1] = -S[0][1] / det;
    S_inv[1][0] = -S[1][0] / det;
    S_inv[1][1] = S[0][0] / det;

    /* K = P H^T S^-1, using the first two columns of P */
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 2; j++) {
            K[i][j] = filter->P[i][0] * S_inv[0][j] + filter->P[i][1] * S_inv[1][j];
        }
    }

    double y[2];
    y[0] = celestial_lat - filter->lat;
    y[1] = wrap_longitude(celestial_lon - filter->lon);

    double inno_norm = y[0] * y[0] + y[1] * y[1];
    if (inno_norm > 10.0) return STARFIX_ERR_EKF_DIVERGENCE;

    /* Gate the measurement before mutating the filter so a rejected update is atomic. */
    filter->lat += K[0][0] * y[0] + K[0][1] * y[1];
    filter->lon += K[1][0] * y[0] + K[1][1] * y[1];
    filter->scale += K[2][0] * y[0] + K[2][1] * y[1];
    filter->bias += K[3][0] * y[0] + K[3][1] * y[1];
    filter->lon = wrap_longitude(filter->lon);

    /* H selects lat/lon; form each coefficient of I - K H as needed */
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            P_next[i][j] = 0.0;
            for (int k = 0; k < 4; k++) {
                double coefficient = (i == k ? 1.0 : 0.0) - (k < 2 ? K[i][k] : 0.0);
                P_next[i][j] += coefficient * filter->P[k][j];
            }
        }
    }

    /* copy back to P and symmetrize */
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            filter->P[i][j] = (P_next[i][j] + P_next[j][i]) / 2.0;
        }
    }
    return STARFIX_SUCCESS;
}
