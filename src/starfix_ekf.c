#include "starfix_ekf.h"

#include <math.h>
#include <stddef.h>

#include "starfix_geometry.h"
#include "starfix_status.h"

static int valid_filter(const starfix_ekf_t* filter) {
    if (!isfinite(filter->lat) || filter->lat <= -90.0 || filter->lat >= 90.0 ||
        !isfinite(filter->lon) || filter->lon < -180.0 || filter->lon > 180.0 ||
        !isfinite(filter->scale) || filter->scale < 1e-9 || !isfinite(filter->bias))
        return 0;
    for (int i = 0; i < 4; i++) {
        if (!isfinite(filter->Q[i]) || filter->Q[i] < 0 || filter->P[i][i] < 0) return 0;
        for (int j = 0; j < 4; j++)
            if (!isfinite(filter->P[i][j])) return 0;
    }
    for (int i = 0; i < 2; i++)
        if (!isfinite(filter->R[i]) || filter->R[i] < 0) return 0;
    return 1;
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
    double cos_lat = fmax(1e-6, cos(starfix_deg_to_rad(init_lat)));
    double bias_sigma = starfix_deg_to_rad(10.0);
    /* initial uncertainty: 0.5 deg great-circle arc, 0.15 speed scale, 10 deg heading */
    *filter = (starfix_ekf_t){
        .lat = init_lat,
        .lon = init_lon,
        .scale = 1.0,
        .P = {{0.5 * 0.5, 0, 0, 0},
              {0, (0.5 * 0.5) / (cos_lat * cos_lat), 0, 0},
              {0, 0, 0.15 * 0.15, 0},
              {0, 0, 0, bias_sigma * bias_sigma}},
        .Q = {1e-12, 1e-12 / (cos_lat * cos_lat), 1e-10, 3e-8},
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
        !valid_filter(filter)) {
        return STARFIX_ERR_INVALID_PARAM;
    }
    starfix_ekf_t next = *filter;
    double corrected_speed = speed_knots / filter->scale;
    double corrected_heading = starfix_deg_to_rad(heading_deg) - filter->bias;

    double dist_nm = corrected_speed * dt_sec / 3600.0;
    double lat_rad = starfix_deg_to_rad(filter->lat);

    double pred_d_lat = dist_nm * cos(corrected_heading) / 60.0;
    double pred_d_lon = dist_nm * sin(corrected_heading) / (60.0 * cos(lat_rad));

    next.lat += pred_d_lat;
    next.lon = starfix_wrap_deg(filter->lon + pred_d_lon);

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

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            double val = 0.0;
            for (int k = 0; k < 4; k++) {
                val += FP[i][k] * F[j][k];
            }
            next.P[i][j] = val;
        }
        next.P[i][i] += filter->Q[i];
    }
    if (!valid_filter(&next)) return STARFIX_ERR_EKF_DIVERGENCE;
    *filter = next;
    return STARFIX_SUCCESS;
}

starfix_status_t starfix_ekf_correct(starfix_ekf_t* filter, double celestial_lat,
                                     double celestial_lon) {
    if (filter == NULL) {
        return STARFIX_ERR_NULL_POINTER;
    }
    if (!isfinite(celestial_lat) || !isfinite(celestial_lon) || celestial_lat < -90.0 ||
        celestial_lat > 90.0 || celestial_lon < -180.0 || celestial_lon > 180.0 ||
        !valid_filter(filter)) {
        return STARFIX_ERR_INVALID_PARAM;
    }
    starfix_ekf_t next = *filter;
    double S[2][2];
    double S_inv[2][2];
    double K[4][2];
    double P_next[4][4];

    /* scale longitude measurement variance by 1 / cos^2(lat) for metric isotropy */
    double cos_lat = fmax(1e-6, cos(starfix_deg_to_rad(filter->lat)));
    double r_lon = filter->R[1] / (cos_lat * cos_lat);
    S[0][0] = filter->P[0][0] + filter->R[0];
    S[0][1] = filter->P[0][1];
    S[1][0] = filter->P[1][0];
    S[1][1] = filter->P[1][1] + r_lon;

    double det = S[0][0] * S[1][1] - S[0][1] * S[1][0];
    if (!isfinite(det)) return STARFIX_ERR_EKF_DIVERGENCE;
    if (det < 1e-12) {
        return STARFIX_ERR_SOLVER_SINGULAR;
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
    y[1] = starfix_wrap_deg(celestial_lon - filter->lon);

    double dy_east = y[1] * cos_lat;
    double inno_norm = y[0] * y[0] + dy_east * dy_east;
    if (inno_norm > 10.0) return STARFIX_ERR_EKF_DIVERGENCE;

    next.lat += K[0][0] * y[0] + K[0][1] * y[1];
    next.lon += K[1][0] * y[0] + K[1][1] * y[1];
    next.scale += K[2][0] * y[0] + K[2][1] * y[1];
    next.bias += K[3][0] * y[0] + K[3][1] * y[1];
    next.lon = starfix_wrap_deg(next.lon);

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

    /* enforce symmetry on P */
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            next.P[i][j] = 0.5 * (P_next[i][j] + P_next[j][i]);
        }
    }
    if (!valid_filter(&next)) return STARFIX_ERR_EKF_DIVERGENCE;
    *filter = next;
    return STARFIX_SUCCESS;
}
