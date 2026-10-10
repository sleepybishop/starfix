#ifndef STARFIX_EKF_H
#define STARFIX_EKF_H

#include "starfix_status.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double lat;     /* latitude in degrees */
    double lon;     /* longitude in degrees */
    double scale;   /* positive speed log scale factor: measured speed / true speed */
    double bias;    /* compass heading bias in radians */
    double P[4][4]; /* state covariance matrix */
    double Q[4];    /* diagonal process noise covariance */
    double R[2];    /* diagonal measurement noise covariance */
} starfix_ekf_t;

/* initialize finite state and covariances; latitude must be strictly between the poles */
starfix_status_t starfix_ekf_init(starfix_ekf_t* filter, double init_lat, double init_lon);

/* propagate by dead reckoning; errors leave the filter unchanged */
starfix_status_t starfix_ekf_predict(starfix_ekf_t* filter, double speed_knots, double heading_deg,
                                     double dt_sec);

/* correct with a celestial fix; errors leave the filter unchanged */
starfix_status_t starfix_ekf_correct(starfix_ekf_t* filter, double celestial_lat,
                                     double celestial_lon);

#ifdef __cplusplus
}
#endif
#endif /* STARFIX_EKF_H */
