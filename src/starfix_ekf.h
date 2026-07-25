#include "starfix_status.h"

#ifndef STARFIX_EKF_H
#define STARFIX_EKF_H

/* starfix_ekf: lightweight c implementation of extended kalman filter
   for sailboat/vehicle gps-denied navigation. fuses dead reckoning
   (speed log and compass) with celestial navigation updates. */

typedef struct {
    double lat;     /* latitude in degrees */
    double lon;     /* longitude in degrees */
    double scale;   /* speed log scale factor */
    double bias;    /* compass heading bias in radians */
    double P[4][4]; /* state covariance matrix */
    double Q[4];    /* diagonal process noise covariance */
    double R[2];    /* diagonal measurement noise covariance */
} starfix_ekf_t;

/* initialize filter state and covariances */
starfix_status_t starfix_ekf_init(starfix_ekf_t* filter, double init_lat, double init_lon);

/* propagate state using dead reckoning speed and heading */
starfix_status_t starfix_ekf_predict(starfix_ekf_t* filter, double speed_knots, double heading_deg,
                                     double dt_sec);

/* correct state using celestial lat/lon fix */
starfix_status_t starfix_ekf_correct(starfix_ekf_t* filter, double celestial_lat,
                                     double celestial_lon, starfix_telemetry_t* telem);

#endif /* STARFIX_EKF_H */
