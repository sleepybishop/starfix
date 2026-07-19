#include <math.h>
#include "starfix_status.h"
#include <stdio.h>
#include <stdlib.h>

#include "starfix_ekf.h"

/* test runner to verify the C EKF implementation under multiple conditions */

int main() {
    starfix_ekf_t filter;

    printf("--- Running EKF Comprehensive Tests ---\n");

    /* 1. NULL Pointer Validation */
    printf("Test 1: NULL pointer robustness... ");
    starfix_ekf_init(NULL, 34.0, -118.0);
    starfix_ekf_predict(NULL, 10.0, 45.0, 10.0);
    starfix_ekf_correct(NULL, 34.0, -118.0, NULL);
    printf("Passed (no crash)\n");

    /* 2. Happy Path Initialization & DR Propagation */
    printf("Test 2: Happy path propagation... ");
    double init_lat = 34.0522;
    double init_lon = -118.2437;
    starfix_ekf_init(&filter, init_lat, init_lon);

    double speed = 10.0;
    double heading = 45.0;
    double dt = 10.0;
    int steps = 360; /* 1 hour */
    int i;
    for (i = 0; i < steps; i++) {
        starfix_ekf_predict(&filter, speed, heading, dt);
    }

    if (fabs(filter.lat - 34.170051) < 1e-4 && fabs(filter.lon - -118.101360) < 1e-4) {
        printf("Passed\n");
        printf("  Lat = %f, Lon = %f\n", filter.lat, filter.lon);
    } else {
        printf("FAILED\n");
        printf("  Got Lat = %f, Lon = %f\n", filter.lat, filter.lon);
    }

    /* 3. Celestial Correction */
    printf("Test 3: Celestial correction... ");
    double celestial_lat = 34.1720;
    double celestial_lon = -118.1050;
    starfix_ekf_correct(&filter, celestial_lat, celestial_lon, NULL);

    if (fabs(filter.lat - 34.171969) < 1e-4 && fabs(filter.lon - -118.104949) < 1e-4) {
        printf("Passed\n");
        printf("  Corrected: Lat = %f, Lon = %f\n", filter.lat, filter.lon);
        printf("  Estimated: Scale = %f, Bias = %f deg\n", filter.scale,
               filter.bias * 180.0 / M_PI);
    } else {
        printf("FAILED\n");
        printf("  Got Lat = %f, Lon = %f\n", filter.lat, filter.lon);
    }

    /* 4. Singular Covariance / Extreme noise check */
    printf("Test 4: Extreme noise correction limits... ");
    /* set measurement noise extremely large */
    filter.R[0] = 1e12;
    filter.R[1] = 1e12;
    double old_lat = filter.lat;
    double old_lon = filter.lon;
    /* correct with a far-away point */
    starfix_ekf_correct(&filter, 40.0, -100.0, NULL);
    /* state should not change much because R is huge */
    if (fabs(filter.lat - old_lat) < 1e-6 && fabs(filter.lon - old_lon) < 1e-6) {
        printf("Passed\n");
    } else {
        printf("FAILED (moved too much: %f, %f)\n", filter.lat - old_lat, filter.lon - old_lon);
    }

    return 0;
}
