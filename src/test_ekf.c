#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "starfix_ekf.h"
#include "starfix_status.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void) {
    int failures = 0;
    starfix_ekf_t filter;

    printf("--- Running EKF Comprehensive Tests ---\n");

    /* 1. NULL Pointer Validation */
    printf("Test 1: NULL pointer robustness... ");
    int null_ok = starfix_ekf_init(NULL, 34.0, -118.0) == STARFIX_ERR_NULL_POINTER &&
                  starfix_ekf_predict(NULL, 10.0, 45.0, 10.0) == STARFIX_ERR_NULL_POINTER &&
                  starfix_ekf_correct(NULL, 34.0, -118.0) == STARFIX_ERR_NULL_POINTER;
    if (!null_ok) failures++;
    printf(null_ok ? "Passed\n" : "FAILED (wrong null-pointer status)\n");

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
        failures++;
        printf("FAILED\n");
        printf("  Got Lat = %f, Lon = %f\n", filter.lat, filter.lon);
    }

    /* 3. Celestial Correction */
    printf("Test 3: Celestial correction... ");
    double celestial_lat = 34.1720;
    double celestial_lon = -118.1050;
    starfix_ekf_correct(&filter, celestial_lat, celestial_lon);

    if (fabs(filter.lat - 34.171969) < 1e-4 && fabs(filter.lon - -118.104949) < 1e-4) {
        printf("Passed\n");
        printf("  Corrected: Lat = %f, Lon = %f\n", filter.lat, filter.lon);
        printf("  Estimated: Scale = %f, Bias = %f deg\n", filter.scale,
               filter.bias * 180.0 / M_PI);
    } else {
        failures++;
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
    starfix_ekf_correct(&filter, 40.0, -100.0);
    /* state should not change much because R is huge */
    if (fabs(filter.lat - old_lat) < 1e-6 && fabs(filter.lon - old_lon) < 1e-6) {
        printf("Passed\n");
    } else {
        failures++;
        printf("FAILED (moved too much: %f, %f)\n", filter.lat - old_lat, filter.lon - old_lon);
    }

    /* 5. A rejected innovation must not partially mutate filter state. */
    printf("Test 5: Divergent update is atomic... ");
    starfix_ekf_init(&filter, 10.0, 20.0);
    starfix_ekf_t before = filter;
    int divergent_status = starfix_ekf_correct(&filter, 80.0, 150.0);
    int state_changed = memcmp(&filter, &before, sizeof(filter)) != 0;
    if (divergent_status == STARFIX_ERR_EKF_DIVERGENCE && !state_changed) {
        printf("Passed\n");
    } else {
        failures++;
        printf("FAILED (status=%d, state changed=%d)\n", divergent_status, state_changed);
    }

    printf("Test 6: Longitude boundaries and crossings... ");
    int wrapped = 1;
    for (int sign = -1; sign <= 1; sign += 2) {
        starfix_ekf_init(&filter, 0, sign * 180.0);
        wrapped &= starfix_ekf_predict(&filter, 0, 0, 1) == STARFIX_SUCCESS;
        wrapped &= starfix_ekf_correct(&filter, 0, sign * 180.0) == STARFIX_SUCCESS;
        wrapped &= fabs(filter.lon - sign * 180.0) < 1e-12;
        starfix_ekf_init(&filter, 0, sign * 179.75);
        wrapped &= starfix_ekf_predict(&filter, 60, sign * 90.0, 3600) == STARFIX_SUCCESS;
        wrapped &= fabs(filter.lon + sign * 179.25) < 1e-12;
        starfix_ekf_init(&filter, 0, sign * 179.9);
        wrapped &= starfix_ekf_correct(&filter, 0, -sign * 179.9) == STARFIX_SUCCESS;
        wrapped &= fabs(filter.lon + sign * 179.9) < 0.01;
    }
    if (!wrapped) failures++;
    printf(wrapped ? "Passed\n" : "FAILED (longitude wrap)\n");

    printf("Test 7: Invalid state rejection is atomic... ");
    int invalid_ok = 1;
    for (int field = 0; field < 7; field++) {
        starfix_ekf_init(&filter, 10, 20);
        double* fields[] = {&filter.lat,     &filter.lon,  &filter.scale, &filter.bias,
                            &filter.P[0][1], &filter.Q[0], &filter.R[0]};
        *fields[field] = NAN;
        before = filter;
        invalid_ok &= starfix_ekf_predict(&filter, 10, 45, 1) == STARFIX_ERR_INVALID_PARAM;
        invalid_ok &= memcmp(&filter, &before, sizeof(filter)) == 0;
        invalid_ok &= starfix_ekf_correct(&filter, 10, 20) == STARFIX_ERR_INVALID_PARAM;
        invalid_ok &= memcmp(&filter, &before, sizeof(filter)) == 0;
    }
    if (!invalid_ok) failures++;
    printf(invalid_ok ? "Passed\n" : "FAILED (nonfinite state accepted or modified)\n");

    printf("Test 8: Numerical failures preserve the filter... ");
    int numerical_ok = 1;
    starfix_ekf_init(&filter, 10, 20);
    before = filter;
    numerical_ok &=
        starfix_ekf_predict(&filter, DBL_MAX, 90, DBL_MAX) == STARFIX_ERR_EKF_DIVERGENCE;
    numerical_ok &= memcmp(&filter, &before, sizeof(filter)) == 0;
    numerical_ok &= starfix_ekf_predict(&filter, 6000, 0, 3600) == STARFIX_ERR_EKF_DIVERGENCE;
    numerical_ok &= memcmp(&filter, &before, sizeof(filter)) == 0;
    filter.P[0][0] = DBL_MAX;
    filter.P[1][1] = DBL_MAX;
    before = filter;
    numerical_ok &= starfix_ekf_correct(&filter, 10, 20) == STARFIX_ERR_EKF_DIVERGENCE;
    numerical_ok &= memcmp(&filter, &before, sizeof(filter)) == 0;
    starfix_ekf_init(&filter, 10, 20);
    filter.P[2][0] = filter.P[0][2] = DBL_MAX;
    before = filter;
    numerical_ok &= starfix_ekf_correct(&filter, 11, 20) == STARFIX_ERR_EKF_DIVERGENCE;
    numerical_ok &= memcmp(&filter, &before, sizeof(filter)) == 0;
    if (!numerical_ok) failures++;
    printf(numerical_ok ? "Passed\n" : "FAILED (numerical failure accepted or modified state)\n");

    return failures ? 1 : 0;
}
