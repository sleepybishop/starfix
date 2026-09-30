#include <math.h>
#include <stdio.h>

#include "starfix_solver.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* A nontrivial proper rotation: celestial (x,y,z) -> camera (z,x,y). */
static const double R[3][3] = {{0, 0, 1}, {1, 0, 0}, {0, 1, 0}};
static int failures = 0;

static starfix_vector3_t observe(double lat, double lon, double gha, double scale) {
    double dec = lat * M_PI / 180.0;
    double ra = (lon + gha) * M_PI / 180.0;
    return (starfix_vector3_t){scale * sin(dec), scale * cos(dec) * cos(ra),
                               scale * cos(dec) * sin(ra)};
}

static void report(const char* name, int passed) {
    printf("%s... %s\n", name, passed ? "Passed" : "FAILED");
    if (!passed) failures++;
}

static int check_fix(double lat, double lon, double gha) {
    starfix_vector3_t w = observe(lat, lon, gha, 3.0);
    double actual_lat, actual_lon;
    starfix_telemetry_t telem = {.solver_iterations = 999};
    int status = starfix_solve_position(1, &w, R, &gha, &actual_lat, &actual_lon, &telem);
    return status == STARFIX_SUCCESS && fabs(actual_lat - lat) < 1e-9 &&
           fabs(remainder(actual_lon - lon, 360.0)) < 1e-9 && actual_lon >= -180.0 &&
           actual_lon < 180.0 && telem.solver_iterations == 0;
}

int main(void) {
    double lat = 12.0, lon = 34.0, gha = 120.0;
    starfix_vector3_t valid = observe(34.0522, -118.2437, gha, 1.0);
    starfix_vector3_t zero = {0, 0, 0};
    starfix_vector3_t invalid = {NAN, 0, 1};
    double invalid_gha = INFINITY;
    int valid_errors = starfix_solve_position(0, &valid, R, &gha, &lat, &lon, NULL) < 0 &&
                       starfix_solve_position(1, NULL, R, &gha, &lat, &lon, NULL) < 0 &&
                       starfix_solve_position(1, &valid, NULL, &gha, &lat, &lon, NULL) < 0 &&
                       starfix_solve_position(1, &valid, R, NULL, &lat, &lon, NULL) < 0 &&
                       starfix_solve_position(1, &valid, R, &gha, NULL, &lon, NULL) < 0 &&
                       starfix_solve_position(1, &zero, R, &gha, &lat, &lon, NULL) < 0 &&
                       starfix_solve_position(1, &invalid, R, &gha, &lat, &lon, NULL) < 0 &&
                       starfix_solve_position(1, &valid, R, &invalid_gha, &lat, &lon, NULL) < 0;
    report("Test 1: Argument validation",
           valid_errors && fabs(lat - 12.0) < 1e-12 && fabs(lon - 34.0) < 1e-12);
    report("Test 2: Happy path position fix solver", check_fix(34.0522, -118.2437, 120));

    int global_ok = 1;
    const double latitudes[] = {-89.999, -63.78, -20, 0, 34.05, 63.78402778, 89.999};
    const double longitudes[] = {-180, -179.999, -118.2437, -18.01597222, 0, 90, 179.999, 180};
    const double times[] = {-720.5, 0, 120, 357.826258, 1080.25};
    for (unsigned int i = 0; i < sizeof(latitudes) / sizeof(latitudes[0]); i++)
        for (unsigned int j = 0; j < sizeof(longitudes) / sizeof(longitudes[0]); j++)
            for (unsigned int k = 0; k < sizeof(times) / sizeof(times[0]); k++)
                global_ok &= check_fix(latitudes[i], longitudes[j], times[k]);
    report("Test 3: Global fixes and longitude wrapping", global_ok);

    double times_multi[3] = {359, 0, 90};
    starfix_vector3_t measurements[3];
    for (int i = 0; i < 3; i++)
        measurements[i] = observe(-35, 179.95, times_multi[i], (double)i + 1);
    int status = starfix_solve_position(3, measurements, R, times_multi, &lat, &lon, NULL);
    report("Test 4: Average time-aligned Earth-fixed directions",
           status == 0 && fabs(lat + 35) < 1e-9 && fabs(lon - 179.95) < 1e-9);
    double pair_times[2] = {0, 0};
    starfix_vector3_t pair[2] = {observe(0, 179.9, 0, 1), observe(0, -179.9, 0, 1)};
    status = starfix_solve_position(2, pair, R, pair_times, &lat, &lon, NULL);
    report("Test 5: Average across antimeridian",
           status == 0 && fabs(lat) < 1e-9 && fabs(fabs(lon) - 180) < 1e-9);

    int poles_ok = 1;
    for (int sign = -1; sign <= 1; sign += 2) {
        starfix_vector3_t pole = {(double)sign, 0, 0};
        status = starfix_solve_position(1, &pole, R, &gha, &lat, &lon, NULL);
        poles_ok &= status == 0 && fabs(lat - 90 * sign) < 1e-9 && fabs(lon) < 1e-9;
    }
    report("Test 6: Pole longitude convention", poles_ok);
    pair[0] = (starfix_vector3_t){1, 0, 0};
    pair[1] = (starfix_vector3_t){-1, 0, 0};
    report("Test 7: Reject cancelling directions",
           starfix_solve_position(2, pair, R, pair_times, &lat, &lon, NULL) ==
               STARFIX_ERR_SOLVER_SINGULAR);
    return failures ? 1 : 0;
}
