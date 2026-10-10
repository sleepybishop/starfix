#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "starfix_attitude.h"
#include "starfix_centroid.h"
#include "starfix_graph.h"
#include "starfix_identify.h"
#include "starfix_nmea.h"

static int failures;

static void check(int passed, const char* name) {
    if (passed) return;
    fprintf(stderr, "FAILED: %s\n", name);
    failures++;
}

int main(void) {
    unsigned char image[4] = {0};
    starfix_centroid_t centroids[5] = {{0}};
    int u[2], v[2];
    check(starfix_find_centroids(NULL, 2, 2, 1, 1, centroids, NULL, 1) == STARFIX_ERR_NULL_POINTER,
          "centroid null image");
    check(
        starfix_find_centroids(image, 0, 2, 1, 1, centroids, NULL, 1) == STARFIX_ERR_INVALID_PARAM,
        "centroid dimensions");
    check(starfix_calibrate_hot_pixels(NULL, 2, 2, 1, 2, u, v) == STARFIX_ERR_NULL_POINTER &&
              starfix_calibrate_hot_pixels(image, 0, 2, 1, 2, u, v) == STARFIX_ERR_INVALID_PARAM &&
              starfix_calibrate_hot_pixels(image, INT_MAX, 2, 1, 2, u, v) ==
                  STARFIX_ERR_INVALID_PARAM,
          "hot-pixel null image and dimensions");
    check(starfix_apply_hot_pixel_mask(image, 2, 2, -1, u, v) == STARFIX_ERR_INVALID_PARAM &&
              starfix_apply_hot_pixel_mask(NULL, 2, 2, 0, u, v) == STARFIX_ERR_NULL_POINTER,
          "hot-pixel count and null image");

    starfix_vector3_t stars[2] = {{1, 0, 0}, {0, 1, 0}};
    double R[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, before[3][3];
    memcpy(before, R, sizeof(R));
    check(starfix_solve_attitude(2, NULL, stars, R) == STARFIX_ERR_NULL_POINTER &&
              starfix_solve_attitude(1, stars, stars, R) == STARFIX_ERR_INVALID_PARAM,
          "attitude pointer and count errors");
    stars[0].x = NAN;
    check(starfix_solve_attitude(2, stars, stars, R) == STARFIX_ERR_DATA_CORRUPTION &&
              memcmp(before, R, sizeof(R)) == 0,
          "nonfinite observations preserve attitude output");
    check(starfix_correct_refraction(1, NULL, stars, 10, 1010) == STARFIX_ERR_NULL_POINTER &&
              starfix_correct_refraction(0, stars, stars, 10, 1010) == STARFIX_ERR_INVALID_PARAM,
          "refraction pointer and count errors");
    double velocity[3] = {0};
    check(starfix_correct_aberration(1, stars, NULL) == STARFIX_ERR_NULL_POINTER &&
              starfix_correct_aberration(0, stars, velocity) == STARFIX_ERR_INVALID_PARAM,
          "aberration pointer and count errors");

    starfix_db_t db = {0};
    check(starfix_load_db("/dev/null/starfix.bin", &db) == STARFIX_ERR_IO && !db.stars &&
              !db.hashes && !db.num_stars && !db.num_hashes && !db.bins,
          "database open error leaves an empty database");
    check(starfix_load_db("/dev/null", &db) == STARFIX_ERR_DATA_CORRUPTION,
          "empty database is malformed data");
    starfix_catalog_star_t catalog[5] = {{0}};
    db.stars = catalog;
    db.num_stars = 5;
    starfix_match_t matches[5];
    starfix_arena_t scratch = {0};
    check(starfix_identify_stars(2, centroids, 10, 10, 10, &db, 5, matches, scratch, NULL) ==
                  STARFIX_ERR_INVALID_PARAM &&
              starfix_identify_stars(5, centroids, 10, 10, 10, &db, 2, matches, scratch, NULL) ==
                  STARFIX_ERR_INVALID_PARAM,
          "identification input and output capacity");
    check(starfix_identify_stars(3, centroids, 10, 10, 10, &db, 5, matches, scratch, NULL) ==
                  STARFIX_ERR_IDENTIFY_NO_MATCH &&
              starfix_identify_stars(4, centroids, 10, 10, 10, &db, 5, matches, scratch, NULL) ==
                  STARFIX_ERR_IDENTIFY_NO_MATCH,
          "blind identification requires an independent fifth star");
    check(starfix_identify_stars(5, centroids, 10, 10, 10, &db, 5, matches, scratch, NULL) ==
              STARFIX_ERR_NULL_POINTER,
          "blind identification requires hashes");
    centroids[0].u = NAN;
    check(starfix_identify_stars(5, centroids, 10, 10, 10, &db, 5, matches, scratch, NULL) ==
              STARFIX_ERR_DATA_CORRUPTION,
          "nonfinite centroid rejected before hashing");
    centroids[0].u = DBL_MAX;
    check(starfix_identify_stars(5, centroids, 1, 1, 179, &db, 5, matches, scratch, NULL) ==
              STARFIX_ERR_DATA_CORRUPTION,
          "projection overflow rejected before hashing");
    check(starfix_propagate_astrometry(catalog, 0, 2026) == STARFIX_ERR_INVALID_PARAM &&
              starfix_propagate_astrometry(NULL, 5, 2026) == STARFIX_ERR_NULL_POINTER,
          "astrometry pointer and count errors");

    starfix_pose_t poses[2];
    starfix_odometry_t odom = {0};
    double scale = 1, bias = 0;
    check(starfix_solve_graph(1, 1, 1, 0, 0, NULL, 0, NULL, poses, &scale, &bias, scratch) ==
                  STARFIX_ERR_NULL_POINTER &&
              starfix_solve_graph(0, 1, 1, 0, 0, &odom, 0, NULL, poses, &scale, &bias, scratch) ==
                  STARFIX_ERR_INVALID_PARAM,
          "graph pointer and count errors");

    char sentence[256], output[256] = "unchanged";
    check(starfix_format_nmea(34, -118, 0, sentence, sizeof(sentence)) == STARFIX_SUCCESS,
          "NMEA complete output");
    size_t length = strlen(sentence);
    check(starfix_format_nmea(34, -118, 0, output, length) == STARFIX_ERR_BUFFER_FULL &&
              strcmp(output, "unchanged") == 0,
          "NMEA short buffer includes room for the terminator and stays unchanged");
    check(starfix_format_nmea(34, -118, 0, output, length + 1) == STARFIX_SUCCESS &&
              strcmp(output, sentence) == 0,
          "NMEA exact capacity");
    check(starfix_format_nmea(34, -118, 0, output, 0) == STARFIX_ERR_BUFFER_FULL &&
              starfix_format_nmea(34, -118, 0, NULL, 0) == STARFIX_ERR_NULL_POINTER &&
              starfix_format_nmea(NAN, -118, 0, output, sizeof(output)) ==
                  STARFIX_ERR_INVALID_PARAM &&
              starfix_format_nmea(91, 0, 0, output, sizeof(output)) == STARFIX_ERR_INVALID_PARAM,
          "NMEA invalid output and coordinates");
    if (!failures) puts("Public contract checks passed");
    return failures ? 1 : 0;
}
