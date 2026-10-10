#ifndef STARFIX_PHOTO_H
#define STARFIX_PHOTO_H

#include <stdint.h>

#include "starfix_attitude.h"
#include "starfix_identify.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double width, height;
    double fov_deg; /* horizontal field of view */
    double gha_aries;
    starfix_vector3_t zenith_camera; /* measured unit vector: right, down, forward */
} starfix_photo_options_t;

typedef struct {
    double rotation[3][3]; /* celestial -> camera */
    double rms_pixels;
    uint32_t inliers;
    int matches;
    double latitude, longitude;
} starfix_photo_result_t;

/* fit identified matches against a caller-precessed catalog; indices refer to the supplied arrays.
 * uses scratch arena space without advancing its cursor; result is unchanged on error.
 * success requires RMS <= 3 pixels, >= 3 inliers, and inlier fraction >= 0.6 */
starfix_status_t starfix_solve_photo(int num_centroids, const starfix_centroid_t* centroids,
                                     uint32_t num_stars, const starfix_catalog_star_t* catalog,
                                     int num_matches, const starfix_match_t* matches,
                                     const starfix_photo_options_t* options,
                                     starfix_photo_result_t* result, starfix_arena_t arena);

#ifdef __cplusplus
}
#endif
#endif
