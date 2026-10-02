#ifndef STARFIX_PHOTO_H
#define STARFIX_PHOTO_H

#include "starfix_attitude.h"
#include "starfix_identify.h"

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

/* Fit already-identified observations and compute a direct position fix.
 * The caller owns/precesses the catalog and can reuse it across FOV trials.
 * Refraction, RANSAC, projection checks, and quality limits are shared by all callers.
 * No file I/O, logging, or heap allocations. The supplied arena is scratch space;
 * its cursor and result are unchanged on failure, and its cursor is also unchanged
 * on success. A successful result includes RMS <= 3 pixels, >= 3 inliers, and an
 * inlier fraction >= 0.6. Match indices refer to the supplied arrays. */
starfix_status_t starfix_solve_photo(int num_centroids, const starfix_centroid_t* centroids,
                                     uint32_t num_stars, const starfix_catalog_star_t* catalog,
                                     int num_matches, const starfix_match_t* matches,
                                     const starfix_photo_options_t* options,
                                     starfix_photo_result_t* result, starfix_arena_t* arena);

#endif
