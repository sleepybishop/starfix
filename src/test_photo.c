#include <math.h>
#include <stdio.h>
#include <string.h>

#include "starfix_photo.h"

int main(void) {
    uint8_t memory[65536];
    starfix_arena_t arena = {memory, memory + sizeof(memory)};
    starfix_photo_options_t options = {
        .width = 400,
        .height = 400,
        .fov_deg = 22.619864948,
        .gha_aries = 120,
        .zenith_camera = {0, 0, 1},
    };
    /* keep inputs ahead of the scratch cursor to exercise nested workspace use */
    starfix_centroid_t* centroids = starfix_alloc_array(&arena, starfix_centroid_t, 6);
    starfix_catalog_star_t* catalog = starfix_alloc_array(&arena, starfix_catalog_star_t, 6);
    starfix_match_t* matches = starfix_alloc_array(&arena, starfix_match_t, 6);
    if (!centroids || !catalog || !matches) return 1;
    uint8_t* scratch = arena.beg;
    for (int i = 0; i < 4; i++) {
        double angle = (double)i * M_PI / 2;
        centroids[i].u = 200 + 50 * cos(angle);
        centroids[i].v = 200 + 50 * sin(angle);
        catalog[i] = starfix_make_catalog_star(i + 1, angle, atan2(1000, 50), 2.0);
        matches[i] = (starfix_match_t){i, i};
    }
    starfix_photo_result_t result;
    for (int i = 0; i < 5; i++) {
        if (starfix_solve_photo(4, centroids, 4, catalog, 4, matches, &options, &result, arena) !=
                STARFIX_SUCCESS ||
            fabs(result.latitude - 90) > 1e-6 || result.rms_pixels > 0.1 || result.inliers != 4 ||
            arena.beg != scratch) {
            fprintf(stderr,
                    "FAILED: repeated photo fits must preserve scratch space and recover pole\n");
            return 1;
        }
    }
    /* verify RANSAC rejects a gross pixel outlier and excludes it from rms_pixels */
    centroids[4] = (starfix_centroid_t){.u = 200.0, .v = 200.0, .flux = 100.0};
    catalog[4] = starfix_make_catalog_star(5, 0.0, M_PI / 2.0, 2.0);
    matches[4] = (starfix_match_t){4, 4};
    centroids[5] = (starfix_centroid_t){
        .u = 200.0 + 50.0 * cos(0.5) + 25.0,
        .v = 200.0 + 50.0 * sin(0.5) - 20.0,
        .flux = 100.0,
    };
    catalog[5] = starfix_make_catalog_star(6, 0.5, atan2(1000, 50), 2.0);
    matches[5] = (starfix_match_t){5, 5};
    if (starfix_solve_photo(6, centroids, 6, catalog, 6, matches, &options, &result, arena) !=
            STARFIX_SUCCESS ||
        result.inliers != 5 || result.rms_pixels > 0.1 || fabs(result.latitude - 90) > 1e-6) {
        fprintf(stderr, "FAILED: RANSAC must reject 25px outlier and exclude it from rms_pixels\n");
        return 1;
    }
    starfix_photo_result_t before = result;
    matches[0].catalog_idx = 4;
    if (starfix_solve_photo(4, centroids, 4, catalog, 4, matches, &options, &result, arena) !=
            STARFIX_ERR_INVALID_PARAM ||
        memcmp(&result, &before, sizeof(result)) != 0 || arena.beg != scratch) {
        fprintf(stderr, "FAILED: invalid match must preserve result and arena\n");
        return 1;
    }
    matches[0].catalog_idx = 0;
    starfix_arena_t small = {scratch, scratch + 1};
    if (starfix_solve_photo(4, centroids, 4, catalog, 4, matches, &options, &result, small) !=
            STARFIX_ERR_MEMORY ||
        memcmp(&result, &before, sizeof(result)) != 0 || small.beg != scratch) {
        fprintf(stderr, "FAILED: insufficient workspace must preserve result and arena\n");
        return 1;
    }
    options.zenith_camera.z = 0;
    if (starfix_solve_photo(4, centroids, 4, catalog, 4, matches, &options, &result, arena) !=
        STARFIX_ERR_INVALID_PARAM) {
        fprintf(stderr, "FAILED: invalid zenith must be rejected\n");
        return 1;
    }
    puts("Photo API workspace reuse and failure handling passed");
    return 0;
}
