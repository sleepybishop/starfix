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
    starfix_centroid_t centroids[4] = {{0}};
    starfix_catalog_star_t catalog[4] = {{0}};
    starfix_match_t matches[4];
    for (int i = 0; i < 4; i++) {
        double angle = (double)i * M_PI / 2;
        centroids[i].u = 200 + 50 * cos(angle);
        centroids[i].v = 200 + 50 * sin(angle);
        catalog[i].ra = angle;
        catalog[i].dec = atan2(1000, 50);
        matches[i] = (starfix_match_t){i, i};
    }
    starfix_photo_result_t result;
    for (int i = 0; i < 5; i++) {
        if (starfix_solve_photo(4, centroids, 4, catalog, 4, matches, &options, &result, &arena) !=
                STARFIX_SUCCESS ||
            fabs(result.latitude - 90) > 1e-6 || result.rms_pixels > 0.1 || result.inliers != 4 ||
            arena.beg != memory) {
            fprintf(stderr,
                    "FAILED: repeated photo fits must preserve scratch space and recover pole\n");
            return 1;
        }
    }
    starfix_photo_result_t before = result;
    matches[0].catalog_idx = 4;
    if (starfix_solve_photo(4, centroids, 4, catalog, 4, matches, &options, &result, &arena) !=
            STARFIX_ERR_INVALID_PARAM ||
        memcmp(&result, &before, sizeof(result)) != 0 || arena.beg != memory) {
        fprintf(stderr, "FAILED: invalid match must preserve result and arena\n");
        return 1;
    }
    matches[0].catalog_idx = 0;
    starfix_arena_t small = {memory, memory + 1};
    if (starfix_solve_photo(4, centroids, 4, catalog, 4, matches, &options, &result, &small) !=
            STARFIX_ERR_MEMORY ||
        memcmp(&result, &before, sizeof(result)) != 0 || small.beg != memory) {
        fprintf(stderr, "FAILED: insufficient workspace must preserve result and arena\n");
        return 1;
    }
    options.zenith_camera.z = 0;
    if (starfix_solve_photo(4, centroids, 4, catalog, 4, matches, &options, &result, &arena) !=
        STARFIX_ERR_INVALID_PARAM) {
        fprintf(stderr, "FAILED: invalid zenith must be rejected\n");
        return 1;
    }
    puts("Photo API workspace reuse and failure handling passed");
    return 0;
}
