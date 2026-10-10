#include "starfix_photo.h"

#include <math.h>

#include "starfix_solver.h"

starfix_status_t starfix_solve_photo(int num_centroids, const starfix_centroid_t* centroids,
                                     uint32_t num_stars, const starfix_catalog_star_t* catalog,
                                     int num_matches, const starfix_match_t* matches,
                                     const starfix_photo_options_t* options,
                                     starfix_photo_result_t* result, starfix_arena_t arena) {
    if (!centroids || !catalog || !matches || !options || !result) return STARFIX_ERR_NULL_POINTER;
    if (num_centroids < 3 || num_stars < 3 || num_matches < 3 || num_matches > num_centroids ||
        !isfinite(options->width) || options->width <= 0 || !isfinite(options->height) ||
        options->height <= 0 || !isfinite(options->fov_deg) || options->fov_deg <= 0 ||
        options->fov_deg >= 180 || !isfinite(options->gha_aries))
        return STARFIX_ERR_INVALID_PARAM;
    const starfix_vector3_t* zenith = &options->zenith_camera;
    double norm = hypot(hypot(zenith->x, zenith->y), zenith->z);
    if (!isfinite(norm) || fabs(norm - 1.0) > 1e-3) return STARFIX_ERR_INVALID_PARAM;

    starfix_vector3_t* camera = starfix_alloc_array(&arena, starfix_vector3_t, num_matches);
    starfix_vector3_t* celestial = starfix_alloc_array(&arena, starfix_vector3_t, num_matches);
    if (!camera || !celestial) return STARFIX_ERR_MEMORY;
    double focal = starfix_focal_length(options->width, options->fov_deg);
    for (int i = 0; i < num_matches; i++) {
        int ci = matches[i].centroid_idx, si = matches[i].catalog_idx;
        if (ci < 0 || ci >= num_centroids || si < 0 || (uint32_t)si >= num_stars)
            return STARFIX_ERR_INVALID_PARAM;
        if (!isfinite(centroids[ci].u) || !isfinite(centroids[ci].v))
            return STARFIX_ERR_DATA_CORRUPTION;
        camera[i] = starfix_camera_ray(centroids[ci].u - options->width / 2.0,
                                       centroids[ci].v - options->height / 2.0, focal);
        double ra = catalog[si].ra, dec = catalog[si].dec;
        if (!isfinite(ra) || !isfinite(dec)) return STARFIX_ERR_DATA_CORRUPTION;
        celestial[i] = starfix_celestial_ray(ra, dec);
    }

    /* A measured zenith permits refraction correction before the single attitude fit. */
    starfix_status_t status = starfix_correct_refraction(num_matches, camera, zenith, 10.0, 1010.0);
    if (status != STARFIX_SUCCESS) return status;
    starfix_photo_result_t fitted = {0};
    int inliers =
        starfix_solve_attitude_ransac(num_matches, camera, celestial, 0.25, fitted.rotation, arena);
    if (inliers < 0) return inliers;
    double sum_sq = 0;
    for (int i = 0; i < num_matches; i++) {
        starfix_vector3_t projected =
            starfix_rotate((const double (*)[3])fitted.rotation, celestial[i]);
        if (projected.z <= 0 || camera[i].z <= 0) return STARFIX_ERR_ATTITUDE_DEGENERATE;
        double du = focal * (projected.x / projected.z - camera[i].x / camera[i].z);
        double dv = focal * (projected.y / projected.z - camera[i].y / camera[i].z);
        sum_sq += du * du + dv * dv;
    }
    fitted.rms_pixels = sqrt(sum_sq / (double)num_matches);
    fitted.inliers = (uint32_t)inliers;
    fitted.matches = num_matches;
    if (!isfinite(fitted.rms_pixels) || fitted.rms_pixels > 3.0 || fitted.inliers < 3 ||
        (double)fitted.inliers / (double)num_matches < 0.6)
        return STARFIX_ERR_RANSAC_INSUFFICIENT_INLIERS;
    status = starfix_solve_position(1, zenith, (const double (*)[3])fitted.rotation,
                                    &options->gha_aries, &fitted.latitude, &fitted.longitude);
    if (status != STARFIX_SUCCESS) return status;
    *result = fitted;
    return STARFIX_SUCCESS;
}
