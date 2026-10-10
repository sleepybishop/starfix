#include "starfix_photo.h"

#include <math.h>

starfix_status_t starfix_solve_photo(int num_centroids, const starfix_centroid_t* centroids,
                                     uint32_t num_stars, const starfix_catalog_star_t* catalog,
                                     int num_matches, const starfix_match_t* matches,
                                     const starfix_photo_options_t* options,
                                     starfix_photo_result_t* result, starfix_arena_t* arena) {
    if (!centroids || !catalog || !matches || !options || !result || !arena)
        return STARFIX_ERR_NULL_POINTER;
    if (num_centroids < 3 || num_stars < 3 || num_matches < 3 || num_matches > num_centroids ||
        !isfinite(options->width) || options->width <= 0 || !isfinite(options->height) ||
        options->height <= 0 || !isfinite(options->fov_deg) || options->fov_deg <= 0 ||
        options->fov_deg >= 180 || !isfinite(options->gha_aries))
        return STARFIX_ERR_INVALID_PARAM;
    const starfix_vector3_t* zenith = &options->zenith_camera;
    double norm = hypot(hypot(zenith->x, zenith->y), zenith->z);
    if (!isfinite(norm) || fabs(norm - 1.0) > 1e-3) return STARFIX_ERR_INVALID_PARAM;

    /* Copy the cursor so repeated FOV fits can reuse the same workspace. */
    starfix_arena_t scratch = *arena;
    starfix_vector3_t* camera = starfix_alloc_array(&scratch, starfix_vector3_t, num_matches);
    starfix_vector3_t* celestial = starfix_alloc_array(&scratch, starfix_vector3_t, num_matches);
    if (!camera || !celestial) return STARFIX_ERR_MEMORY;
    double focal = options->width / (2.0 * tan(options->fov_deg * M_PI / 360.0));
    for (int i = 0; i < num_matches; i++) {
        int ci = matches[i].centroid_idx, si = matches[i].catalog_idx;
        if (ci < 0 || ci >= num_centroids || si < 0 || (uint32_t)si >= num_stars)
            return STARFIX_ERR_INVALID_PARAM;
        if (!isfinite(centroids[ci].u) || !isfinite(centroids[ci].v))
            return STARFIX_ERR_DATA_CORRUPTION;
        double x = (centroids[ci].u - options->width / 2.0) / focal;
        double y = (centroids[ci].v - options->height / 2.0) / focal;
        norm = hypot(hypot(x, y), 1.0);
        camera[i] = (starfix_vector3_t){x / norm, y / norm, 1.0 / norm};
        double ra = catalog[si].ra, dec = catalog[si].dec;
        if (!isfinite(ra) || !isfinite(dec)) return STARFIX_ERR_DATA_CORRUPTION;
        celestial[i] = (starfix_vector3_t){cos(dec) * cos(ra), cos(dec) * sin(ra), sin(dec)};
    }

    /* A measured zenith permits refraction correction before the single attitude fit. */
    starfix_status_t status = starfix_correct_refraction(num_matches, camera, zenith, 10.0, 1010.0);
    if (status != STARFIX_SUCCESS) return status;
    starfix_photo_result_t fitted = {0};
    starfix_quaternion_t q;
    starfix_telemetry_t telemetry = {0};
    status = starfix_solve_attitude_ransac(num_matches, camera, celestial, 0.25, &q,
                                           fitted.rotation, &scratch, &telemetry);
    if (status != STARFIX_SUCCESS) return status;
    double sum_sq = 0;
    for (int i = 0; i < num_matches; i++) {
        double projected[3];
        for (int row = 0; row < 3; row++) {
            projected[row] = fitted.rotation[row][0] * celestial[i].x +
                             fitted.rotation[row][1] * celestial[i].y +
                             fitted.rotation[row][2] * celestial[i].z;
        }
        if (projected[2] <= 0 || camera[i].z <= 0) return STARFIX_ERR_ATTITUDE_DEGENERATE;
        double du = focal * (projected[0] / projected[2] - camera[i].x / camera[i].z);
        double dv = focal * (projected[1] / projected[2] - camera[i].y / camera[i].z);
        sum_sq += du * du + dv * dv;
    }
    fitted.rms_pixels = sqrt(sum_sq / (double)num_matches);
    fitted.inliers = telemetry.ransac_inliers;
    fitted.matches = num_matches;
    if (!isfinite(fitted.rms_pixels) || fitted.rms_pixels > 3.0 || fitted.inliers < 3 ||
        (double)fitted.inliers / (double)num_matches < 0.6)
        return STARFIX_ERR_RANSAC_INSUFFICIENT_INLIERS;
    status = starfix_solve_position(1, zenith, (const double (*)[3])fitted.rotation,
                                    &options->gha_aries, &fitted.latitude, &fitted.longitude, NULL);
    if (status != STARFIX_SUCCESS) return status;
    *result = fitted;
    return STARFIX_SUCCESS;
}
