#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "starfix_arena.h"
#include "starfix_attitude.h"
#include "starfix_centroid.h"
#include "starfix_cli_io.h"
#include "starfix_identify.h"
#include "starfix_solver.h"

/* 10 MB statically allocated memory pool for the entire pipeline */
static uint8_t mempool[10 * 1024 * 1024];
static starfix_arena_t arena;

#define RESET_ARENA()                          \
    do {                                       \
        arena.beg = mempool;                   \
        arena.end = mempool + sizeof(mempool); \
    } while (0)

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* actions */
static int do_centroid(int threshold) {
    printf("Loading image from data/mock_image.pgm...\n");
    int w = 0, h = 0;
    unsigned char* img = read_pgm("data/mock_image.pgm", &w, &h);
    if (img == NULL) {
        printf("Error loading PGM image.\n");
        return -1;
    }

    starfix_centroid_t centroids[1000];
    starfix_telemetry_t telem = {0};
    starfix_status_t count_status = starfix_find_centroids(img, w, h, (unsigned char)threshold,
                                                           1000, centroids, NULL, 1, &telem);
    int count = (int)telem.num_stars_detected;
    free(img);

    if (count_status != STARFIX_SUCCESS) {
        printf("Centroid extraction failed with code %d\n", (int)count_status);
        return -1;
    }

    printf("Detected %d star centroids.\n", count);
    write_centroids_json("data/detected_centroids.json", centroids, count);
    return 0;
}

static int do_identify(int width, int height, double fov_deg, double observation_year) {
    starfix_catalog_star_t* catalog = NULL;
    starfix_hash_entry_t* hash_table = NULL;
    uint32_t num_stars = 0, num_entries = 0, bin_factor = 0;

    printf("Loading database from data/starfix_db.bin...\n");
    if (starfix_load_db("data/starfix_db.bin", &catalog, &num_stars, &hash_table, &num_entries,
                        &bin_factor) != 0) {
        if (starfix_load_db("../data/starfix_db.bin", &catalog, &num_stars, &hash_table,
                            &num_entries, &bin_factor) != 0) {
            printf("Error loading binary database.\n");
            return -1;
        }
    }

    printf("Propagating catalog precession to year %.2f...\n", observation_year);
    starfix_propagate_astrometry(catalog, num_stars, observation_year, NULL);

    starfix_centroid_t centroids[1000];
    int cent_count = read_centroids_json("data/detected_centroids.json", centroids, 1000);
    if (cent_count <= 0) {
        printf("Error loading detected centroids.\n");
        free(catalog);
        free(hash_table);
        return -1;
    }

    printf("Using top %d brightest centroids for pattern matching...\n",
           (cent_count < 8) ? cent_count : 8);
    starfix_match_t matches[20];
    starfix_telemetry_t telem = {0};
    RESET_ARENA();
    starfix_status_t match_status = starfix_identify_stars(
        cent_count, centroids, width, height, fov_deg, num_stars, catalog, num_entries, hash_table,
        bin_factor, 20, matches, &arena, &telem, NULL);
    int match_count = (int)telem.identify_matches;

    if (match_status != STARFIX_SUCCESS) {
        printf("Identification failed with code %d\n", (int)match_status);
        free(catalog);
        free(hash_table);
        return -1;
    }

    printf("Identified Stars:\n");
    printf("Centroid Index  | HIP ID     | RA (deg)   | Dec (deg)\n");
    printf("-----------------------------------------------------\n");
    int i;
    for (i = 0; i < match_count; i++) {
        int c_idx = matches[i].centroid_idx;
        int cat_idx = matches[i].catalog_idx;
        printf("%-15d | %-10d | %10.4f | %10.4f\n", c_idx, catalog[cat_idx].hip,
               starfix_catalog_ra(&catalog[cat_idx]) * 180.0 / M_PI,
               starfix_catalog_dec(&catalog[cat_idx]) * 180.0 / M_PI);
    }

    write_identified_json("data/identified_stars.json", matches, match_count, catalog);
    printf("Saved %d identified stars to data/identified_stars.json\n", match_count);

    free(catalog);
    free(hash_table);
    return 0;
}

static int do_solve_photo(double gha_aries, double g_x, double g_y, double g_z, double fov_deg,
                          double observation_year, int attitude_only, double crop_size) {
    if (!isfinite(gha_aries) || !isfinite(fov_deg) || !isfinite(observation_year) ||
        fov_deg <= 0.0 || fov_deg >= 180.0 || !isfinite(crop_size) || crop_size <= 0.0) {
        printf("Invalid photo attitude parameters.\n");
        return 1;
    }
    {
        double gravity_norm = sqrt(g_x * g_x + g_y * g_y + g_z * g_z);
        if (!isfinite(gravity_norm) || fabs(gravity_norm - 1.0) > 1e-3) {
            printf("Camera zenith vector must be normalized.\n");
            return 1;
        }
    }

    starfix_catalog_star_t* catalog = NULL;
    starfix_hash_entry_t* hash_table = NULL;
    uint32_t num_stars = 0, num_entries = 0, bin_factor = 0;

    if (starfix_load_db("data/starfix_db.bin", &catalog, &num_stars, &hash_table, &num_entries,
                        &bin_factor) != 0) {
        if (starfix_load_db("../data/starfix_db.bin", &catalog, &num_stars, &hash_table,
                            &num_entries, &bin_factor) != 0) {
            printf("Failed to load db (tried data/starfix_db.bin and ../data/starfix_db.bin).\n");
            return 1;
        }
    }
    if (starfix_propagate_astrometry(catalog, num_stars, observation_year, NULL) != 0) {
        printf("Failed to propagate catalog astrometry.\n");
        free(catalog);
        free(hash_table);
        return 1;
    }

    /* load detected centroids */
    starfix_centroid_t centroids[1000];
    int cent_count = read_centroids_json("data/detected_centroids.json", centroids, 1000);
    if (cent_count <= 0) {
        cent_count = read_centroids_json("../data/detected_centroids.json", centroids, 1000);
    }
    if (cent_count <= 0) {
        printf("Error: Could not load centroids from data/detected_centroids.json\n");
        free(catalog);
        free(hash_table);
        return 1;
    }

    /* load identified stars */
    starfix_match_t matches[1000];
    int match_count = read_identified_json("data/identified_stars.json", matches, 1000);
    if (match_count <= 0) {
        match_count = read_identified_json("../data/identified_stars.json", matches, 1000);
    }
    if (match_count < 3) {
        printf(
            "Error: Need at least 3 identified stars in data/identified_stars.json (found: %d)\n",
            match_count);
        free(catalog);
        free(hash_table);
        return 1;
    }

    double fov_rad = fov_deg * M_PI / 180.0;
    double f_len = crop_size / (2.0 * tan(fov_rad / 2.0));

    starfix_vector3_t* cam_vecs =
        (starfix_vector3_t*)malloc((size_t)match_count * sizeof(starfix_vector3_t));
    starfix_vector3_t* raw_cam_vecs =
        (starfix_vector3_t*)malloc((size_t)match_count * sizeof(starfix_vector3_t));
    starfix_vector3_t* cel_vecs =
        (starfix_vector3_t*)malloc((size_t)match_count * sizeof(starfix_vector3_t));
    if (cam_vecs == NULL || raw_cam_vecs == NULL || cel_vecs == NULL) {
        free(catalog);
        free(hash_table);
        free(cam_vecs);
        free(raw_cam_vecs);
        free(cel_vecs);
        return 1;
    }

    int i;
    for (i = 0; i < match_count; i++) {
        int c_idx = matches[i].centroid_idx;
        int cat_idx = matches[i].catalog_idx;
        if (c_idx < 0 || c_idx >= cent_count || cat_idx < 0 || (uint32_t)cat_idx >= num_stars) {
            free(catalog);
            free(hash_table);
            free(cam_vecs);
            free(raw_cam_vecs);
            free(cel_vecs);
            return 1;
        }

        /* convert centroids to camera frame unit vectors */
        double cx = (centroids[c_idx].u - crop_size / 2.0) / f_len;
        double cy = (centroids[c_idx].v - crop_size / 2.0) / f_len;
        double cz = 1.0;
        double norm = sqrt(cx * cx + cy * cy + cz * cz);
        cam_vecs[i].x = cx / norm;
        cam_vecs[i].y = cy / norm;
        cam_vecs[i].z = cz / norm;
        raw_cam_vecs[i] = cam_vecs[i];

        /* convert catalog coordinates to celestial frame unit vectors */
        double dec = starfix_catalog_dec(&catalog[cat_idx]);
        double ra = starfix_catalog_ra(&catalog[cat_idx]);
        cel_vecs[i].x = cos(dec) * cos(ra);
        cel_vecs[i].y = cos(dec) * sin(ra);
        cel_vecs[i].z = sin(dec);
    }

    /* pass 1: estimate a rough camera attitude matrix R */
    double R[3][3];
    starfix_quaternion_t q_est;
    starfix_telemetry_t attitude_telem = {0};
    RESET_ARENA();
    starfix_status_t status = starfix_solve_attitude_ransac(match_count, cam_vecs, cel_vecs, 0.25,
                                                            &q_est, R, &arena, &attitude_telem);

    if (status == STARFIX_SUCCESS) {
        printf("Applying atmospheric refraction correction (measured zenith)...\n");
        {
            starfix_vector3_t z_cam = {g_x, g_y, g_z};
            for (i = 0; i < match_count; i++) cam_vecs[i] = raw_cam_vecs[i];
            if (starfix_correct_refraction(match_count, cam_vecs, &z_cam, 10.0, 1010.0) !=
                STARFIX_SUCCESS) {
                status = STARFIX_ERR_INVALID_PARAM;
            }
            attitude_telem = (starfix_telemetry_t){0};
            RESET_ARENA();
            if (status == STARFIX_SUCCESS)
                status = starfix_solve_attitude_ransac(match_count, cam_vecs, cel_vecs, 0.25,
                                                       &q_est, R, &arena, &attitude_telem);
        }
    }

    /* calculate RMS against the corrected observations. every matched star must project in front
     * of the camera; otherwise a skipped projection could make an invalid solution look perfect */
    double sum_sq_err = 0.0;
    double rms_err = INFINITY;
    if (status == STARFIX_SUCCESS) {
        for (i = 0; i < match_count; i++) {
            double vx = R[0][0] * cel_vecs[i].x + R[0][1] * cel_vecs[i].y + R[0][2] * cel_vecs[i].z;
            double vy = R[1][0] * cel_vecs[i].x + R[1][1] * cel_vecs[i].y + R[1][2] * cel_vecs[i].z;
            double vz = R[2][0] * cel_vecs[i].x + R[2][1] * cel_vecs[i].y + R[2][2] * cel_vecs[i].z;
            if (vz <= 0.0 || cam_vecs[i].z <= 0.0) {
                status = STARFIX_ERR_INVALID_PARAM;
                break;
            }
            double u_proj = crop_size / 2.0 + (vx / vz) * f_len;
            double v_proj = crop_size / 2.0 + (vy / vz) * f_len;
            double u_observed = crop_size / 2.0 + (cam_vecs[i].x / cam_vecs[i].z) * f_len;
            double v_observed = crop_size / 2.0 + (cam_vecs[i].y / cam_vecs[i].z) * f_len;
            double du = u_proj - u_observed;
            double dv = v_proj - v_observed;
            sum_sq_err += du * du + dv * dv;
        }
        if (status == STARFIX_SUCCESS) {
            rms_err = sqrt(sum_sq_err / (double)match_count);
            printf("RMS_ERROR: %.9f pixels\n", rms_err);
            printf("INLIERS: %u/%d\n", attitude_telem.ransac_inliers, match_count);
        }
    }

    free(cam_vecs);
    free(raw_cam_vecs);
    free(cel_vecs);

    if (status != STARFIX_SUCCESS || !isfinite(rms_err) || rms_err > 3.0 ||
        attitude_telem.ransac_inliers < 3 ||
        (double)attitude_telem.ransac_inliers / (double)match_count < 0.6) {
        printf("Attitude failed quality checks (status %d, RMS %.6f, inliers %u/%d)\n", status,
               rms_err, attitude_telem.ransac_inliers, match_count);
        free(catalog);
        free(hash_table);
        return 1;
    }

    /* Convert R to Euler Angles */
    double R_flat[9] = {R[0][0], R[0][1], R[0][2], R[1][0], R[1][1],
                        R[1][2], R[2][0], R[2][1], R[2][2]};
    double ra_deg, dec_deg, roll_deg;
    starfix_rotation_to_euler(R_flat, &ra_deg, &dec_deg, &roll_deg);

    printf("\nEstimated Attitude:\n");
    printf("  Right Ascension : %.4f degrees\n", ra_deg);
    printf("  Declination     : %.4f degrees\n", dec_deg);
    printf("  Roll            : %.4f degrees\n", roll_deg);
    printf("ROTATION_MATRIX: %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f\n", R[0][0],
           R[0][1], R[0][2], R[1][0], R[1][1], R[1][2], R[2][0], R[2][1], R[2][2]);

    if (attitude_only) {
        free(catalog);
        free(hash_table);
        return 0;
    }

    /* Solve absolute position using command line parameters */
    double solved_lat = 0.0;
    double solved_lon = 0.0;

    starfix_vector3_t zenith_meas[1];
    zenith_meas[0].x = g_x;
    zenith_meas[0].y = g_y;
    zenith_meas[0].z = g_z;

    RESET_ARENA();
    starfix_telemetry_t telem = {0};
    starfix_status_t solve_status = starfix_solve_position(
        1, zenith_meas, (const double (*)[3])R, &gha_aries, &solved_lat, &solved_lon, &telem);

    if (solve_status == STARFIX_SUCCESS) {
        printf("\n>>> CELESTIAL POSITION SOLVED! <<<\n");
        printf("  Solved Latitude  : %f degrees North\n", solved_lat);
        printf("  Solved Longitude : %f degrees East\n", solved_lon);

    } else {
        printf("Position fix solver failed with status %d\n", solve_status);
    }

    free(catalog);
    free(hash_table);
    return (solve_status == STARFIX_SUCCESS) ? 0 : 1;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        printf("StarFix Native C Command-Line Interface\n");
        printf("Usage: %s [option]\n", argv[0]);
        printf("Options:\n");
        printf("  --centroid       Extract sub-pixel star centroids from raw PGM image\n");
        printf("  --identify [width height fov observation_year]\n");
        printf("                   Match detected centroids to database catalog quads\n");
        printf("  --fit-photo <g_x> <g_y> <g_z> [fov_deg] [observation_year] [crop_size]\n");
        printf("                   Fit photo attitude using a measured camera zenith\n");
        printf(
            "  --solve-photo    Solve terrestrial position from a real photo (uses: <gha_aries> "
            "<g_x> <g_y> <g_z> [fov_deg] [observation_year] [crop_size])\n");
        return 1;
    }

    if (strcmp(argv[1], "--centroid") == 0) {
        int threshold = (argc >= 3) ? atoi(argv[2]) : 35;
        return do_centroid(threshold);
    } else if (strcmp(argv[1], "--identify") == 0) {
        int w = (argc >= 3) ? atoi(argv[2]) : 1024;
        int h = (argc >= 4) ? atoi(argv[3]) : 1024;
        double fov = (argc >= 5) ? atof(argv[4]) : 12.0;
        double observation_year;
        if (argc >= 6) {
            observation_year = atof(argv[5]);
        } else {
            time_t now = time(NULL);
            struct tm* now_utc = gmtime(&now);
            if (now_utc == NULL) return 1;
            observation_year = 1900.0 + now_utc->tm_year + now_utc->tm_yday / 365.25;
        }
        return do_identify(w, h, fov, observation_year);
    } else if (strcmp(argv[1], "--solve-photo") == 0) {
        if (argc < 6 || argc > 9) {
            printf(
                "Usage: %s --solve-photo <gha_aries> <g_x> <g_y> <g_z> "
                "[fov_deg] [observation_year] [crop_size]\n",
                argv[0]);
            return 1;
        }
        double gha = atof(argv[2]);
        double gx = atof(argv[3]);
        double gy = atof(argv[4]);
        double gz = atof(argv[5]);
        double fov_deg = (argc >= 7) ? atof(argv[6]) : 11.28;
        double observation_year;
        if (argc >= 8) {
            observation_year = atof(argv[7]);
        } else {
            time_t now = time(NULL);
            struct tm* now_utc = gmtime(&now);
            if (now_utc == NULL) return 1;
            observation_year = 1900.0 + now_utc->tm_year + now_utc->tm_yday / 365.25;
        }
        double crop_size = (argc >= 9) ? atof(argv[8]) : 560.0;
        return do_solve_photo(gha, gx, gy, gz, fov_deg, observation_year, 0, crop_size);
    } else if (strcmp(argv[1], "--fit-photo") == 0) {
        if (argc < 5) {
            printf(
                "Usage: %s --fit-photo <g_x> <g_y> <g_z> [fov_deg] [observation_year] "
                "[crop_size]\n",
                argv[0]);
            return 1;
        }
        double gx = atof(argv[2]), gy = atof(argv[3]), gz = atof(argv[4]);
        double fov_deg = (argc >= 6) ? atof(argv[5]) : 11.28;
        double observation_year;
        if (argc >= 7) {
            observation_year = atof(argv[6]);
        } else {
            time_t now = time(NULL);
            struct tm* now_utc = gmtime(&now);
            if (now_utc == NULL) return 1;
            observation_year = 1900.0 + now_utc->tm_year + now_utc->tm_yday / 365.25;
        }
        double crop_size = (argc >= 8) ? atof(argv[7]) : 560.0;
        return do_solve_photo(0.0, gx, gy, gz, fov_deg, observation_year, 1, crop_size);
    } else {
        printf("Invalid option: %s\n", argv[1]);
        return 1;
    }
}
