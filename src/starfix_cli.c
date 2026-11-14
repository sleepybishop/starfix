#include <limits.h>
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
#include "starfix_photo.h"

static uint8_t mempool[65536];

static int numeric_argument(const char* text, const char* end) {
    if (end && !*end) return 1;
    fprintf(stderr, "Invalid numeric argument: %s\n", text);
    return 0;
}

static int load_catalog(starfix_db_t* db, double year) {
    if (!isfinite(year)) return STARFIX_ERR_INVALID_PARAM;
    int status = starfix_load_db("data/starfix_db.bin", db);
    if (status < 0) status = starfix_load_db("../data/starfix_db.bin", db);
    if (status < 0) {
        fprintf(stderr, "Could not load data/starfix_db.bin or ../data/starfix_db.bin.\n");
        return status;
    }
    return starfix_propagate_astrometry(db->stars, db->num_stars, year);
}

static int do_centroid(int threshold) {
    printf("Loading image from data/mock_image.pgm...\n");
    int w = 0, h = 0;
    unsigned char* img = read_pgm("data/mock_image.pgm", &w, &h);
    if (img == NULL) {
        printf("Error loading PGM image.\n");
        return -1;
    }

    starfix_centroid_t centroids[1000];
    int count =
        starfix_find_centroids(img, w, h, (unsigned char)threshold, 1000, centroids, NULL, 1);
    free(img);

    if (count < 0) {
        printf("Centroid extraction failed with code %d\n", count);
        return -1;
    }

    printf("Detected %d star centroids.\n", count);
    if (write_centroids_json("data/detected_centroids.json", centroids, count) < 0) {
        fprintf(stderr, "Could not write data/detected_centroids.json.\n");
        return 1;
    }
    return 0;
}

static int do_identify(int width, int height, double fov_deg, double observation_year) {
    starfix_db_t db = {0};
    int result = -1;
    if (load_catalog(&db, observation_year) != 0) goto done;

    starfix_centroid_t centroids[1000];
    int cent_count = read_centroids_json("data/detected_centroids.json", centroids, 1000);
    if (cent_count <= 0) {
        printf("Error loading detected centroids.\n");
        goto done;
    }

    printf("Using top %d brightest centroids for pattern matching...\n",
           (cent_count < 12) ? cent_count : 12);
    starfix_match_t matches[20];
    int match_count = starfix_identify_stars(cent_count, centroids, width, height, fov_deg, &db, 20,
                                             matches, NULL);

    if (match_count < 0) {
        printf("Identification failed with code %d\n", match_count);
        goto done;
    }

    printf("Identified Stars:\n");
    printf("Centroid Index  | HIP ID     | RA (deg)   | Dec (deg)\n");
    printf("-----------------------------------------------------\n");
    for (int i = 0; i < match_count; i++) {
        int c_idx = matches[i].centroid_idx;
        int cat_idx = matches[i].catalog_idx;
        printf("%-15d | %-10d | %10.4f | %10.4f\n", c_idx, db.stars[cat_idx].hip,
               db.stars[cat_idx].ra * 180.0 / M_PI, db.stars[cat_idx].dec * 180.0 / M_PI);
    }

    if (write_identified_json("data/identified_stars.json", matches, match_count, db.stars) < 0) {
        fprintf(stderr, "Could not write data/identified_stars.json.\n");
        goto done;
    }
    printf("Saved %d identified stars to data/identified_stars.json\n", match_count);

    result = 0;
done:
    starfix_free_db(&db);
    return result;
}

static void print_solution(const starfix_photo_result_t* fitted, int attitude_only) {
    const double (*R)[3] = (const double (*)[3])fitted->rotation;
    printf("RMS_ERROR: %.9f pixels\n", fitted->rms_pixels);
    printf("INLIERS: %u/%d\n", fitted->inliers, fitted->matches);

    double ra_deg, dec_deg, roll_deg;
    starfix_rotation_to_euler(R, &ra_deg, &dec_deg, &roll_deg);

    printf("\nEstimated Attitude:\n");
    printf("  Right Ascension : %.4f degrees\n", ra_deg);
    printf("  Declination     : %.4f degrees\n", dec_deg);
    printf("  Roll            : %.4f degrees\n", roll_deg);
    printf("ROTATION_MATRIX: %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f\n", R[0][0],
           R[0][1], R[0][2], R[1][0], R[1][1], R[1][2], R[2][0], R[2][1], R[2][2]);

    if (!attitude_only) {
        printf("\n>>> CELESTIAL POSITION SOLVED! <<<\n");
        printf("  Solved Latitude  : %f degrees North\n", fitted->latitude);
        printf("  Solved Longitude : %f degrees East\n", fitted->longitude);
    }
}

static int do_solve_photo(const starfix_photo_options_t* options, double year, int attitude_only) {
    starfix_db_t db = {0};
    int result = 1;
    if (load_catalog(&db, year) != 0) goto done;

    starfix_centroid_t centroids[1000];
    int cent_count = read_centroids_json("data/detected_centroids.json", centroids, 1000);
    if (cent_count == -1) {
        cent_count = read_centroids_json("../data/detected_centroids.json", centroids, 1000);
    }
    if (cent_count <= 0) {
        printf("Error: Could not load centroids from data/detected_centroids.json\n");
        goto done;
    }

    starfix_match_t matches[1000];
    int match_count = read_identified_json("data/identified_stars.json", matches, 1000);
    if (match_count == -1) {
        match_count = read_identified_json("../data/identified_stars.json", matches, 1000);
    }
    if (match_count < 3) {
        printf(
            "Error: Need at least 3 identified stars in data/identified_stars.json (found: %d)\n",
            match_count);
        goto done;
    }

    starfix_photo_result_t fitted;
    starfix_arena_t arena = {mempool, mempool + sizeof(mempool)};
    starfix_status_t status = starfix_solve_photo(cent_count, centroids, db.num_stars, db.stars,
                                                  match_count, matches, options, &fitted, arena);
    if (status != STARFIX_SUCCESS) {
        fprintf(stderr, "Photo solution failed quality or input checks (status %d)\n", status);
        goto done;
    }
    print_solution(&fitted, attitude_only);
    result = 0;
done:
    starfix_free_db(&db);
    return result;
}

static double eval_sweep_fov(int cent_count, const starfix_centroid_t* centroids,
                             const starfix_db_t* db, int match_count,
                             const starfix_match_t* matches, const starfix_photo_options_t* base,
                             double fov_deg, uint32_t min_inliers, starfix_arena_t arena,
                             starfix_photo_result_t* out_fit) {
    starfix_photo_options_t opt = *base;
    opt.fov_deg = fov_deg;
    if (starfix_solve_photo(cent_count, centroids, db->num_stars, db->stars, match_count, matches,
                            &opt, out_fit, arena) != STARFIX_SUCCESS ||
        out_fit->inliers < min_inliers) {
        return HUGE_VAL;
    }
    return out_fit->rms_pixels;
}

static int do_sweep_photo(const starfix_photo_options_t* base, double year) {
    starfix_db_t db = {0};
    int result = 1;
    if (load_catalog(&db, year) != 0) goto done;

    starfix_centroid_t centroids[1000];
    int cent_count = read_centroids_json("data/detected_centroids.json", centroids, 1000);
    if (cent_count == -1)
        cent_count = read_centroids_json("../data/detected_centroids.json", centroids, 1000);
    if (cent_count < 5) goto done;

    double low_search = fmax(2.0, base->fov_deg * 0.75);
    double high_search = fmin(90.0, base->fov_deg * 1.25);
    double step = (high_search - low_search) / 40.0;

    starfix_match_t best_matches[20];
    int best_match_count = 0;
    starfix_photo_result_t best_fit = {0};
    double best_fov = base->fov_deg;
    starfix_arena_t arena = {mempool, mempool + sizeof(mempool)};

    for (int i = 0; i <= 40; i++) {
        starfix_photo_options_t opt = *base;
        opt.fov_deg = low_search + i * step;
        starfix_match_t matches[20];
        int match_count = starfix_identify_stars(cent_count, centroids, opt.width, opt.height,
                                                 opt.fov_deg, &db, 20, matches, NULL);
        if (match_count < 3) continue;
        starfix_photo_result_t fit;
        if (starfix_solve_photo(cent_count, centroids, db.num_stars, db.stars, match_count, matches,
                                &opt, &fit, arena) != STARFIX_SUCCESS)
            continue;
        if (best_match_count == 0 || fit.inliers > best_fit.inliers ||
            (fit.inliers == best_fit.inliers && fit.rms_pixels < best_fit.rms_pixels)) {
            best_fit = fit;
            best_fov = opt.fov_deg;
            best_match_count = match_count;
            memcpy(best_matches, matches, (size_t)match_count * sizeof(*matches));
        }
    }
    if (best_match_count == 0) goto done;

    uint32_t target_inliers = best_fit.inliers;
    double half_width = fmax(step, best_fov * 0.025);
    double low = fmax(1.0, best_fov - half_width);
    double high = fmin(100.0, best_fov + half_width);
    double ratio = (sqrt(5.0) - 1.0) / 2.0;
    double x1 = high - ratio * (high - low);
    double x2 = low + ratio * (high - low);
    starfix_photo_result_t fit1 = {0}, fit2 = {0};
    double rms1 = eval_sweep_fov(cent_count, centroids, &db, best_match_count, best_matches, base,
                                 x1, target_inliers, arena, &fit1);
    double rms2 = eval_sweep_fov(cent_count, centroids, &db, best_match_count, best_matches, base,
                                 x2, target_inliers, arena, &fit2);
    for (int iter = 0; iter < 18; iter++) {
        if (rms1 <= rms2) {
            high = x2;
            x2 = x1;
            fit2 = fit1;
            rms2 = rms1;
            x1 = high - ratio * (high - low);
            rms1 = eval_sweep_fov(cent_count, centroids, &db, best_match_count, best_matches, base,
                                  x1, target_inliers, arena, &fit1);
        } else {
            low = x1;
            x1 = x2;
            fit1 = fit2;
            rms1 = rms2;
            x2 = low + ratio * (high - low);
            rms2 = eval_sweep_fov(cent_count, centroids, &db, best_match_count, best_matches, base,
                                  x2, target_inliers, arena, &fit2);
        }
    }
    if (rms1 < best_fit.rms_pixels) {
        best_fit = fit1;
        best_fov = x1;
    }
    if (rms2 < best_fit.rms_pixels) {
        best_fit = fit2;
        best_fov = x2;
    }
    if (write_identified_json("data/identified_stars.json", best_matches, best_match_count,
                              db.stars) < 0)
        goto done;
    printf("SOLVED_FOV: %.9f\n", best_fov);
    print_solution(&best_fit, 0);
    result = 0;
done:
    starfix_free_db(&db);
    return result;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        printf("StarFix Native C Command-Line Interface\n");
        printf("Usage: %s [option]\n", argv[0]);
        printf("Options:\n");
        printf("  --centroid [threshold]  Extract star centroids; threshold is 0 to 255\n");
        printf("  --identify [width height fov observation_year]\n");
        printf("                   Match detected centroids to database catalog quads\n");
        printf("  --fit-photo <g_x> <g_y> <g_z> [fov_deg] [observation_year] [crop_size]\n");
        printf("                   Fit photo attitude using a measured camera zenith\n");
        printf(
            "  --solve-photo    Solve terrestrial position from matched centroids (uses: "
            "<gha_aries> <g_x> <g_y> <g_z> [fov_deg] [observation_year] [crop_size])\n");
        printf(
            "  --sweep-photo    Sweep FOV, identify stars, and solve terrestrial position "
            "(uses: <gha_aries> <g_x> <g_y> <g_z> [fov_deg] [observation_year] [crop_size])\n");
        return 1;
    }

    if (strcmp(argv[1], "--centroid") == 0) {
        if (argc > 3) {
            fprintf(stderr, "Usage: %s --centroid [threshold]\n", argv[0]);
            return 1;
        }
        int threshold = 35;
        if (argc == 3 && !numeric_argument(argv[2], parse_int(argv[2], 0, 255, &threshold)))
            return 1;
        return do_centroid(threshold);
    } else if (strcmp(argv[1], "--identify") == 0) {
        if (argc > 6) {
            fprintf(stderr, "Usage: %s --identify [width height fov observation_year]\n", argv[0]);
            return 1;
        }
        int w = 1024, h = 1024;
        double fov = 12.0, observation_year = current_year_utc();
        if ((argc >= 3 && !numeric_argument(argv[2], parse_int(argv[2], 1, INT_MAX, &w))) ||
            (argc >= 4 && !numeric_argument(argv[3], parse_int(argv[3], 1, INT_MAX, &h))) ||
            (argc >= 5 && !numeric_argument(argv[4], parse_double(argv[4], &fov))) ||
            (argc >= 6 && !numeric_argument(argv[5], parse_double(argv[5], &observation_year))))
            return 1;
        return do_identify(w, h, fov, observation_year);
    } else if (strcmp(argv[1], "--solve-photo") == 0 || strcmp(argv[1], "--fit-photo") == 0 ||
               strcmp(argv[1], "--sweep-photo") == 0) {
        int sweep = strcmp(argv[1], "--sweep-photo") == 0;
        int attitude_only = strcmp(argv[1], "--fit-photo") == 0;
        int first = attitude_only ? 2 : 3;
        if (argc < first + 3 || argc > first + 6) {
            printf("Usage: %s %s %s<g_x> <g_y> <g_z> [fov_deg] [observation_year] [crop_size]\n",
                   argv[0], argv[1], attitude_only ? "" : "<gha_aries> ");
            return 1;
        }
        double year = current_year_utc();
        starfix_photo_options_t options = {.width = 560.0, .fov_deg = 11.28};
        double* fields[] = {&options.gha_aries,
                            &options.zenith_camera.x,
                            &options.zenith_camera.y,
                            &options.zenith_camera.z,
                            &options.fov_deg,
                            &year,
                            &options.width};
        int offset = attitude_only ? 1 : 0;
        for (int i = 2; i < argc; i++)
            if (!numeric_argument(argv[i], parse_double(argv[i], fields[i - 2 + offset]))) return 1;
        options.height = options.width;
        return sweep ? do_sweep_photo(&options, year)
                     : do_solve_photo(&options, year, attitude_only);
    } else {
        printf("Invalid option: %s\n", argv[1]);
        return 1;
    }
}
