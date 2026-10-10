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

static uint8_t mempool[10 * 1024 * 1024];

static int numeric_argument(const char* text, const char* end) {
    if (end && !*end) return 1;
    fprintf(stderr, "Invalid numeric argument: %s\n", text);
    return 0;
}

static double current_year(void) {
    time_t now = time(NULL);
    struct tm* utc = gmtime(&now);
    return utc ? 1900.0 + utc->tm_year + utc->tm_yday / 365.25 : NAN;
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
    starfix_arena_t arena = {mempool, mempool + sizeof(mempool)};
    int match_count = starfix_identify_stars(cent_count, centroids, width, height, fov_deg, &db, 20,
                                             matches, arena, NULL);

    if (match_count < 0) {
        printf("Identification failed with code %d\n", match_count);
        goto done;
    }

    printf("Identified Stars:\n");
    printf("Centroid Index  | HIP ID     | RA (deg)   | Dec (deg)\n");
    printf("-----------------------------------------------------\n");
    int i;
    for (i = 0; i < match_count; i++) {
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
    const double (*R)[3] = (const double (*)[3])fitted.rotation;
    printf("RMS_ERROR: %.9f pixels\n", fitted.rms_pixels);
    printf("INLIERS: %u/%d\n", fitted.inliers, fitted.matches);

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
        printf("  Solved Latitude  : %f degrees North\n", fitted.latitude);
        printf("  Solved Longitude : %f degrees East\n", fitted.longitude);
    }
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
            "  --solve-photo    Solve terrestrial position from a real photo (uses: <gha_aries> "
            "<g_x> <g_y> <g_z> [fov_deg] [observation_year] [crop_size])\n");
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
        double fov = 12.0, observation_year = current_year();
        if ((argc >= 3 && !numeric_argument(argv[2], parse_int(argv[2], 1, INT_MAX, &w))) ||
            (argc >= 4 && !numeric_argument(argv[3], parse_int(argv[3], 1, INT_MAX, &h))) ||
            (argc >= 5 && !numeric_argument(argv[4], parse_double(argv[4], &fov))) ||
            (argc >= 6 && !numeric_argument(argv[5], parse_double(argv[5], &observation_year))))
            return 1;
        return do_identify(w, h, fov, observation_year);
    } else if (strcmp(argv[1], "--solve-photo") == 0 || strcmp(argv[1], "--fit-photo") == 0) {
        int attitude_only = strcmp(argv[1], "--fit-photo") == 0;
        int first = attitude_only ? 2 : 3;
        if (argc < first + 3 || argc > first + 6) {
            printf("Usage: %s %s %s<g_x> <g_y> <g_z> [fov_deg] [observation_year] [crop_size]\n",
                   argv[0], argv[1], attitude_only ? "" : "<gha_aries> ");
            return 1;
        }
        double year = current_year();
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
        return do_solve_photo(&options, year, attitude_only);
    } else {
        printf("Invalid option: %s\n", argv[1]);
        return 1;
    }
}
