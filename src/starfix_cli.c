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

static double current_year(void) {
    time_t now = time(NULL);
    struct tm* utc = gmtime(&now);
    return utc ? 1900.0 + utc->tm_year + utc->tm_yday / 365.25 : NAN;
}

static int load_catalog(starfix_db_t* db, double year) {
    if (!isfinite(year)) return STARFIX_ERR_INVALID_PARAM;
    if (starfix_load_db("data/starfix_db.bin", db) != 0 &&
        starfix_load_db("../data/starfix_db.bin", db) != 0) {
        fprintf(stderr, "Could not load data/starfix_db.bin or ../data/starfix_db.bin.\n");
        return STARFIX_ERR_DATA_CORRUPTION;
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
    write_centroids_json("data/detected_centroids.json", centroids, count);
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

    write_identified_json("data/identified_stars.json", matches, match_count, db.stars);
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
    if (cent_count <= 0) {
        cent_count = read_centroids_json("../data/detected_centroids.json", centroids, 1000);
    }
    if (cent_count <= 0) {
        printf("Error: Could not load centroids from data/detected_centroids.json\n");
        goto done;
    }

    starfix_match_t matches[1000];
    int match_count = read_identified_json("data/identified_stars.json", matches, 1000);
    if (match_count <= 0) {
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
        double observation_year = argc >= 6 ? atof(argv[5]) : current_year();
        return do_identify(w, h, fov, observation_year);
    } else if (strcmp(argv[1], "--solve-photo") == 0 || strcmp(argv[1], "--fit-photo") == 0) {
        int attitude_only = strcmp(argv[1], "--fit-photo") == 0;
        int first = attitude_only ? 2 : 3;
        if (argc < first + 3 || argc > first + 6) {
            printf("Usage: %s %s %s<g_x> <g_y> <g_z> [fov_deg] [observation_year] [crop_size]\n",
                   argv[0], argv[1], attitude_only ? "" : "<gha_aries> ");
            return 1;
        }
        double year = argc > first + 4 ? atof(argv[first + 4]) : current_year();
        double size = argc > first + 5 ? atof(argv[first + 5]) : 560.0;
        starfix_photo_options_t options = {
            .width = size,
            .height = size,
            .fov_deg = argc > first + 3 ? atof(argv[first + 3]) : 11.28,
            .gha_aries = attitude_only ? 0.0 : atof(argv[2]),
            .zenith_camera = {atof(argv[first]), atof(argv[first + 1]), atof(argv[first + 2])},
        };
        return do_solve_photo(&options, year, attitude_only);
    } else {
        printf("Invalid option: %s\n", argv[1]);
        return 1;
    }
}
