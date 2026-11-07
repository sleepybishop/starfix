#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "starfix_arena.h"
#include "starfix_attitude.h"
#include "starfix_cli_io.h"
#include "starfix_identify.h"
#include "starfix_status.h"

static uint8_t mempool[10 * 1024 * 1024];
static starfix_arena_t arena = {mempool, mempool + sizeof(mempool)};

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int test_catalog_record_formats(void) {
    const char* path = "/tmp/starfix_coordinate_formats.bin";
    unsigned char bytes[40] = {
        [0] = 'S', [1] = 'F',   [2] = 'I',   [3] = 'X',
        [4] = 1,   [8] = 1,     [12] = 50, /* one star, one hash entry */
        [16] = 42, [21] = 0x80, [23] = 0x80, [24] = 128,
    };
    const unsigned char formats[] = {0, 1, 1, 2, 0};
    const unsigned char ra_lows[] = {0, 0, 52, 0, 52};
    const unsigned char dec_lows[] = {0, 0, 86, 0, 86};
    int passed = 1;
    for (size_t i = 0; i < sizeof(formats); i++) {
        bytes[25] = ra_lows[i];
        bytes[26] = dec_lows[i];
        bytes[27] = formats[i];
        FILE* file = fopen(path, "wb");
        if (!file) return 0;
        int written = fwrite(bytes, 1, sizeof(bytes), file) == sizeof(bytes);
        if (fclose(file) != 0) written = 0;
        if (!written) {
            remove(path);
            return 0;
        }
        starfix_db_t loaded = {0};
        int status = starfix_load_db(path, &loaded);
        if (i < 3) {
            double ra = i == 0 ? 32768.0 * (2.0 * M_PI) / 65535.0
                               : (8388608.0 + ra_lows[i]) * (2.0 * M_PI) / 16777215.0;
            double dec = i == 0 ? 32768.0 * M_PI / 65535.0 - M_PI / 2.0
                                : (8388608.0 + dec_lows[i]) * M_PI / 16777215.0 - M_PI / 2.0;
            if (status != STARFIX_SUCCESS || loaded.num_stars != 1 ||
                fabs(loaded.stars[0].ra0 - ra) > 1e-14 || fabs(loaded.stars[0].dec0 - dec) > 1e-14)
                passed = 0;
        } else if (status != STARFIX_ERR_DATA_CORRUPTION || loaded.stars || loaded.hashes ||
                   loaded.num_stars || loaded.num_hashes || loaded.bins) {
            passed = 0;
        }
        starfix_free_db(&loaded);
    }
    remove(path);
    return passed;
}

static int identify_synthetic_field(const double x_deg[5], const double y_deg[5],
                                    const double dx_deg[5], const double dy_deg[5]) {
    const double focal = starfix_focal_length(1024.0, 24.0);
    starfix_vector3_t vectors[5];
    starfix_catalog_star_t catalog[5] = {{0}};
    starfix_centroid_t centroids[5] = {{0}};
    for (int i = 0; i < 5; i++) {
        double x = x_deg[i] * M_PI / 180.0, y = y_deg[i] * M_PI / 180.0;
        vectors[i] = (starfix_vector3_t){sin(x) * cos(y), sin(y), cos(x) * cos(y)};
        catalog[i].ra = atan2(vectors[i].y, vectors[i].x);
        catalog[i].dec = asin(vectors[i].z);
        /* Reflection across an almost straight star pattern is a small perturbation;
         * for a genuinely two-dimensional field it cannot fit a proper rotation. */
        double observed_x = x + dx_deg[i] * M_PI / 180.0;
        double observed_y = -y + dy_deg[i] * M_PI / 180.0;
        double z = cos(observed_x) * cos(observed_y);
        centroids[i].u = 512.0 + focal * sin(observed_x) * cos(observed_y) / z;
        centroids[i].v = 512.0 + focal * sin(observed_y) / z;
        centroids[i].flux = 100.0 - i;
    }

    starfix_hash_entry_t hashes[5] = {{0}};
    for (int omitted = 0; omitted < 5; omitted++) {
        starfix_hash_entry_t* entry = &hashes[omitted];
        int at = 0;
        for (int i = 0; i < 5; i++)
            if (i != omitted) entry->stars[at++] = (uint16_t)i;
        double angles[6];
        at = 0;
        for (int a = 0; a < 3; a++) {
            for (int b = a + 1; b < 4; b++) {
                double dot = starfix_dot(vectors[entry->stars[a]], vectors[entry->stars[b]]);
                angles[at++] = acos(fmax(-1.0, fmin(1.0, dot)));
            }
        }
        for (int a = 1; a < 6; a++) {
            for (int b = a; b > 0 && angles[b] < angles[b - 1]; b--) {
                double swap = angles[b];
                angles[b] = angles[b - 1];
                angles[b - 1] = swap;
            }
        }
        for (int i = 0; i < 5; i++)
            entry->hash_key |= (uint32_t)(angles[i] / angles[5] * 63.0) << (6 * i);
    }
    for (int a = 1; a < 5; a++) {
        for (int b = a; b > 0 && hashes[b].hash_key < hashes[b - 1].hash_key; b--) {
            starfix_hash_entry_t swap = hashes[b];
            hashes[b] = hashes[b - 1];
            hashes[b - 1] = swap;
        }
    }
    starfix_db_t db = {catalog, hashes, 5, 5, 63};
    starfix_match_t matches[5];
    int count = starfix_identify_stars(5, centroids, 1024, 1024, 24, &db, 5, matches, arena, NULL);
    for (int i = 0; i < count; i++)
        if (matches[i].centroid_idx != matches[i].catalog_idx) return 0;
    return count;
}

int main(void) {
    int failures = 0;
    printf("--- Running Star Identification (TETRA) Comprehensive Tests ---\n");

    /* 1. NULL parameter validations */
    printf("Test 1: Sizing and NULL validations... ");
    starfix_match_t matches[10];
    int ret1 = starfix_identify_stars(0, NULL, 1024, 1024, 12.0, NULL, 10, matches, arena, NULL);
    if (ret1 == STARFIX_ERR_NULL_POINTER) {
        printf("Passed\n");
    } else {
        failures++;
        printf("FAILED (got %d)\n", ret1);
    }

    /* 2. Database loading */
    printf("Test 2: Loading binary database... ");
    starfix_db_t db = {0};

    int db_status = starfix_load_db("data/starfix_db.bin", &db);
    if (db_status == 0 && db.num_stars > 0 && db.num_hashes > 0 && db.bins > 0) {
        printf("Passed\n");
        printf("  Loaded %d catalog stars and %d hash entries (bin factor: %d)\n", db.num_stars,
               db.num_hashes, db.bins);
    } else {
        failures++;
        printf("FAILED (load status: %d, stars: %d, entries: %d)\n", db_status, db.num_stars,
               db.num_hashes);
        return -1;
    }

    /* 2b. corrupt database indices must be rejected before catalog access */
    printf("Test 2b: Rejecting corrupt database indices... ");
    const char* corrupt_path = "/tmp/starfix_corrupt_db.bin";
    unsigned char corrupt_db[40] = {
        'S', 'F', 'I', 'X', 1, 0, 0, 0,             /* one catalog star */
        1,   0,   0,   0,                           /* one hash entry */
        50,  0,   0,   0,                           /* bin factor */
        0,   0,   0,   0,   0, 0, 0, 0, 0, 0, 0, 0, /* catalog record */
        0,   0,   0,   0,   1, 0, 0, 0, 0, 0, 0, 0  /* invalid star index 1 */
    };
    FILE* corrupt_file = fopen(corrupt_path, "wb");
    starfix_db_t corrupt = {.num_stars = 99, .num_hashes = 99, .bins = 99};
    int corrupt_status = STARFIX_ERR_MEMORY;
    if (corrupt_file != NULL) {
        if (fwrite(corrupt_db, 1, sizeof(corrupt_db), corrupt_file) == sizeof(corrupt_db)) {
            fclose(corrupt_file);
            corrupt_file = NULL;
            corrupt_status = starfix_load_db(corrupt_path, &corrupt);
        }
        if (corrupt_file != NULL) fclose(corrupt_file);
        remove(corrupt_path);
    }
    if (corrupt_status == STARFIX_ERR_DATA_CORRUPTION && corrupt.stars == NULL &&
        corrupt.hashes == NULL && corrupt.num_stars == 0 && corrupt.num_hashes == 0 &&
        corrupt.bins == 0) {
        printf("Passed\n");
    } else {
        failures++;
        printf("FAILED (status=%d)\n", corrupt_status);
        starfix_free_db(&corrupt);
    }

    printf("Test 2c: Tagged 24-bit and legacy 16-bit coordinates... ");
    if (test_catalog_record_formats()) {
        printf("Passed\n");
    } else {
        failures++;
        printf("FAILED\n");
    }

    printf("Test 2d: Noisy near-collinear seeds and mirrored fields... ");
    const double near_x[5] = {-10, -4.5, 0.3, 5.5, 10};
    const double near_y[5] = {0.04, -0.04, 0.03, -0.03, 0.02};
    const double dx[5] = {0.003, -0.002, 0.004, -0.004, 0.001};
    const double dy[5] = {0.004, -0.006, 0.007, 0.001, -0.003};
    const double wide_x[5] = {-8, -3, 1, 5, 9};
    const double wide_y[5] = {-5, 6, -3, 4, -1};
    int near_status = identify_synthetic_field(near_x, near_y, dx, dy);
    int mirror_status = identify_synthetic_field(wide_x, wide_y, dx, dy);
    if (near_status == 5 && mirror_status == STARFIX_ERR_IDENTIFY_NO_MATCH) {
        printf("Passed\n");
    } else {
        failures++;
        printf("FAILED (near-collinear=%d, mirrored=%d)\n", near_status, mirror_status);
    }

    /* 3. Lost-in-Space Star Identification */
    printf("Test 3: Lost-in-Space Star Identification... ");

    starfix_centroid_t centroids[50];
    int count = read_centroids_json("data/mock_centroids.json", centroids, 50);
    if (count < 0) {
        failures++;
        printf("FAILED (cannot read data/mock_centroids.json)\n");
        starfix_free_db(&db);
        return -1;
    }
    for (int i = 0; i < count; i++) centroids[i].flux = 1000.0 - i * 10.0;

    printf("  Parsed %d mock centroids\n", count);
    int di;
    for (di = 0; di < (count < 5 ? count : 5); di++) {
        printf("    c[%d]: u=%f v=%f flux=%f\n", di, centroids[di].u, centroids[di].v,
               centroids[di].flux);
    }

    if (count < 4) {
        failures++;
        printf("FAILED (too few mock centroids parsed: %d)\n", count);
        starfix_free_db(&db);
        return -1;
    }

    /* Run identification */
    starfix_match_t out_matches[20];
    int ret = starfix_identify_stars(count, centroids, 1024, 1024, 12.0, &db, 20, out_matches,
                                     arena, NULL);

    if (ret >= 3) {
        int lis_unique = 1;
        int unique_i, unique_j;
        for (unique_i = 0; unique_i < ret; unique_i++) {
            for (unique_j = unique_i + 1; unique_j < ret; unique_j++) {
                if (out_matches[unique_i].centroid_idx == out_matches[unique_j].centroid_idx ||
                    out_matches[unique_i].catalog_idx == out_matches[unique_j].catalog_idx) {
                    lis_unique = 0;
                }
            }
        }
        printf(lis_unique ? "Passed\n" : "FAILED (lost-in-space matches were not unique)\n");
        if (!lis_unique) return -1;
        printf("  Successfully matched %d stars in C!\n", ret);
        int i;
        for (i = 0; i < ret; i++) {
            int c_idx = out_matches[i].centroid_idx;
            int cat_idx = out_matches[i].catalog_idx;
            printf("    Centroid %d matches Catalog Star Index %d (HIP: %d, Name: %s)\n", c_idx,
                   cat_idx, db.stars[cat_idx].hip,
                   (db.stars[cat_idx].hip == 26241) ? "Hatysa" : "Unknown");
        }
        printf("Test 3c: Blind fallback and output bounds... ");
        const double wrong_hint[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
        starfix_match_t limited[4] = {{-1, -1}, {-1, -1}, {-1, -1}, {-1, -1}};
        int fallback_status = starfix_identify_stars(count, centroids, 1024, 1024, 12.0, &db, 3,
                                                     limited, arena, wrong_hint);
        int fallback_ok = fallback_status == 3 &&
                          memcmp(limited, out_matches, 3 * sizeof(*limited)) == 0 &&
                          limited[3].centroid_idx == -1 && limited[3].catalog_idx == -1;
        int quad_status =
            starfix_identify_stars(4, centroids, 1024, 1024, 12.0, &db, 3, limited, arena, NULL);
        if (!fallback_ok || quad_status != STARFIX_ERR_IDENTIFY_NO_MATCH) {
            failures++;
            printf("FAILED (fallback=%d, quad=%d)\n", fallback_status, quad_status);
            return -1;
        }
        printf("Passed\n");

        /* 3b. Tracking Mode Check (with attitude hint) */
        printf("Test 3b: Tracking Mode with Attitude Hint... ");
        starfix_vector3_t w_vecs[20], v_vecs[20];
        int m_cnt = ret;
        for (i = 0; i < m_cnt; i++) {
            int c_i = out_matches[i].centroid_idx;
            int cat_i = out_matches[i].catalog_idx;
            double f_len = 1024.0 / (2.0 * tan((12.0 * M_PI / 180.0) / 2.0));
            double cx = (centroids[c_i].u - 512.0) / f_len;
            double cy = (centroids[c_i].v - 512.0) / f_len;
            double norm = sqrt(cx * cx + cy * cy + 1.0);
            w_vecs[i].x = cx / norm;
            w_vecs[i].y = cy / norm;
            w_vecs[i].z = 1.0 / norm;

            double dec_s = db.stars[cat_i].dec;
            double ra_s = db.stars[cat_i].ra;
            v_vecs[i].x = cos(dec_s) * cos(ra_s);
            v_vecs[i].y = cos(dec_s) * sin(ra_s);
            v_vecs[i].z = sin(dec_s);
        }
        double R_hint[3][3];
        starfix_solve_attitude(m_cnt, w_vecs, v_vecs, R_hint);

        starfix_match_t track_matches[20];
        int track_ret = starfix_identify_stars(count, centroids, 1024, 1024, 12.0, &db, 20,
                                               track_matches, arena, (const double (*)[3])R_hint);
        if (track_ret >= 3) {
            int unique = 1;
            int a, b;
            for (a = 0; a < track_ret; a++) {
                for (b = a + 1; b < track_ret; b++) {
                    if (track_matches[a].centroid_idx == track_matches[b].centroid_idx) unique = 0;
                }
            }
            if (unique) {
                printf("Passed (%d unique matches in tracking mode!)\n", track_ret);
            } else {
                failures++;
                printf("FAILED (tracking mode reused a centroid)\n");
                return -1;
            }
        } else {
            failures++;
            printf("FAILED (tracking mode matched %d stars)\n", track_ret);
            return -1;
        }
    } else {
        failures++;
        printf("FAILED (matched only %d stars)\n", ret);
        return -1;
    }

    /* 4. Stellar Precession Propagation Check */
    printf("Test 4: Stellar Precession Propagation... ");
    starfix_catalog_star_t dummy_cat[1] = {{0}};
    dummy_cat[0].hip = 12345;
    dummy_cat[0].ra0 = 0.0;
    dummy_cat[0].dec0 = 0.0;
    dummy_cat[0].ra = 0.0;
    dummy_cat[0].dec = 0.0;
    dummy_cat[0].mag = 2.0;

    int prec_status = starfix_propagate_astrometry(dummy_cat, 1, 2100.0);

    /* Expected offsets using exact IAU 1976 Precession Matrix:
       T = 1.0 */
    double expected_d_dec = 0.009755;
    double expected_d_ra = 0.022339;

    if (prec_status == 0 && fabs(dummy_cat[0].dec - expected_d_dec) < 1e-4 &&
        fabs(dummy_cat[0].ra - expected_d_ra) < 1e-4) {
        printf("Passed\n");
    } else {
        failures++;
        printf("FAILED (status=%d, ra=%f (exp=%f), dec=%f (exp=%f))\n", prec_status,
               dummy_cat[0].ra, expected_d_ra, dummy_cat[0].dec, expected_d_dec);
    }

    starfix_free_db(&db);
    starfix_free_db(&db);
    if (db.stars || db.hashes || db.num_stars || db.num_hashes || db.bins ||
        starfix_load_db(NULL, &db) != STARFIX_ERR_NULL_POINTER ||
        starfix_load_db("data/starfix_db.bin", &db) != STARFIX_SUCCESS) {
        failures++;
        printf("FAILED (database cleanup and reload)\n");
        return -1;
    }
    starfix_free_db(&db);
    return failures ? 1 : 0;
}
