#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "starfix_arena.h"
#include "starfix_status.h"

static uint8_t mempool[10 * 1024 * 1024];
static starfix_arena_t arena = {mempool, mempool + sizeof(mempool)};
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "starfix_attitude.h"
#include "starfix_identify.h"

/* test runner to verify the C implementation of the TETRA Star ID module */

int main() {
    printf("--- Running Star Identification (TETRA) Comprehensive Tests ---\n");

    /* 1. NULL parameter validations */
    printf("Test 1: Sizing and NULL validations... ");
    starfix_match_t matches[10];
    int ret1 = starfix_identify_stars(0, NULL, 1024, 1024, 12.0, NULL, 10, matches, arena, NULL);
    if (ret1 == -1) {
        printf("Passed\n");
    } else {
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
        printf("FAILED (load status: %d, stars: %d, entries: %d)\n", db_status, db.num_stars,
               db.num_hashes);
        return -1;
    }

    /* 2b. Corrupt database indices must be rejected before catalog access. */
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
        printf("FAILED (status=%d)\n", corrupt_status);
        starfix_free_db(&corrupt);
    }

    /* 3. Lost-in-Space Star Identification */
    printf("Test 3: Lost-in-Space Star Identification... ");

    /* Parse mock centroids from data/mock_centroids.json using simple string parsing */
    FILE* f = fopen("data/mock_centroids.json", "r");
    if (f == NULL) {
        printf("FAILED (cannot open data/mock_centroids.json)\n");
        starfix_free_db(&db);
        return -1;
    }

    fseek(f, 0, SEEK_END);
    long sz = ftell(f);
    fseek(f, 0, SEEK_SET);
    char* buf = (char*)malloc((size_t)sz + 1);
    if (buf == NULL) {
        fclose(f);
        starfix_free_db(&db);
        return -1;
    }
    if (fread(buf, 1, (size_t)sz, f) != (size_t)sz) {
        free(buf);
        fclose(f);
        starfix_free_db(&db);
        return -1;
    }
    buf[sz] = '\0';
    fclose(f);

    starfix_centroid_t centroids[50];
    int count = 0;
    char* obj_start = buf;
    while (count < 50) {
        obj_start = strchr(obj_start, '{');
        if (!obj_start) break;
        char* obj_end = strchr(obj_start, '}');
        if (!obj_end) break;

        double u = -1.0;
        double v = -1.0;

        char* u_ptr = strstr(obj_start, "\"u_true\"");
        if (u_ptr && u_ptr < obj_end) {
            char* p_val = strchr(u_ptr, ':');
            if (p_val) u = atof(p_val + 1);
        }

        char* v_ptr = strstr(obj_start, "\"v_true\"");
        if (v_ptr && v_ptr < obj_end) {
            char* p_val = strchr(v_ptr, ':');
            if (p_val) v = atof(p_val + 1);
        }

        if (u >= 0.0 && v >= 0.0) {
            centroids[count].u = u;
            centroids[count].v = v;
            centroids[count].flux = 1000.0 - count * 10.0;
            count++;
        }
        obj_start = obj_end + 1;
    }
    free(buf);

    printf("  Parsed %d mock centroids\n", count);
    int di;
    for (di = 0; di < (count < 5 ? count : 5); di++) {
        printf("    c[%d]: u=%f v=%f flux=%f\n", di, centroids[di].u, centroids[di].v,
               centroids[di].flux);
    }

    if (count < 4) {
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
        starfix_quaternion_t q_sol;
        double R_hint[3][3];
        starfix_solve_attitude(m_cnt, w_vecs, v_vecs, &q_sol, R_hint);

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
                printf("FAILED (tracking mode reused a centroid)\n");
                return -1;
            }
        } else {
            printf("FAILED (tracking mode matched %d stars)\n", track_ret);
            return -1;
        }
    } else {
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

    int prec_status = starfix_propagate_astrometry(dummy_cat, 1, 2100.0, NULL);

    /* Expected offsets using exact IAU 1976 Precession Matrix:
       T = 1.0 */
    double expected_d_dec = 0.009755;
    double expected_d_ra = 0.022339;

    if (prec_status == 0 && fabs(dummy_cat[0].dec - expected_d_dec) < 1e-4 &&
        fabs(dummy_cat[0].ra - expected_d_ra) < 1e-4) {
        printf("Passed\n");
    } else {
        printf("FAILED (status=%d, ra=%f (exp=%f), dec=%f (exp=%f))\n", prec_status,
               dummy_cat[0].ra, expected_d_ra, dummy_cat[0].dec, expected_d_dec);
    }

    starfix_free_db(&db);
    starfix_free_db(&db);
    if (db.stars || db.hashes || db.num_stars || db.num_hashes || db.bins ||
        starfix_load_db(NULL, &db) != STARFIX_ERR_NULL_POINTER ||
        starfix_load_db("data/starfix_db.bin", &db) != STARFIX_SUCCESS) {
        printf("FAILED (database cleanup and reload)\n");
        return -1;
    }
    starfix_free_db(&db);
    return 0;
}
