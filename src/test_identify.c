#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "starfix_attitude.h"
#include "starfix_cli_io.h"
#include "starfix_identify.h"
#include "starfix_status.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int test_catalog_record_proper_motion(void) {
    const char* path = "/tmp/starfix_pm_record.bin";
    /* 16-byte header + 18-byte star record + 12-byte hash entry = 46 bytes */
    unsigned char bytes[46] = {
        [0] = 'S',   [1] = 'F',   [2] = 'I',   [3] = 'X',   [4] = 1,     [8] = 1, [12] = 50,
        [16] = 42,   [23] = 0x80, [27] = 0x80, [28] = 0xC4, [29] = 0x09, /* mag = 2.500 */
        [30] = 0x88, [31] = 0x13, /* pm_ra_q = +5000 (+1000 mas/yr) */
        [32] = 0xF0, [33] = 0xD8, /* pm_dec_q = -10000 (-2000 mas/yr) */
    };
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
    double exp_ra = 2147483648.0 * (2.0 * M_PI) / 4294967295.0;
    double exp_dec = 2147483648.0 * M_PI / 4294967295.0 - M_PI / 2.0;
    int passed = status == STARFIX_SUCCESS && loaded.num_stars == 1 && loaded.stars[0].hip == 42 &&
                 fabs(loaded.stars[0].mag - 2.5) < 1e-12 &&
                 fabs(loaded.stars[0].ra0 - exp_ra) < 1e-14 &&
                 fabs(loaded.stars[0].dec0 - exp_dec) < 1e-14;
    if (passed) {
        starfix_catalog_star_t zero_pm = starfix_make_catalog_star(42, exp_ra, exp_dec, 2.5);
        starfix_propagate_astrometry(&zero_pm, 1, 2100.0);
        starfix_propagate_astrometry(loaded.stars, 1, 2100.0);
        double arcsec_rad = M_PI / (180.0 * 3600.0);
        passed &= fabs((loaded.stars[0].ra - zero_pm.ra) - 100.0 * arcsec_rad) < 2e-6;
        passed &= fabs((loaded.stars[0].dec - zero_pm.dec) + 200.0 * arcsec_rad) < 2e-6;
    }
    starfix_free_db(&loaded);

    /* truncated 45-byte file must fail with STARFIX_ERR_DATA_CORRUPTION */
    file = fopen(path, "wb");
    if (!file) return 0;
    written = fwrite(bytes, 1, sizeof(bytes) - 1, file) == sizeof(bytes) - 1;
    if (fclose(file) != 0) written = 0;
    if (written) {
        passed &= starfix_load_db(path, &loaded) == STARFIX_ERR_DATA_CORRUPTION && !loaded.stars;
    } else {
        passed = 0;
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
        catalog[i] = starfix_make_catalog_star(i + 1, atan2(vectors[i].y, vectors[i].x),
                                               asin(vectors[i].z), 2.0);
        /* reflection across an almost straight star pattern is a small perturbation;
         * for a genuinely two-dimensional field it cannot fit a proper rotation */
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
    int count = starfix_identify_stars(5, centroids, 1024, 1024, 24, &db, 5, matches, NULL);
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
    int ret1 = starfix_identify_stars(0, NULL, 1024, 1024, 12.0, NULL, 10, matches, NULL);
    if (ret1 == STARFIX_ERR_NULL_POINTER) {
        printf("Passed\n");
    } else {
        failures++;
        printf("FAILED (got %d)\n", ret1);
    }

    /* 2. database loading and dense asterism coverage */
    printf("Test 2: Loading binary database... ");
    starfix_db_t db = {0};

    int db_status = starfix_load_db("data/starfix_db.bin", &db);
    int belt_found = 0;
    for (uint32_t s_i = 0; db_status == 0 && s_i < db.num_stars; s_i++) {
        int32_t hip = db.stars[s_i].hip;
        if (hip == 26727 || hip == 26311 || hip == 25930) belt_found++;
    }
    if (db_status == 0 && db.num_stars > 0 && db.num_hashes > 0 && db.bins > 0 && belt_found == 3) {
        printf("Passed\n");
        printf("  Loaded %d catalog stars and %d hash entries (bin factor: %d)\n", db.num_stars,
               db.num_hashes, db.bins);
    } else {
        failures++;
        printf("FAILED (load status: %d, stars: %d, entries: %d, belt: %d)\n", db_status,
               db.num_stars, db.num_hashes, belt_found);
        return -1;
    }

    /* 2b. corrupt database indices must be rejected before catalog access */
    printf("Test 2b: Rejecting corrupt database indices... ");
    const char* corrupt_path = "/tmp/starfix_corrupt_db.bin";
    unsigned char corrupt_db[46] = {
        [0] = 'S', [1] = 'F', [2] = 'I', [3] = 'X',
        [4] = 1,   [8] = 1,   [12] = 50, /* one catalog star, one hash entry */
        [38] = 1,                        /* invalid star index 1 in hash entry */
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

    printf("Test 2c: 32-bit coordinates and proper motion record... ");
    if (test_catalog_record_proper_motion()) {
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

    /* run identification */
    starfix_match_t out_matches[20];
    int ret =
        starfix_identify_stars(count, centroids, 1024, 1024, 12.0, &db, 20, out_matches, NULL);

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
        int fallback_status =
            starfix_identify_stars(count, centroids, 1024, 1024, 12.0, &db, 3, limited, wrong_hint);
        int fallback_ok = fallback_status == 3 &&
                          memcmp(limited, out_matches, 3 * sizeof(*limited)) == 0 &&
                          limited[3].centroid_idx == -1 && limited[3].catalog_idx == -1;
        int quad_status =
            starfix_identify_stars(4, centroids, 1024, 1024, 12.0, &db, 3, limited, NULL);
        if (!fallback_ok || quad_status != STARFIX_ERR_IDENTIFY_NO_MATCH) {
            failures++;
            printf("FAILED (fallback=%d, quad=%d)\n", fallback_status, quad_status);
            return -1;
        }
        printf("Passed\n");

        /* 3b. tracking mode check (with attitude hint) */
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
            v_vecs[i] = db.stars[cat_i].unit;
        }
        double R_hint[3][3];
        starfix_solve_attitude(m_cnt, w_vecs, v_vecs, R_hint);

        starfix_match_t track_matches[20];
        int track_ret = starfix_identify_stars(count, centroids, 1024, 1024, 12.0, &db, 20,
                                               track_matches, (const double (*)[3])R_hint);
        if (track_ret >= 3) {
            int unique = 1;
            int a, b;
            for (a = 0; a < track_ret; a++) {
                for (b = a + 1; b < track_ret; b++) {
                    if (track_matches[a].centroid_idx == track_matches[b].centroid_idx ||
                        track_matches[a].catalog_idx == track_matches[b].catalog_idx)
                        unique = 0;
                }
            }
            if (unique) {
                printf("Passed (%d unique matches in tracking mode!)\n", track_ret);
            } else {
                failures++;
                printf("FAILED (tracking mode reused a centroid or catalog star)\n");
                return -1;
            }
        } else {
            failures++;
            printf("FAILED (tracking mode matched %d stars)\n", track_ret);
            return -1;
        }

        /* 3d. wide-aspect corner tracking mode beyond horizontal half-FOV */
        printf("Test 3d: Wide-aspect corner tracking... ");
        starfix_centroid_t corner_c[4];
        starfix_catalog_star_t corner_cat[4];
        const double cu[4] = {30.0, 1890.0, 30.0, 1890.0};
        const double cv[4] = {30.0, 30.0, 1050.0, 1050.0};
        double corner_focal = starfix_focal_length(1920.0, 12.0);
        for (int k = 0; k < 4; k++) {
            corner_c[k] = (starfix_centroid_t){.u = cu[k], .v = cv[k], .flux = 100.0};
            starfix_vector3_t ray = starfix_camera_ray(cu[k] - 960.0, cv[k] - 540.0, corner_focal);
            corner_cat[k] = starfix_make_catalog_star(k + 1, atan2(ray.y, ray.x), asin(ray.z), 2.0);
        }
        starfix_db_t corner_db = {corner_cat, NULL, 4, 0, 50};
        const double eye[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
        starfix_match_t corner_m[4];
        int corner_ret =
            starfix_identify_stars(4, corner_c, 1920, 1080, 12.0, &corner_db, 4, corner_m, eye);
        if (corner_ret == 4) {
            printf("Passed\n");
        } else {
            failures++;
            printf("FAILED (matched %d/4 corner stars)\n", corner_ret);
        }
    } else {
        failures++;
        printf("FAILED (matched only %d stars)\n", ret);
        return -1;
    }

    /* 4. stellar precession propagation check */
    printf("Test 4: Stellar Precession Propagation... ");
    starfix_catalog_star_t dummy_cat[1] = {starfix_make_catalog_star(12345, 0.0, 0.0, 2.0)};

    int prec_status = starfix_propagate_astrometry(dummy_cat, 1, 2100.0);

    /* expected offsets using IAU 1976 precession + nutation + aberration at T = 1.0 */
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
