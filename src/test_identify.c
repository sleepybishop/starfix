#include <math.h>
#include "starfix_status.h"
#include <stdio.h>
#include <stdlib.h>

#include "starfix_arena.h"

static uint8_t mempool[10 * 1024 * 1024];
static starfix_arena_t arena;
#define RESET_ARENA() \
    (arena.beg = mempool, arena.end = mempool + sizeof(mempool))
#include <string.h>

#include "starfix_identify.h"

/* test runner to verify the C implementation of the TETRA Star ID module */

int main() {
    starfix_telemetry_t telem = {0};
    (void)telem;
    printf("--- Running Star Identification (TETRA) Comprehensive Tests ---\n");

    /* 1. NULL parameter validations */
    printf("Test 1: Sizing and NULL validations... ");
    starfix_match_t matches[10];
    int ret1 = (RESET_ARENA(), starfix_identify_stars(0, NULL, 1024, 1024, 12.0, 0, NULL, 0, NULL,
                                                      50, 10, matches, &arena, NULL));
    if (ret1 == -1) {
        printf("Passed\n");
    } else {
        printf("FAILED (got %d)\n", ret1);
    }

    /* 2. Database loading */
    printf("Test 2: Loading binary database... ");
    starfix_catalog_star_t* catalog = NULL;
    starfix_hash_entry_t* hash_table = NULL;
    uint32_t num_stars = 0;
    uint32_t num_entries = 0;
    uint32_t bin_factor = 0;

    int db_status = starfix_load_db("data/starfix_db.bin", &catalog, &num_stars, &hash_table,
                                    &num_entries, &bin_factor);
    if (db_status == 0 && num_stars > 0 && num_entries > 0 && bin_factor > 0) {
        printf("Passed\n");
        printf("  Loaded %d catalog stars and %d hash entries (bin factor: %d)\n", num_stars,
               num_entries, bin_factor);
    } else {
        printf("FAILED (load status: %d, stars: %d, entries: %d)\n", db_status, num_stars,
               num_entries);
        return -1;
    }

    /* 3. Lost-in-Space Star Identification */
    printf("Test 3: Lost-in-Space Star Identification... ");

    /* Parse mock centroids from data/mock_centroids.json using simple string parsing */
    FILE* f = fopen("data/mock_centroids.json", "r");
    if (f == NULL) {
        printf("FAILED (cannot open data/mock_centroids.json)\n");
        free(catalog);
        free(hash_table);
        return -1;
    }

    fseek(f, 0, SEEK_END);
    long sz = ftell(f);
    fseek(f, 0, SEEK_SET);
    char* buf = (char*)malloc((size_t)sz + 1);
    if (buf == NULL) {
        fclose(f);
        free(catalog);
        free(hash_table);
        return -1;
    }
    if (fread(buf, 1, (size_t)sz, f) != (size_t)sz) {
        free(buf);
        fclose(f);
        free(catalog);
        free(hash_table);
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
        free(catalog);
        free(hash_table);
        return -1;
    }

    /* Run identification */
    starfix_match_t out_matches[20];
    starfix_status_t match_status =
        (RESET_ARENA(),
         starfix_identify_stars(count, centroids, 1024, 1024, 12.0, num_stars, catalog, num_entries,
                                hash_table, bin_factor, 20, out_matches, &arena, NULL));

    if (match_status == STARFIX_SUCCESS && telem.identify_matches >= 3) {
        printf("Passed\n");
        printf("  Successfully matched %d stars in C!\n", (int)telem.identify_matches);
        int i;
        for (i = 0; i < (int)telem.identify_matches; i++) {
            int c_idx = out_matches[i].centroid_idx;
            int cat_idx = out_matches[i].catalog_idx;
            printf("    Centroid %d matches Catalog Star Index %d (HIP: %d, Name: %s)\n", c_idx,
                   cat_idx, catalog[cat_idx].hip,
                   (catalog[cat_idx].hip == 26241) ? "Hatysa" : "Unknown");
        }
    } else {
        printf("FAILED (matched only %d stars)\n", (int)telem.identify_matches);
    }

    /* 4. Stellar Precession Propagation Check */
    printf("Test 4: Stellar Precession Propagation... ");
    starfix_catalog_star_t dummy_cat[1];
    dummy_cat[0].hip = 12345;
    starfix_catalog_set_ra(&dummy_cat[0], 0.0);
    starfix_catalog_set_dec(&dummy_cat[0], 0.0);
    starfix_catalog_set_mag(&dummy_cat[0], 2.0);

    int prec_status = starfix_propagate_precession(dummy_cat, 1, 2100.0);

    /* Expected offsets:
       T = 1.0
       d_dec = 20.0431 * PI / (3600.0 * 180.0) ~ 9.717e-5
       d_ra = 46.1244 * PI / (3600.0 * 180.0) ~ 2.236e-4 */
    double expected_d_dec = 20.0431 * M_PI / (3600.0 * 180.0);
    double expected_d_ra = 46.1244 * M_PI / (3600.0 * 180.0);

    if (prec_status == 0 && fabs(starfix_catalog_dec(&dummy_cat[0]) - expected_d_dec) < 1e-4 &&
        fabs(starfix_catalog_ra(&dummy_cat[0]) - expected_d_ra) < 1e-4) {
        printf("Passed\n");
    } else {
        printf("FAILED (status=%d, ra=%f (exp=%f), dec=%f (exp=%f))\n", prec_status,
               starfix_catalog_ra(&dummy_cat[0]), expected_d_ra, starfix_catalog_dec(&dummy_cat[0]),
               expected_d_dec);
    }

    free(catalog);
    free(hash_table);
    return 0;
}
