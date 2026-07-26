#include "starfix_identify.h"

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "starfix_status.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* helper to compress 5D hash code into a single 30-bit integer */
static uint32_t compress_hash_key(const int* code) {
    uint32_t key_val = 0;
    int i;
    for (i = 0; i < 5; i++) {
        int val = code[i];
        if (val < 0) val = 0;
        if (val > 63) val = 63;
        key_val += (uint32_t)val * (1U << (6 * i));
    }
    return key_val;
}

int starfix_load_db(const char* db_path, starfix_catalog_star_t** catalog_stars,
                    uint32_t* num_stars, starfix_hash_entry_t** hash_entries, uint32_t* num_entries,
                    uint32_t* bin_factor) {
    if (db_path == NULL || catalog_stars == NULL || num_stars == NULL || hash_entries == NULL ||
        num_entries == NULL || bin_factor == NULL) {
        return STARFIX_ERR_NULL_POINTER;
    }

    FILE* f = fopen(db_path, "rb");
    if (f == NULL) {
        return STARFIX_ERR_MEMORY;
    }

    char magic[4];
    if (fread(magic, 1, 4, f) != 4 || memcmp(magic, "SFIX", 4) != 0) {
        fclose(f);
        return -3;
    }

    uint32_t n_stars, n_entries, b_factor;
    if (fread(&n_stars, 4, 1, f) != 1 || fread(&n_entries, 4, 1, f) != 1 ||
        fread(&b_factor, 4, 1, f) != 1) {
        fclose(f);
        return -4;
    }

    *catalog_stars =
        (starfix_catalog_star_t*)malloc((size_t)n_stars * sizeof(starfix_catalog_star_t));
    *hash_entries = (starfix_hash_entry_t*)malloc((size_t)n_entries * sizeof(starfix_hash_entry_t));

    if (*catalog_stars == NULL || *hash_entries == NULL) {
        if (*catalog_stars != NULL) free(*catalog_stars);
        if (*hash_entries != NULL) free(*hash_entries);
        fclose(f);
        return -5;
    }

    char* cat_buf = (char*)malloc((size_t)n_stars * 12);
    if (cat_buf == NULL || fread(cat_buf, 12, n_stars, f) != n_stars ||
        fread(*hash_entries, sizeof(starfix_hash_entry_t), n_entries, f) != n_entries) {
        if (cat_buf != NULL) free(cat_buf);
        free(*catalog_stars);
        free(*hash_entries);
        *catalog_stars = NULL;
        *hash_entries = NULL;
        fclose(f);
        return -5;
    }

    uint32_t i;
    for (i = 0; i < n_stars; i++) {
        int32_t hip;
        uint16_t ra_q, dec_q;
        uint8_t mag_q;
        memcpy(&hip, cat_buf + i * 12, 4);
        memcpy(&ra_q, cat_buf + i * 12 + 4, 2);
        memcpy(&dec_q, cat_buf + i * 12 + 6, 2);
        memcpy(&mag_q, cat_buf + i * 12 + 8, 1);

        (*catalog_stars)[i].hip = hip;
        (*catalog_stars)[i].ra = ra_q * (2.0 * M_PI) / 65535.0;
        (*catalog_stars)[i].dec = dec_q * M_PI / 65535.0 - M_PI / 2.0;
        (*catalog_stars)[i].mag = mag_q * 10.0 / 255.0 - 2.0;
    }
    free(cat_buf);

    *num_stars = n_stars;
    *num_entries = n_entries;
    *bin_factor = b_factor;

    fclose(f);
    return 0;
}

/* binary search helper to find the first entry matching a hash key */
static int find_first_hash_entry(uint32_t key, const starfix_hash_entry_t* entries,
                                 uint32_t num_entries) {
    int low = 0;
    int high = (int)num_entries - 1;
    int result = -1;
    while (low <= high) {
        int mid = low + (high - low) / 2;
        if (entries[mid].hash_key == key) {
            result = mid;
            high = mid - 1; /* search left for first match */
        } else if (entries[mid].hash_key < key) {
            low = mid + 1;
        } else {
            high = mid - 1;
        }
    }
    return result;
}

/* sorting helper for angles */
static void sort_angles(double* angles, int* indices, int n) {
    int i, j;

    for (i = 0; i < n - 1; i++) {
        for (j = 0; j < n - i - 1; j++) {
            if (angles[j] > angles[j + 1]) {
                double temp_a = angles[j];
                angles[j] = angles[j + 1];
                angles[j + 1] = temp_a;
                int temp_i = indices[j];
                indices[j] = indices[j + 1];
                indices[j + 1] = temp_i;
            }
        }
    }
}

starfix_status_t starfix_identify_stars(
    int num_centroids, const starfix_centroid_t* centroids, double width, double height,
    double fov_deg, uint32_t num_stars, const starfix_catalog_star_t* catalog_stars,
    uint32_t num_entries, const starfix_hash_entry_t* hash_entries, uint32_t bin_factor,
    int max_matches, starfix_match_t* matches, starfix_arena_t* arena, starfix_telemetry_t* telem,
    const double (*attitude_hint)[3]) {
    int i, j;

    if (num_centroids < 3 || centroids == NULL || catalog_stars == NULL || hash_entries == NULL ||
        matches == NULL || max_matches <= 0) {
        return STARFIX_ERR_NULL_POINTER;
    }

    double fov_rad = fov_deg * M_PI / 180.0;
    double focal_length = width / (2.0 * tan(fov_rad / 2.0));

    /* convert up to 20 centroids for tracking mode */
    int max_track_centroids = (num_centroids < 20) ? num_centroids : 20;
    double all_cam_vecs[20][3];
    for (i = 0; i < max_track_centroids; i++) {
        double cx = (centroids[i].u - width / 2.0) / focal_length;
        double cy = -(centroids[i].v - height / 2.0) / focal_length;
        double cz = 1.0;
        double norm = sqrt(cx * cx + cy * cy + cz * cz);
        all_cam_vecs[i][0] = cx / norm;
        all_cam_vecs[i][1] = cy / norm;
        all_cam_vecs[i][2] = cz / norm;
    }

    if (attitude_hint != NULL) {
        /* TRACKING MODE */
        int num_matches = 0;
        double match_tol_rad = 0.5 * M_PI / 180.0;
        double cos_match_tol = cos(match_tol_rad);
        double cos_fov_limit = cos(fov_rad / 2.0 + 1.0 * M_PI / 180.0);

        for (j = 0; j < (int)num_stars && num_matches < max_matches; j++) {
            double dec_s = starfix_catalog_dec(&catalog_stars[j]);
            double ra_s = starfix_catalog_ra(&catalog_stars[j]);
            double v_icrs[3] = {cos(dec_s) * cos(ra_s), cos(dec_s) * sin(ra_s), sin(dec_s)};

            double v_cam[3];
            v_cam[0] = attitude_hint[0][0] * v_icrs[0] + attitude_hint[0][1] * v_icrs[1] +
                       attitude_hint[0][2] * v_icrs[2];
            v_cam[1] = attitude_hint[1][0] * v_icrs[0] + attitude_hint[1][1] * v_icrs[1] +
                       attitude_hint[1][2] * v_icrs[2];
            v_cam[2] = attitude_hint[2][0] * v_icrs[0] + attitude_hint[2][1] * v_icrs[1] +
                       attitude_hint[2][2] * v_icrs[2];

            if (v_cam[2] < cos_fov_limit) continue;

            int best_c = -1;
            double best_dot = -1.0;
            for (i = 0; i < max_track_centroids; i++) {
                double dot = v_cam[0] * all_cam_vecs[i][0] + v_cam[1] * all_cam_vecs[i][1] +
                             v_cam[2] * all_cam_vecs[i][2];
                if (dot > best_dot) {
                    best_dot = dot;
                    best_c = i;
                }
            }

            if (best_c >= 0 && best_dot >= cos_match_tol) {
                matches[num_matches].centroid_idx = best_c;
                matches[num_matches].catalog_idx = j;
                num_matches++;
            }
        }

        if (num_matches >= 3) {
            if (telem) {
                telem->identify_matches = (uint32_t)num_matches;
            }
            return STARFIX_SUCCESS;
        }
    }

    if (num_centroids < 4) return STARFIX_ERR_NULL_POINTER;

    /* 1. select up to 8 brightest centroids for LIS hashing */
    int n_test = (num_centroids < 8) ? num_centroids : 8;
    double cam_vecs[8][3];
    for (i = 0; i < n_test; i++) {
        cam_vecs[i][0] = all_cam_vecs[i][0];
        cam_vecs[i][1] = all_cam_vecs[i][1];
        cam_vecs[i][2] = all_cam_vecs[i][2];
    }

    /* 2. allocate voting table (size: n_test x num_stars) */
    int** votes = starfix_alloc_array(arena, int*, (size_t)n_test);
    if (votes == NULL) return STARFIX_ERR_MEMORY;
    for (i = 0; i < n_test; i++) {
        votes[i] = starfix_alloc_array(arena, int, (size_t)num_stars);
        if (votes[i] == NULL) {
            return STARFIX_ERR_MEMORY;
        }
    }

    /* all 24 permutations of matching [0, 1, 2, 3] camera stars to catalog stars */
    static const int perms[24][4] = {
        {0, 1, 2, 3}, {0, 1, 3, 2}, {0, 2, 1, 3}, {0, 2, 3, 1}, {0, 3, 1, 2}, {0, 3, 2, 1},
        {1, 0, 2, 3}, {1, 0, 3, 2}, {1, 2, 0, 3}, {1, 2, 3, 0}, {1, 3, 0, 2}, {1, 3, 2, 0},
        {2, 0, 1, 3}, {2, 0, 3, 1}, {2, 1, 0, 3}, {2, 1, 3, 0}, {2, 3, 0, 1}, {2, 3, 1, 0},
        {3, 0, 1, 2}, {3, 0, 2, 1}, {3, 1, 0, 2}, {3, 1, 2, 0}, {3, 2, 0, 1}, {3, 2, 1, 0}};

    double angle_tolerance = 0.5 * M_PI / 180.0; /* 0.5 degrees tolerance for real sensors */

    /* 3. generate all combinations of 4 centroids (quads) */
    int a, b, c, d;
    for (a = 0; a < n_test; a++) {
        for (b = a + 1; b < n_test; b++) {
            for (c = b + 1; c < n_test; c++) {
                for (d = c + 1; d < n_test; d++) {
                    int quad_idx[4] = {a, b, c, d};

                    /* calculate 6 pairwise angles in camera frame */
                    double cam_angles[6];
                    int cam_angle_pairs[6][2];
                    int pair_idx = 0;

                    int p1, p2;
                    for (p1 = 0; p1 < 4; p1++) {
                        for (p2 = p1 + 1; p2 < 4; p2++) {
                            int idx1 = quad_idx[p1];
                            int idx2 = quad_idx[p2];
                            double dot = cam_vecs[idx1][0] * cam_vecs[idx2][0] +
                                         cam_vecs[idx1][1] * cam_vecs[idx2][1] +
                                         cam_vecs[idx1][2] * cam_vecs[idx2][2];
                            if (dot > 1.0) dot = 1.0;
                            if (dot < -1.0) dot = -1.0;
                            cam_angles[pair_idx] = acos(dot);
                            cam_angle_pairs[pair_idx][0] = p1;
                            cam_angle_pairs[pair_idx][1] = p2;
                            pair_idx++;
                        }
                    }

                    /* sort angles */
                    int sorted_idx[6] = {0, 1, 2, 3, 4, 5};
                    sort_angles(cam_angles, sorted_idx, 6);

                    double largest_angle = cam_angles[5];
                    if (largest_angle < 1e-5) continue;

                    /* compute ratios and hash coordinates */
                    double ratios[5];
                    int hash_code[5];
                    for (i = 0; i < 5; i++) {
                        ratios[i] = cam_angles[i] / largest_angle;
                        hash_code[i] = (int)(ratios[i] * bin_factor);
                    }

                    /* generate list of 11 hash keys to check (exact + +/-1 adjacent bins) */
                    uint32_t search_keys[11];
                    search_keys[0] = compress_hash_key(hash_code);

                    int key_count = 1;
                    int dk, dim;
                    for (dk = -1; dk <= 1; dk += 2) {
                        for (dim = 0; dim < 5; dim++) {
                            int adj_code[5];
                            memcpy(adj_code, hash_code, 5 * sizeof(int));
                            adj_code[dim] += dk;
                            if (adj_code[dim] < 0) adj_code[dim] = 0;
                            if (adj_code[dim] > (int)bin_factor) adj_code[dim] = (int)bin_factor;

                            search_keys[key_count++] = compress_hash_key(adj_code);
                        }
                    }

                    /* check each search key in the sorted hash database */
                    int k;
                    for (k = 0; k < key_count; k++) {
                        /* remove duplicates among adjacent keys */
                        int dup = 0;
                        int prev_k;
                        for (prev_k = 0; prev_k < k; prev_k++) {
                            if (search_keys[prev_k] == search_keys[k]) {
                                dup = 1;
                                break;
                            }
                        }
                        if (dup) continue;

                        int entry_idx =
                            find_first_hash_entry(search_keys[k], hash_entries, num_entries);
                        if (entry_idx < 0) continue;

                        /* scan all matching database entries */
                        int hash_scans = 0;
                        while (entry_idx < (int)num_entries &&
                               hash_entries[entry_idx].hash_key == search_keys[k] &&
                               hash_scans < 100) {
                            hash_scans++;
                            const starfix_hash_entry_t* cand = &hash_entries[entry_idx];

                            /* extract candidate star unit vectors */
                            double cat_vecs[4][3];
                            for (i = 0; i < 4; i++) {
                                int c_idx = cand->stars[i];
                                double dec_s = starfix_catalog_dec(&catalog_stars[c_idx]);
                                double ra_s = starfix_catalog_ra(&catalog_stars[c_idx]);
                                cat_vecs[i][0] = cos(dec_s) * cos(ra_s);
                                cat_vecs[i][1] = cos(dec_s) * sin(ra_s);
                                cat_vecs[i][2] = sin(dec_s);
                            }

                            /* calculate 6 pairwise angles in catalog frame */
                            double cat_angles[4][4];
                            int c1, c2;
                            for (c1 = 0; c1 < 4; c1++) {
                                for (c2 = 0; c2 < 4; c2++) {
                                    if (c1 == c2) {
                                        cat_angles[c1][c2] = 0.0;
                                    } else {
                                        double dot = cat_vecs[c1][0] * cat_vecs[c2][0] +
                                                     cat_vecs[c1][1] * cat_vecs[c2][1] +
                                                     cat_vecs[c1][2] * cat_vecs[c2][2];
                                        if (dot > 1.0) dot = 1.0;
                                        if (dot < -1.0) dot = -1.0;
                                        cat_angles[c1][c2] = acos(dot);
                                    }
                                }
                            }

                            /* test all 24 matching permutations */
                            int p;
                            for (p = 0; p < 24; p++) {
                                int valid = 1;
                                int pair;
                                for (pair = 0; pair < 6; pair++) {
                                    /* get camera frame angle and pair mapping indices */
                                    double c_angle = cam_angles[pair];
                                    int u_p1 = cam_angle_pairs[sorted_idx[pair]][0];
                                    int u_p2 = cam_angle_pairs[sorted_idx[pair]][1];

                                    /* map camera frame indices to catalog indices under permutation
                                     * p */
                                    int cat_p1 = perms[p][u_p1];
                                    int cat_p2 = perms[p][u_p2];

                                    double c_cat_angle = cat_angles[cat_p1][cat_p2];

                                    if (fabs(c_angle - c_cat_angle) > angle_tolerance) {
                                        valid = 0;
                                        break;
                                    }
                                }

                                if (valid) {
                                    /* permutation matched! vote for the matched pairs */
                                    int match_i;
                                    for (match_i = 0; match_i < 4; match_i++) {
                                        int centroid_idx = quad_idx[match_i];
                                        int catalog_idx = cand->stars[perms[p][match_i]];
                                        if (catalog_idx >= 0 && (uint32_t)catalog_idx < num_stars) {
                                            votes[centroid_idx][catalog_idx]++;
                                        }
                                    }
                                }
                            }
                            entry_idx++;
                        }
                    }
                }
            }
        }
    }

    /* 4. extract final matched pairs based on vote thresholds */
    int num_matches = 0;
    for (i = 0; i < n_test; i++) {
        int max_v = 0;
        int best_cat_idx = -1;
        for (j = 0; j < (int)num_stars; j++) {
            if (votes[i][j] > max_v) {
                max_v = votes[i][j];
                best_cat_idx = j;
            }
        }

        /* standard threshold: require at least 4 votes for verification */
        if (max_v >= 4 && best_cat_idx >= 0) {
            if (num_matches < max_matches) {
                matches[num_matches].centroid_idx = i;
                matches[num_matches].catalog_idx = best_cat_idx;
                num_matches++;
            }
        }
    }

    if (telem) {
        telem->identify_matches = (uint32_t)num_matches;
    }
    return STARFIX_SUCCESS;
}

int starfix_propagate_precession(starfix_catalog_star_t* catalog, uint32_t num_stars, double year) {
    if (catalog == NULL || num_stars == 0) {
        return STARFIX_ERR_NULL_POINTER;
    }

    double T = (year - 2000.0) / 100.0; /* Julian centuries */

    /* IAU 1976 Precession angles */
    double zeta_rad = (2306.2181 * T) * M_PI / (3600.0 * 180.0);
    double z_rad = (2306.2181 * T) * M_PI / (3600.0 * 180.0);
    double theta_rad = (2004.3109 * T) * M_PI / (3600.0 * 180.0);

    double c_zeta = cos(zeta_rad), s_zeta = sin(zeta_rad);
    double c_z = cos(z_rad), s_z = sin(z_rad);
    double c_th = cos(theta_rad), s_th = sin(theta_rad);

    double P[3][3];
    P[0][0] = c_zeta * c_th * c_z - s_zeta * s_z;
    P[0][1] = -s_zeta * c_th * c_z - c_zeta * s_z;
    P[0][2] = -s_th * c_z;
    P[1][0] = c_zeta * c_th * s_z + s_zeta * c_z;
    P[1][1] = -s_zeta * c_th * s_z + c_zeta * c_z;
    P[1][2] = -s_th * s_z;
    P[2][0] = c_zeta * s_th;
    P[2][1] = -s_zeta * s_th;
    P[2][2] = c_th;

    uint32_t i;
    for (i = 0; i < num_stars; i++) {
        double ra0 = starfix_catalog_ra(&catalog[i]);
        double dec0 = starfix_catalog_dec(&catalog[i]);

        double v0[3];
        v0[0] = cos(dec0) * cos(ra0);
        v0[1] = cos(dec0) * sin(ra0);
        v0[2] = sin(dec0);

        double v1[3];
        v1[0] = P[0][0] * v0[0] + P[0][1] * v0[1] + P[0][2] * v0[2];
        v1[1] = P[1][0] * v0[0] + P[1][1] * v0[1] + P[1][2] * v0[2];
        v1[2] = P[2][0] * v0[0] + P[2][1] * v0[1] + P[2][2] * v0[2];

        double ra_new = atan2(v1[1], v1[0]);
        double dec_new = asin(v1[2]);

        /* Wrap RA to [0, 2pi] */
        while (ra_new < 0.0) ra_new += 2.0 * M_PI;
        while (ra_new >= 2.0 * M_PI) ra_new -= 2.0 * M_PI;

        starfix_catalog_set_ra(&catalog[i], ra_new);
        starfix_catalog_set_dec(&catalog[i], dec_new);
    }

    return 0;
}
