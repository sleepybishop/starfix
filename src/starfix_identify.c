#include "starfix_identify.h"

#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "starfix_status.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "starfix_astrometry.h"

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

static uint16_t decode_le16(const unsigned char* p) {
    return (uint16_t)((uint16_t)p[0] | ((uint16_t)p[1] << 8));
}

static uint32_t decode_le32(const unsigned char* p) {
    return (uint32_t)p[0] | ((uint32_t)p[1] << 8) | ((uint32_t)p[2] << 16) | ((uint32_t)p[3] << 24);
}

int starfix_load_db(const char* db_path, starfix_catalog_star_t** catalog_stars,
                    uint32_t* num_stars, starfix_hash_entry_t** hash_entries, uint32_t* num_entries,
                    uint32_t* bin_factor) {
    if (db_path == NULL || catalog_stars == NULL || num_stars == NULL || hash_entries == NULL ||
        num_entries == NULL || bin_factor == NULL) {
        return STARFIX_ERR_NULL_POINTER;
    }

    *catalog_stars = NULL;
    *hash_entries = NULL;
    *num_stars = 0;
    *num_entries = 0;
    *bin_factor = 0;

    FILE* f = fopen(db_path, "rb");
    if (f == NULL) {
        return STARFIX_ERR_MEMORY;
    }

    if (fseek(f, 0, SEEK_END) != 0) {
        fclose(f);
        return STARFIX_ERR_DATA_CORRUPTION;
    }
    long file_size = ftell(f);
    if (file_size < 16 || fseek(f, 0, SEEK_SET) != 0) {
        fclose(f);
        return STARFIX_ERR_DATA_CORRUPTION;
    }

    unsigned char header[16];
    if (fread(header, 1, sizeof(header), f) != sizeof(header) || memcmp(header, "SFIX", 4) != 0) {
        fclose(f);
        return STARFIX_ERR_DATA_CORRUPTION;
    }

    uint32_t n_stars = decode_le32(header + 4);
    uint32_t n_entries = decode_le32(header + 8);
    uint32_t b_factor = decode_le32(header + 12);
    uint64_t expected_size = 16ULL + (uint64_t)n_stars * 12ULL + (uint64_t)n_entries * 12ULL;
    if (n_stars == 0 || n_stars > 65536U || n_entries == 0 || n_entries > (uint32_t)INT_MAX ||
        b_factor == 0 || b_factor > 63U || expected_size != (uint64_t)file_size) {
        fclose(f);
        return STARFIX_ERR_DATA_CORRUPTION;
    }

    *catalog_stars =
        (starfix_catalog_star_t*)malloc((size_t)n_stars * sizeof(starfix_catalog_star_t));
    *hash_entries = (starfix_hash_entry_t*)malloc((size_t)n_entries * sizeof(starfix_hash_entry_t));

    if (*catalog_stars == NULL || *hash_entries == NULL) {
        if (*catalog_stars != NULL) free(*catalog_stars);
        if (*hash_entries != NULL) free(*hash_entries);
        *catalog_stars = NULL;
        *hash_entries = NULL;
        fclose(f);
        return STARFIX_ERR_MEMORY;
    }

    uint32_t i;
    for (i = 0; i < n_stars; i++) {
        unsigned char record[12];
        if (fread(record, 1, sizeof(record), f) != sizeof(record)) goto corrupt;
        int32_t hip;
        uint32_t hip_bits = decode_le32(record);
        uint16_t ra_q = decode_le16(record + 4);
        uint16_t dec_q = decode_le16(record + 6);
        uint8_t mag_q = record[8];
        memcpy(&hip, &hip_bits, sizeof(hip));

        (*catalog_stars)[i].hip = hip;
        (*catalog_stars)[i].ra0 = ra_q * (2.0 * M_PI) / 65535.0;
        (*catalog_stars)[i].dec0 = dec_q * M_PI / 65535.0 - M_PI / 2.0;
        (*catalog_stars)[i].ra = (*catalog_stars)[i].ra0;
        (*catalog_stars)[i].dec = (*catalog_stars)[i].dec0;
        (*catalog_stars)[i].mag = mag_q * 10.0 / 255.0 - 2.0;
    }

    uint32_t previous_key = 0;
    for (i = 0; i < n_entries; i++) {
        unsigned char record[12];
        if (fread(record, 1, sizeof(record), f) != sizeof(record)) goto corrupt;
        starfix_hash_entry_t* entry = &(*hash_entries)[i];
        entry->hash_key = decode_le32(record);
        int j;
        for (j = 0; j < 4; j++) {
            entry->stars[j] = decode_le16(record + 4 + 2 * j);
            if ((uint32_t)entry->stars[j] >= n_stars) goto corrupt;
        }
        if (i > 0 && entry->hash_key < previous_key) goto corrupt;
        previous_key = entry->hash_key;
    }

    *num_stars = n_stars;
    *num_entries = n_entries;
    *bin_factor = b_factor;

    fclose(f);
    return 0;

corrupt:
    free(*catalog_stars);
    free(*hash_entries);
    *catalog_stars = NULL;
    *hash_entries = NULL;
    fclose(f);
    return STARFIX_ERR_DATA_CORRUPTION;
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

    if (num_centroids < 3 || centroids == NULL || catalog_stars == NULL || matches == NULL) {
        return STARFIX_ERR_NULL_POINTER;
    }
    if (width <= 0.0 || height <= 0.0 || fov_deg <= 0.0 || fov_deg >= 180.0 || !isfinite(width) ||
        !isfinite(height) || !isfinite(fov_deg) || num_stars == 0 ||
        num_stars > (uint32_t)INT_MAX || max_matches <= 0) {
        return STARFIX_ERR_INVALID_PARAM;
    }

    double fov_rad = fov_deg * M_PI / 180.0;
    double focal_length = width / (2.0 * tan(fov_rad / 2.0));

    /* convert up to 20 centroids for tracking mode */
    int max_track_centroids = (num_centroids < 20) ? num_centroids : 20;
    double all_cam_vecs[20][3];
    for (i = 0; i < max_track_centroids; i++) {
        double cx = (centroids[i].u - width / 2.0) / focal_length;
        double cy = (centroids[i].v - height / 2.0) / focal_length;
        double cz = 1.0;
        double norm = sqrt(cx * cx + cy * cy + cz * cz);
        all_cam_vecs[i][0] = cx / norm;
        /* Image rows increase downward; camera-frame +Y increases upward. */
        all_cam_vecs[i][1] = -cy / norm;
        all_cam_vecs[i][2] = cz / norm;
    }

    if (attitude_hint != NULL) {
        /* TRACKING MODE */
        int num_matches = 0;
        int centroid_used[20] = {0};
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
                if (centroid_used[best_c]) continue;
                matches[num_matches].centroid_idx = best_c;
                matches[num_matches].catalog_idx = j;
                centroid_used[best_c] = 1;
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
    if (hash_entries == NULL || num_entries == 0 || bin_factor == 0 || bin_factor > 63U ||
        arena == NULL) {
        return STARFIX_ERR_INVALID_PARAM;
    }

    /* Include enough detections to tolerate bright hot pixels, merged stars, and edge artifacts. */
    int n_test = (num_centroids < 12) ? num_centroids : 12;
    double cam_vecs[12][3];
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
                            int candidate_valid = 1;
                            for (i = 0; i < 4; i++) {
                                int c_idx = cand->stars[i];
                                if (c_idx < 0 || (uint32_t)c_idx >= num_stars) {
                                    candidate_valid = 0;
                                    break;
                                }
                                double dec_s = starfix_catalog_dec(&catalog_stars[c_idx]);
                                double ra_s = starfix_catalog_ra(&catalog_stars[c_idx]);
                                cat_vecs[i][0] = cos(dec_s) * cos(ra_s);
                                cat_vecs[i][1] = cos(dec_s) * sin(ra_s);
                                cat_vecs[i][2] = sin(dec_s);
                            }
                            if (!candidate_valid) {
                                entry_idx++;
                                continue;
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

    /* 4. Extract a one-to-one assignment from the vote table. A catalog star and an image
     * centroid may each appear at most once in a physically valid attitude solution. */
    int num_matches = 0;
    int centroid_used[12] = {0};
    unsigned char* catalog_used = starfix_alloc_array(arena, unsigned char, (size_t)num_stars);
    if (catalog_used == NULL) return STARFIX_ERR_MEMORY;

    while (num_matches < n_test && num_matches < max_matches) {
        int best_votes = 3;
        int best_centroid = -1;
        int best_catalog = -1;
        for (i = 0; i < n_test; i++) {
            if (centroid_used[i]) continue;
            for (j = 0; j < (int)num_stars; j++) {
                if (!catalog_used[j] && votes[i][j] > best_votes) {
                    best_votes = votes[i][j];
                    best_centroid = i;
                    best_catalog = j;
                }
            }
        }
        if (best_centroid < 0) break;
        matches[num_matches].centroid_idx = best_centroid;
        matches[num_matches].catalog_idx = best_catalog;
        centroid_used[best_centroid] = 1;
        catalog_used[best_catalog] = 1;
        num_matches++;
    }

    if (telem) {
        telem->identify_matches = (uint32_t)num_matches;
    }
    return (num_matches >= 3) ? STARFIX_SUCCESS : STARFIX_ERR_IDENTIFY_NO_MATCH;
}

int starfix_propagate_astrometry(starfix_catalog_star_t* catalog, uint32_t num_stars, double year,
                                 double* cached_year) {
    if (catalog == NULL || num_stars == 0) {
        return STARFIX_ERR_NULL_POINTER;
    }
    if (!isfinite(year)) return STARFIX_ERR_INVALID_PARAM;

    if (cached_year && fabs(year - *cached_year) < (1.0 / 365.25)) {
        return 0; /* cached coordinates are still valid (within ~1 day) */
    }

    double mjd = starfix_year_to_mjd(year);
    uint32_t i;
    for (i = 0; i < num_stars; i++) {
        double ra = catalog[i].ra0;
        double dec = catalog[i].dec0;

        starfix_precess(mjd, &ra, &dec);
        starfix_nutation(mjd, &ra, &dec);
        starfix_aberration(mjd, &ra, &dec);

        starfix_catalog_set_ra(&catalog[i], ra);
        starfix_catalog_set_dec(&catalog[i], dec);
    }

    if (cached_year) {
        *cached_year = year;
    }
    return 0;
}
