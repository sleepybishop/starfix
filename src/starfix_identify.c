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
#include "starfix_attitude.h"

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
    if (!f) return STARFIX_ERR_MEMORY;
    starfix_catalog_star_t* catalog = NULL;
    starfix_hash_entry_t* entries = NULL;
    int status = STARFIX_ERR_DATA_CORRUPTION;
    if (fseek(f, 0, SEEK_END) != 0) goto done;
    long file_size = ftell(f);
    if (file_size < 16 || fseek(f, 0, SEEK_SET) != 0) goto done;

    unsigned char header[16];
    if (fread(header, 1, sizeof(header), f) != sizeof(header) || memcmp(header, "SFIX", 4) != 0)
        goto done;
    uint32_t n_stars = decode_le32(header + 4);
    uint32_t n_entries = decode_le32(header + 8);
    uint32_t b_factor = decode_le32(header + 12);
    uint64_t expected_size = 16ULL + (uint64_t)n_stars * 12ULL + (uint64_t)n_entries * 12ULL;
    if (n_stars == 0 || n_stars > 65536U || n_entries == 0 || n_entries > (uint32_t)INT_MAX ||
        b_factor == 0 || b_factor > 63U || expected_size != (uint64_t)file_size)
        goto done;

    catalog = malloc((size_t)n_stars * sizeof(*catalog));
    entries = malloc((size_t)n_entries * sizeof(*entries));
    if (!catalog || !entries) {
        status = STARFIX_ERR_MEMORY;
        goto done;
    }

    uint32_t i;
    for (i = 0; i < n_stars; i++) {
        unsigned char record[12];
        if (fread(record, 1, sizeof(record), f) != sizeof(record)) goto done;
        int32_t hip;
        uint32_t hip_bits = decode_le32(record);
        uint16_t ra_q = decode_le16(record + 4);
        uint16_t dec_q = decode_le16(record + 6);
        uint8_t mag_q = record[8];
        memcpy(&hip, &hip_bits, sizeof(hip));

        catalog[i].hip = hip;
        catalog[i].ra0 = ra_q * (2.0 * M_PI) / 65535.0;
        catalog[i].dec0 = dec_q * M_PI / 65535.0 - M_PI / 2.0;
        catalog[i].ra = catalog[i].ra0;
        catalog[i].dec = catalog[i].dec0;
        catalog[i].mag = mag_q * 10.0 / 255.0 - 2.0;
    }

    uint32_t previous_key = 0;
    for (i = 0; i < n_entries; i++) {
        unsigned char record[12];
        if (fread(record, 1, sizeof(record), f) != sizeof(record)) goto done;
        starfix_hash_entry_t* entry = &entries[i];
        entry->hash_key = decode_le32(record);
        int j;
        for (j = 0; j < 4; j++) {
            entry->stars[j] = decode_le16(record + 4 + 2 * j);
            if ((uint32_t)entry->stars[j] >= n_stars) goto done;
        }
        if (i > 0 && entry->hash_key < previous_key) goto done;
        previous_key = entry->hash_key;
    }

    *catalog_stars = catalog;
    *hash_entries = entries;
    *num_stars = n_stars;
    *num_entries = n_entries;
    *bin_factor = b_factor;
    status = STARFIX_SUCCESS;
done:
    if (status < 0) {
        free(catalog);
        free(entries);
    }
    fclose(f);
    return status;
}

static int find_first_hash_entry(uint32_t key, const starfix_hash_entry_t* entries,
                                 uint32_t num_entries) {
    int low = 0, high = (int)num_entries;
    while (low < high) {
        int mid = low + (high - low) / 2;
        if (entries[mid].hash_key < key)
            low = mid + 1;
        else
            high = mid;
    }
    return low < (int)num_entries && entries[low].hash_key == key ? low : -1;
}

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

/* verify a proper rotation and unique extra matches for at most 20 observations */
static int verify_quad(const starfix_vector3_t w_seed[4], const starfix_vector3_t v_seed[4],
                       int num_observed, const starfix_vector3_t* observed, uint32_t num_stars,
                       const starfix_vector3_t* catalog_vectors, double verify_cos,
                       double field_cos, starfix_match_t verified[20], double* error) {
    int i, j, seed;
    double R[3][3];
    starfix_quaternion_t q;
    if (starfix_solve_attitude(4, w_seed, v_seed, &q, R) != STARFIX_SUCCESS) return 0;
    /* reject mirrored seeds before searching for additional matches */
    for (seed = 0; seed < 4; seed++) {
        starfix_vector3_t v = v_seed[seed], w = w_seed[seed];
        double dot = w.x * (R[0][0] * v.x + R[0][1] * v.y + R[0][2] * v.z) +
                     w.y * (R[1][0] * v.x + R[1][1] * v.y + R[1][2] * v.z) +
                     w.z * (R[2][0] * v.x + R[2][1] * v.y + R[2][2] * v.z);
        if (dot < verify_cos) return 0;
    }
    int nearest[20];
    double nearest_dot[20];
    for (i = 0; i < num_observed; i++) {
        nearest[i] = -1;
        nearest_dot[i] = verify_cos;
    }
    for (j = 0; j < (int)num_stars; j++) {
        starfix_vector3_t v = catalog_vectors[j];
        double z = R[2][0] * v.x + R[2][1] * v.y + R[2][2] * v.z;
        if (z < field_cos) continue;
        double x = R[0][0] * v.x + R[0][1] * v.y + R[0][2] * v.z;
        double y = R[1][0] * v.x + R[1][1] * v.y + R[1][2] * v.z;
        for (i = 0; i < num_observed; i++) {
            double dot = x * observed[i].x + y * observed[i].y + z * observed[i].z;
            if (dot > nearest_dot[i]) {
                nearest_dot[i] = dot;
                nearest[i] = j;
            }
        }
    }
    int count = 0;
    *error = 0.0;
    for (i = 0; i < num_observed; i++) {
        if (nearest[i] < 0) continue;
        int duplicate = 0;
        for (j = 0; j < num_observed; j++) {
            if (j != i && nearest[j] == nearest[i] &&
                (nearest_dot[j] > nearest_dot[i] || (nearest_dot[j] >= nearest_dot[i] && j < i))) {
                duplicate = 1;
                break;
            }
        }
        if (duplicate) continue;
        verified[count++] = (starfix_match_t){i, nearest[i]};
        *error += 1.0 - nearest_dot[i];
    }
    return count;
}

typedef struct {
    int count;
    double error;
    starfix_match_t matches[20];
} match_set_t;

static const int pairs[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};

static starfix_vector3_t catalog_vector(const starfix_catalog_star_t* star) {
    return (starfix_vector3_t){cos(star->dec) * cos(star->ra), cos(star->dec) * sin(star->ra),
                               sin(star->dec)};
}

static double separation(starfix_vector3_t a, starfix_vector3_t b) {
    double dot = a.x * b.x + a.y * b.y + a.z * b.z;
    if (dot > 1.0) dot = 1.0;
    if (dot < -1.0) dot = -1.0;
    return acos(dot);
}

static int track_stars(int count, const starfix_vector3_t* camera, uint32_t num_stars,
                       const starfix_catalog_star_t* catalog, const double R[3][3], double fov,
                       int capacity, starfix_match_t* matches) {
    int used[20] = {0}, found = 0;
    double tolerance = cos(0.5 * M_PI / 180.0);
    double field = cos(fov / 2.0 + M_PI / 180.0);
    for (int j = 0; j < (int)num_stars && found < capacity; j++) {
        starfix_vector3_t v = catalog_vector(&catalog[j]);
        double x = R[0][0] * v.x + R[0][1] * v.y + R[0][2] * v.z;
        double y = R[1][0] * v.x + R[1][1] * v.y + R[1][2] * v.z;
        double z = R[2][0] * v.x + R[2][1] * v.y + R[2][2] * v.z;
        if (z < field) continue;
        int nearest = -1;
        double best_dot = -1.0;
        for (int i = 0; i < count; i++) {
            double dot = x * camera[i].x + y * camera[i].y + z * camera[i].z;
            if (dot > best_dot) {
                best_dot = dot;
                nearest = i;
            }
        }
        if (nearest < 0 || best_dot < tolerance || used[nearest]) continue;
        matches[found++] = (starfix_match_t){nearest, j};
        used[nearest] = 1;
    }
    return found;
}

static int next_quad(int indices[4], int count) {
    for (int i = 3; i >= 0; i--) {
        if (indices[i] >= count - 4 + i) continue;
        indices[i]++;
        for (int j = i + 1; j < 4; j++) indices[j] = indices[j - 1] + 1;
        return 1;
    }
    return 0;
}

/* exact hash first, then distinct neighbors one bin away, in search order */
static int quad_keys(const double angles[6], uint32_t bins, uint32_t keys[11]) {
    int code[5];
    for (int i = 0; i < 5; i++) code[i] = (int)((angles[i] / angles[5]) * bins);
    keys[0] = compress_hash_key(code);
    int count = 1;
    for (int step = -1; step <= 1; step += 2) {
        for (int dim = 0; dim < 5; dim++) {
            int adjacent[5];
            memcpy(adjacent, code, sizeof(adjacent));
            adjacent[dim] += step;
            if (adjacent[dim] < 0) adjacent[dim] = 0;
            if (adjacent[dim] > (int)bins) adjacent[dim] = (int)bins;
            uint32_t key = compress_hash_key(adjacent);
            int i = 0;
            while (i < count && keys[i] != key) i++;
            if (i == count) keys[count++] = key;
        }
    }
    return count;
}

static void match_quad(const starfix_vector3_t w[4], const double angles[6], const int order[6],
                       const starfix_vector3_t v[4], int count, const starfix_vector3_t* camera,
                       uint32_t num_stars, const starfix_vector3_t* catalog, double field_cos,
                       match_set_t* best) {
    static const int perms[24][4] = {
        {0, 1, 2, 3}, {0, 1, 3, 2}, {0, 2, 1, 3}, {0, 2, 3, 1}, {0, 3, 1, 2}, {0, 3, 2, 1},
        {1, 0, 2, 3}, {1, 0, 3, 2}, {1, 2, 0, 3}, {1, 2, 3, 0}, {1, 3, 0, 2}, {1, 3, 2, 0},
        {2, 0, 1, 3}, {2, 0, 3, 1}, {2, 1, 0, 3}, {2, 1, 3, 0}, {2, 3, 0, 1}, {2, 3, 1, 0},
        {3, 0, 1, 2}, {3, 0, 2, 1}, {3, 1, 0, 2}, {3, 1, 2, 0}, {3, 2, 0, 1}, {3, 2, 1, 0}};

    double catalog_angles[4][4] = {{0}};
    for (int i = 0; i < 6; i++) {
        int a = pairs[i][0], b = pairs[i][1];
        catalog_angles[a][b] = catalog_angles[b][a] = separation(v[a], v[b]);
    }
    double tolerance = 0.5 * M_PI / 180.0;
    double verify_cos = cos(0.12 * M_PI / 180.0);
    for (int p = 0; p < 24; p++) {
        int i = 0;
        for (; i < 6; i++) {
            int a = perms[p][pairs[order[i]][0]], b = perms[p][pairs[order[i]][1]];
            if (fabs(angles[i] - catalog_angles[a][b]) > tolerance) break;
        }
        if (i < 6) continue;
        starfix_vector3_t seed[4];
        for (i = 0; i < 4; i++) seed[i] = v[perms[p][i]];
        match_set_t fit;
        fit.count = verify_quad(w, seed, count, camera, num_stars, catalog, verify_cos, field_cos,
                                fit.matches, &fit.error);
        /* the fifth star independently validates the seed */
        if (fit.count >= 5 &&
            (fit.count > best->count || (fit.count == best->count && fit.error < best->error))) {
            best->count = fit.count;
            best->error = fit.error;
            memcpy(best->matches, fit.matches, (size_t)fit.count * sizeof(*fit.matches));
        }
    }
}

starfix_status_t starfix_identify_stars(
    int num_centroids, const starfix_centroid_t* centroids, double width, double height,
    double fov_deg, uint32_t num_stars, const starfix_catalog_star_t* catalog_stars,
    uint32_t num_entries, const starfix_hash_entry_t* hash_entries, uint32_t bin_factor,
    int max_matches, starfix_match_t* matches, starfix_arena_t* arena, starfix_telemetry_t* telem,
    const double (*attitude_hint)[3]) {
    if (num_centroids < 3 || centroids == NULL || catalog_stars == NULL || matches == NULL) {
        return STARFIX_ERR_NULL_POINTER;
    }
    if (width <= 0.0 || height <= 0.0 || fov_deg <= 0.0 || fov_deg >= 180.0 || !isfinite(width) ||
        !isfinite(height) || !isfinite(fov_deg) || num_stars == 0 ||
        num_stars > (uint32_t)INT_MAX || max_matches <= 0) {
        return STARFIX_ERR_INVALID_PARAM;
    }

    int observed = num_centroids < 20 ? num_centroids : 20;
    double fov = fov_deg * M_PI / 180.0;
    double focal = width / (2.0 * tan(fov / 2.0));
    starfix_vector3_t camera[20];
    for (int i = 0; i < observed; i++) {
        double x = (centroids[i].u - width / 2.0) / focal;
        double y = (centroids[i].v - height / 2.0) / focal;
        double norm = sqrt(x * x + y * y + 1.0);
        camera[i] = (starfix_vector3_t){x / norm, y / norm, 1.0 / norm};
    }
    if (attitude_hint) {
        int count = track_stars(observed, camera, num_stars, catalog_stars, attitude_hint, fov,
                                max_matches, matches);
        if (count >= 3) {
            if (telem) telem->identify_matches = (uint32_t)count;
            return STARFIX_SUCCESS;
        }
    }
    if (num_centroids < 4) return STARFIX_ERR_NULL_POINTER;
    if (!hash_entries || num_entries == 0 || bin_factor == 0 || bin_factor > 63U || !arena)
        return STARFIX_ERR_INVALID_PARAM;

    starfix_vector3_t* catalog = starfix_alloc_array(arena, starfix_vector3_t, num_stars);
    if (!catalog) return STARFIX_ERR_MEMORY;
    for (uint32_t i = 0; i < num_stars; i++) catalog[i] = catalog_vector(&catalog_stars[i]);
    double corner = atan(hypot(width, height) / (2.0 * focal));
    double field_cos = cos(corner + 0.5 * M_PI / 180.0);
    match_set_t best = {.count = 0, .error = HUGE_VAL};

    /* seed from the brightest twelve detections; verify against up to twenty */
    int seeds = num_centroids < 12 ? num_centroids : 12;
    int indices[4] = {0, 1, 2, 3};
    do {
        starfix_vector3_t w[4];
        for (int i = 0; i < 4; i++) w[i] = camera[indices[i]];
        double angles[6];
        int order[6] = {0, 1, 2, 3, 4, 5};
        for (int i = 0; i < 6; i++) angles[i] = separation(w[pairs[i][0]], w[pairs[i][1]]);
        sort_angles(angles, order, 6);
        if (angles[5] < 1e-5) continue;

        uint32_t keys[11];
        int key_count = quad_keys(angles, bin_factor, keys);
        for (int k = 0; k < key_count; k++) {
            int first = find_first_hash_entry(keys[k], hash_entries, num_entries);
            if (first < 0) continue;
            for (int at = first;
                 at < (int)num_entries && at - first < 100 && hash_entries[at].hash_key == keys[k];
                 at++) {
                const starfix_hash_entry_t* entry = &hash_entries[at];
                starfix_vector3_t v[4];
                int i = 0;
                for (; i < 4; i++) {
                    if (entry->stars[i] >= num_stars) break;
                    v[i] = catalog[entry->stars[i]];
                }
                if (i < 4) continue;
                match_quad(w, angles, order, v, observed, camera, num_stars, catalog, field_cos,
                           &best);
            }
        }
    } while (next_quad(indices, seeds));

    int count = best.count < max_matches ? best.count : max_matches;
    if (count > 0) memcpy(matches, best.matches, (size_t)count * sizeof(*matches));
    if (telem) telem->identify_matches = (uint32_t)count;
    return count >= 3 ? STARFIX_SUCCESS : STARFIX_ERR_IDENTIFY_NO_MATCH;
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

        catalog[i].ra = ra;
        catalog[i].dec = dec;
    }

    if (cached_year) {
        *cached_year = year;
    }
    return 0;
}
