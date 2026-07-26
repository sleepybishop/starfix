#ifndef STARFIX_IDENTIFY_H
#define STARFIX_IDENTIFY_H

#include <stdint.h>

#include "starfix_arena.h"
#include "starfix_centroid.h"
#include "starfix_status.h"

/* starfix_identify: C API for lost-in-space star identification using TETRA.
   loads a flat binary database and matches detected sub-pixel camera centroids
   to catalog star indices by hashing 4-star combinations (quads). */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct {
    int32_t hip; /* HIP ID (-1 if unknown) */
    double ra;   /* Precessed RA in radians (0 to 2pi) */
    double dec;  /* Precessed Dec in radians (-pi/2 to pi/2) */
    double mag;  /* Magnitude */
    double ra0;  /* Base J2000 RA in radians */
    double dec0; /* Base J2000 Dec in radians */
} starfix_catalog_star_t;

/* Inline helper functions for quantization/dequantization */
static inline double starfix_catalog_ra(const starfix_catalog_star_t* star) { return star->ra; }

static inline double starfix_catalog_dec(const starfix_catalog_star_t* star) { return star->dec; }

static inline double starfix_catalog_mag(const starfix_catalog_star_t* star) { return star->mag; }

static inline void starfix_catalog_set_ra(starfix_catalog_star_t* star, double ra) {
    star->ra = ra;
}

static inline void starfix_catalog_set_dec(starfix_catalog_star_t* star, double dec) {
    star->dec = dec;
}

static inline void starfix_catalog_set_mag(starfix_catalog_star_t* star, double mag) {
    star->mag = mag;
}

typedef struct {
    uint32_t hash_key;
    uint16_t stars[4];
} starfix_hash_entry_t;

typedef struct {
    int centroid_idx; /* index of camera centroid */
    int catalog_idx;  /* index in catalog_stars array */
} starfix_match_t;

/* loads the binary database:
   - db_path: path to starfix_db.bin
   - catalog_stars: output pointer to allocated star catalog array
   - num_stars: output pointer to number of stars in catalog
   - hash_entries: output pointer to allocated hash entries array
   - num_entries: output pointer to number of hash entries
   - bin_factor: output pointer to bin factor used in hashing
   - returns: 0 on success, negative value on error */
int starfix_load_db(const char* db_path, starfix_catalog_star_t** catalog_stars,
                    uint32_t* num_stars, starfix_hash_entry_t** hash_entries, uint32_t* num_entries,
                    uint32_t* bin_factor);

/* propagates the entire catalog's J2000 coordinates to the specified year using IAU 1976 precession
   matrix. if cached_year is provided, it skips computation if the year hasn't meaningfully changed
   (>1 day). returns 0 on success. */
int starfix_propagate_precession(starfix_catalog_star_t* catalog, uint32_t num_stars, double year,
                                 double* cached_year);

/* performs lost-in-space star identification using TETRA:
   - num_centroids: number of detected centroids
   - centroids: array of detected centroids
   - width: sensor width in pixels
   - height: sensor height in pixels
   - fov_deg: camera field-of-view in degrees
   - num_stars: number of stars in database catalog
   - catalog_stars: star catalog array
   - num_entries: number of database hash entries
   - hash_entries: sorted database hash entries array
   - bin_factor: database bin factor
   - max_matches: maximum number of matches to output
   - matches: output array of matched pairs
   - returns: number of successfully identified stars, or negative value on error */
starfix_status_t starfix_identify_stars(
    int num_centroids, const starfix_centroid_t* centroids, double width, double height,
    double fov_deg, uint32_t num_stars, const starfix_catalog_star_t* catalog_stars,
    uint32_t num_entries, const starfix_hash_entry_t* hash_entries, uint32_t bin_factor,
    int max_matches, starfix_match_t* matches, starfix_arena_t* arena, starfix_telemetry_t* telem,
    const double (*attitude_hint)[3]);

#endif /* STARFIX_IDENTIFY_H */
