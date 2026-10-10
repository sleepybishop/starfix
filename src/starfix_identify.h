#ifndef STARFIX_IDENTIFY_H
#define STARFIX_IDENTIFY_H

#include <stdint.h>

#include "starfix_arena.h"
#include "starfix_centroid.h"
#include "starfix_status.h"

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

/* load catalog and sorted hashes; caller frees both arrays; returns zero or negative error */
int starfix_load_db(const char* db_path, starfix_catalog_star_t** catalog_stars,
                    uint32_t* num_stars, starfix_hash_entry_t** hash_entries, uint32_t* num_entries,
                    uint32_t* bin_factor);

/* propagate J2000 coordinates to year; optional cached_year skips changes < one day */
int starfix_propagate_astrometry(starfix_catalog_star_t* catalog, uint32_t num_stars, double year,
                                 double* cached_year);

/* match pixel centroids to catalog indices; hash_entries must be sorted, matches holds max_matches.
 * fov_deg is horizontal; camera axes are +X right, +Y down, +Z forward.
 * blind matching needs a quad plus an independent fifth star; tracking can use three matches.
 * returns zero or a negative error; telem->identify_matches gives the count */
starfix_status_t starfix_identify_stars(
    int num_centroids, const starfix_centroid_t* centroids, double width, double height,
    double fov_deg, uint32_t num_stars, const starfix_catalog_star_t* catalog_stars,
    uint32_t num_entries, const starfix_hash_entry_t* hash_entries, uint32_t bin_factor,
    int max_matches, starfix_match_t* matches, starfix_arena_t* arena, starfix_telemetry_t* telem,
    const double (*attitude_hint)[3]);

#endif /* STARFIX_IDENTIFY_H */
