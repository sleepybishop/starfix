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

typedef struct {
    uint32_t hash_key;
    uint16_t stars[4];
} starfix_hash_entry_t;

typedef struct {
    starfix_catalog_star_t* stars;
    starfix_hash_entry_t* hashes;
    uint32_t num_stars, num_hashes, bins;
} starfix_db_t;

typedef struct {
    int centroid_idx; /* index of camera centroid */
    int catalog_idx;  /* index in db.stars */
} starfix_match_t;

/* load into an empty database; failure leaves it zeroed; returns zero or negative error */
int starfix_load_db(const char* db_path, starfix_db_t* db);

/* free both arrays and zero the database */
void starfix_free_db(starfix_db_t* db);

/* propagate J2000 coordinates to year */
int starfix_propagate_astrometry(starfix_catalog_star_t* catalog, uint32_t num_stars, double year);

/* return up to max_matches matches (>= 3), or a negative error; hashes must be sorted.
 * fov_deg is horizontal; camera axes are +X right, +Y down, +Z forward.
 * blind matching needs a quad plus an independent fifth star; tracking needs three matches */
int starfix_identify_stars(int num_centroids, const starfix_centroid_t* centroids, double width,
                           double height, double fov_deg, const starfix_db_t* db, int max_matches,
                           starfix_match_t* matches, starfix_arena_t arena,
                           const double (*attitude_hint)[3]);

#endif /* STARFIX_IDENTIFY_H */
