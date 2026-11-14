#ifndef STARFIX_IDENTIFY_H
#define STARFIX_IDENTIFY_H

#include <stdint.h>

#include "starfix_centroid.h"
#include "starfix_geometry.h"
#include "starfix_status.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int32_t hip;            /* HIP ID (-1 if unknown) */
    double ra;              /* propagated RA in radians (0 to 2*pi) */
    double dec;             /* propagated dec in radians (-pi/2 to pi/2) */
    double mag;             /* visual magnitude */
    double ra0;             /* base J2000 RA in radians */
    double dec0;            /* base J2000 dec in radians */
    double pm_ra_rad_yr;    /* tangent east proper motion in rad/yr */
    double pm_dec_rad_yr;   /* tangent north proper motion in rad/yr */
    starfix_vector3_t unit; /* propagated celestial unit vector */
} starfix_catalog_star_t;

static inline starfix_catalog_star_t starfix_make_catalog_star(int32_t hip, double ra, double dec,
                                                               double mag) {
    starfix_catalog_star_t star = {
        .hip = hip,
        .ra = ra,
        .dec = dec,
        .mag = mag,
        .ra0 = ra,
        .dec0 = dec,
        .pm_ra_rad_yr = 0.0,
        .pm_dec_rad_yr = 0.0,
        .unit = starfix_celestial_ray(ra, dec),
    };
    return star;
}

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
starfix_status_t starfix_load_db(const char* db_path, starfix_db_t* db);

/* free both arrays and zero the database */
void starfix_free_db(starfix_db_t* db);

/* propagate J2000 coordinates and proper motions to year */
starfix_status_t starfix_propagate_astrometry(starfix_catalog_star_t* catalog, uint32_t num_stars,
                                              double year);

/* return up to max_matches (>= 3) matches, or a negative error; blind hashes must be sorted.
 * fov_deg is horizontal; camera axes are +X right, +Y down, +Z forward.
 * blind matching needs a quad plus an independent fifth star; tracking needs three matches */
int starfix_identify_stars(int num_centroids, const starfix_centroid_t* centroids, double width,
                           double height, double fov_deg, const starfix_db_t* db, int max_matches,
                           starfix_match_t* matches, const double (*attitude_hint)[3]);

#ifdef __cplusplus
}
#endif
#endif /* STARFIX_IDENTIFY_H */
