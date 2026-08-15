#ifndef STARFIX_ASTROMETRY_H
#define STARFIX_ASTROMETRY_H

#ifdef __cplusplus
extern "C" {
#endif

/* Converts a fractional year to Modified Julian Date (MJD) */
double starfix_year_to_mjd(double year);

/* Precesses coordinates from J2000 to the given MJD
 * ra, dec: in radians (modified in place) */
void starfix_precess(double mjd, double* ra, double* dec);

/* Applies IAU 1980 nutation (wobble) to coordinates
 * ra, dec: in radians (modified in place) */
void starfix_nutation(double mjd, double* ra, double* dec);

/* Applies annual aberration (starlight displacement) to coordinates
 * ra, dec: in radians (modified in place) */
void starfix_aberration(double mjd, double* ra, double* dec);

/* Converts true altitude to apparent altitude (refraction)
 * pr: atmospheric pressure in mbar
 * tr: temperature in degrees Celsius
 * ta: true altitude in radians
 * aa: apparent altitude in radians (output) */
void starfix_refract(double pr, double tr, double ta, double* aa);

/* Converts apparent altitude to true altitude (unrefract)
 * pr: atmospheric pressure in mbar
 * tr: temperature in degrees Celsius
 * aa: apparent altitude in radians
 * ta: true altitude in radians (output) */
void starfix_unrefract(double pr, double tr, double aa, double* ta);

#ifdef __cplusplus
}
#endif
#endif /* STARFIX_ASTROMETRY_H */
