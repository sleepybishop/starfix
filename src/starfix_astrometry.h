#ifndef STARFIX_ASTROMETRY_H
#define STARFIX_ASTROMETRY_H

#ifdef __cplusplus
extern "C" {
#endif

/* Converts a fractional year to Modified Julian Date (MJD) */
double starfix_year_to_mjd(double year);

/* precess J2000 ra/dec in radians to MJD, in place */
void starfix_precess(double mjd, double* ra, double* dec);

/* apply IAU 1980 nutation to ra/dec in radians, in place */
void starfix_nutation(double mjd, double* ra, double* dec);

/* apply annual aberration to ra/dec in radians, in place */
void starfix_aberration(double mjd, double* ra, double* dec);

/* true to apparent altitude in radians; pr is mbar, tr is celsius */
void starfix_refract(double pr, double tr, double ta, double* aa);

/* apparent to true altitude in radians; pr is mbar, tr is celsius */
void starfix_unrefract(double pr, double tr, double aa, double* ta);

#ifdef __cplusplus
}
#endif
#endif /* STARFIX_ASTROMETRY_H */
