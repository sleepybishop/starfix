#ifndef STARFIX_ASTROMETRY_H
#define STARFIX_ASTROMETRY_H

#ifdef __cplusplus
extern "C" {
#endif

/* Converts a fractional year to Modified Julian Date (MJD) */
double starfix_year_to_mjd(double year);

/* precess, nutate, and aberrate a J2000 unit vector for MJD */
void starfix_astrometry_frame(double mjd, double P[3][3], double aberration[3]);

/* apparent to true altitude in radians; pr is mbar, tr is celsius */
void starfix_unrefract(double pr, double tr, double aa, double* ta);

#ifdef __cplusplus
}
#endif
#endif /* STARFIX_ASTROMETRY_H */
