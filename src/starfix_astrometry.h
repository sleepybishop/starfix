#ifndef STARFIX_ASTROMETRY_H
#define STARFIX_ASTROMETRY_H

#ifdef __cplusplus
extern "C" {
#endif

/* convert fractional year to MJD */
double starfix_year_to_mjd(double year);

/* convert calendar UTC date-time to MJD */
double starfix_utc_to_mjd(int year, int month, int day, int hour, int minute, double second);

/* compute greenwich apparent sidereal time (GAST) in [0, 360) degrees from MJD UTC */
double starfix_gast_deg(double mjd_utc);

/* precess, nutate, and aberrate a J2000 unit vector for MJD */
void starfix_astrometry_frame(double mjd, double P[3][3], double aberration[3]);

/* apparent to true altitude in radians; pr is mbar, tr is celsius */
void starfix_unrefract(double pr, double tr, double aa, double* ta);

/* true to apparent altitude in radians; pr is mbar, tr is celsius */
void starfix_refract(double pr, double tr, double ta, double* aa);

#ifdef __cplusplus
}
#endif
#endif /* STARFIX_ASTROMETRY_H */
