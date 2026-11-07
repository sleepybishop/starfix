#ifndef STARFIX_NMEA_H
#define STARFIX_NMEA_H

#include <stddef.h>
#include <time.h>

#include "starfix_status.h"

#ifdef __cplusplus
extern "C" {
#endif

/* format degree lat/lon (north/east positive) and unix UTC time as GPGGA/GPRMC with CRLF.
 * out_size includes the terminator (256 recommended); short buffers return BUFFER_FULL.
 * returns zero or a negative error; output is unchanged on error */
starfix_status_t starfix_format_nmea(double lat_deg, double lon_deg, time_t time_utc,
                                     char* out_buffer, size_t out_size);

/* format GPGGA/GPRMC with explicit fix quality (1 = star fix, 6 = dead reckoning),
 * matched star count (0..99), and optional speed (knots) / course (degrees true);
 * pass NAN for speed_knots or course_deg when unavailable to omit those fields */
starfix_status_t starfix_format_nmea_fix(double lat_deg, double lon_deg, time_t time_utc,
                                         int fix_quality, int num_stars, double speed_knots,
                                         double course_deg, char* out_buffer, size_t out_size);

#ifdef __cplusplus
}
#endif
#endif /* STARFIX_NMEA_H */
