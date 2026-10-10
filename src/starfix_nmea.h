#ifndef STARFIX_NMEA_H
#define STARFIX_NMEA_H

#include <stddef.h>
#include <stdint.h>
#include <time.h>

#include "starfix_status.h"

#ifdef __cplusplus
extern "C" {
#endif

/* format degree lat/lon (north/east positive) and unix UTC time as GPGGA/GPRMC with CRLF.
 * out_size is in bytes (256 recommended); returns zero or a negative error */
int starfix_format_nmea(double lat_deg, double lon_deg, time_t time_utc, char* out_buffer,
                        size_t out_size);

#ifdef __cplusplus
}
#endif
#endif /* STARFIX_NMEA_H */
