#ifndef STARFIX_NMEA_H
#define STARFIX_NMEA_H

#include <stddef.h>
#include <stdint.h>
#include <time.h>

#include "starfix_status.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Formats latitude, longitude, and time into NMEA 0183 standard sentences.
   The resulting string will contain both $GPGGA and $GPRMC sentences separated by \r\n.
   - lat_deg: solved latitude in degrees (positive for North)
   - lon_deg: solved longitude in degrees (positive for East)
   - time_utc: Unix timestamp of the observation
   - out_buffer: allocated character buffer to receive the NMEA string
   - out_size: size of the output buffer (recommend at least 256 bytes)
   - returns: 0 on success, negative on error. */
int starfix_format_nmea(double lat_deg, double lon_deg, time_t time_utc, char* out_buffer,
                        size_t out_size);

#ifdef __cplusplus
}
#endif
#endif /* STARFIX_NMEA_H */
