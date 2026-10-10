#include "starfix_nmea.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

static void append_checksum(char* sentence, size_t max_len) {
    uint8_t checksum = 0;
    int i = 1; /* skip the starting '$' */

    while (sentence[i] != '\0' && sentence[i] != '*') {
        checksum ^= (uint8_t)sentence[i];
        i++;
    }

    if (sentence[i] == '*') {
        snprintf(&sentence[i + 1], max_len - (size_t)i - 1, "%02X\r\n", checksum);
    } else {
        snprintf(&sentence[i], max_len - (size_t)i, "*%02X\r\n", checksum);
    }
}

/* converts degrees format (DD.DDDDD) to NMEA format (DDMM.MMMM) */
static void deg_to_nmea(double deg, int is_lat, char* out_str, size_t out_size, char* dir) {
    /* round in ten-thousandths of a minute before splitting degrees and minutes */
    unsigned ticks = (unsigned)lround(fabs(deg) * 600000.0);
    *dir = is_lat ? (deg >= 0.0 ? 'N' : 'S') : (deg >= 0.0 ? 'E' : 'W');
    snprintf(out_str, out_size, "%0*u%02u.%04u", is_lat ? 2 : 3, ticks / 600000, ticks / 10000 % 60,
             ticks % 10000);
}

starfix_status_t starfix_format_nmea(double lat_deg, double lon_deg, time_t time_utc,
                                     char* out_buffer, size_t out_size) {
    if (!out_buffer) return STARFIX_ERR_NULL_POINTER;
    if (!isfinite(lat_deg) || !isfinite(lon_deg) || lat_deg < -90 || lat_deg > 90 ||
        lon_deg < -180 || lon_deg > 180)
        return STARFIX_ERR_INVALID_PARAM;

    struct tm* tm_info = gmtime(&time_utc);
    if (tm_info == NULL) {
        return STARFIX_ERR_INVALID_PARAM;
    }

    char lat_str[16];
    char lon_str[16];
    char lat_dir, lon_dir;

    deg_to_nmea(lat_deg, 1, lat_str, sizeof(lat_str), &lat_dir);
    deg_to_nmea(lon_deg, 0, lon_str, sizeof(lon_str), &lon_dir);

    char time_str[32];
    snprintf(time_str, sizeof(time_str), "%02d%02d%02d", tm_info->tm_hour, tm_info->tm_min,
             tm_info->tm_sec);

    char date_str[32];
    snprintf(date_str, sizeof(date_str), "%02d%02d%02d", tm_info->tm_mday, tm_info->tm_mon + 1,
             (tm_info->tm_year + 1900) % 100);

    char gga[128];
    /* 1 = GPS fix, 08 = satellites (mock), 1.0 = HDOP, 0.0 = altitude */
    snprintf(gga, sizeof(gga), "$GPGGA,%s.00,%s,%c,%s,%c,1,08,1.0,0.0,M,0.0,M,,*", time_str,
             lat_str, lat_dir, lon_str, lon_dir);
    append_checksum(gga, sizeof(gga));

    char rmc[128];
    /* A = active, 0.0 = speed, 0.0 = course */
    snprintf(rmc, sizeof(rmc), "$GPRMC,%s.00,A,%s,%c,%s,%c,0.0,0.0,%s,0.0,E*", time_str, lat_str,
             lat_dir, lon_str, lon_dir, date_str);
    append_checksum(rmc, sizeof(rmc));

    size_t gga_len = strlen(gga), rmc_len = strlen(rmc);
    if (out_size <= gga_len + rmc_len) return STARFIX_ERR_BUFFER_FULL;
    memcpy(out_buffer, gga, gga_len);
    memcpy(out_buffer + gga_len, rmc, rmc_len + 1);
    return STARFIX_SUCCESS;
}
