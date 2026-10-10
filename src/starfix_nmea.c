#include "starfix_nmea.h"

#include <math.h>
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
    double abs_deg = fabs(deg);
    int d = (int)abs_deg;
    double m = (abs_deg - (double)d) * 60.0;

    if (is_lat) {
        *dir = (deg >= 0.0) ? 'N' : 'S';
        snprintf(out_str, out_size, "%02d%07.4f", d, m);
    } else {
        *dir = (deg >= 0.0) ? 'E' : 'W';
        snprintf(out_str, out_size, "%03d%07.4f", d, m);
    }
}

int starfix_format_nmea(double lat_deg, double lon_deg, time_t time_utc, char* out_buffer,
                        size_t out_size) {
    if (out_buffer == NULL || out_size == 0) {
        return STARFIX_ERR_NULL_POINTER;
    }

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
    /* $GPGGA,hhmmss.ss,llll.ll,a,yyyyy.yy,a,x,xx,x.x,x.x,M,x.x,M,x.x,xxxx*hh */
    /* 1 = GPS fix, 08 = satellites (mock), 1.0 = HDOP, 0.0 = altitude */
    snprintf(gga, sizeof(gga), "$GPGGA,%s.00,%s,%c,%s,%c,1,08,1.0,0.0,M,0.0,M,,*", time_str,
             lat_str, lat_dir, lon_str, lon_dir);
    append_checksum(gga, sizeof(gga));

    char rmc[128];
    /* $GPRMC,hhmmss.ss,A,llll.ll,a,yyyyy.yy,a,x.x,x.x,ddmmyy,x.x,a*hh */
    /* A = active, 0.0 = speed, 0.0 = course */
    snprintf(rmc, sizeof(rmc), "$GPRMC,%s.00,A,%s,%c,%s,%c,0.0,0.0,%s,0.0,E*", time_str, lat_str,
             lat_dir, lon_str, lon_dir, date_str);
    append_checksum(rmc, sizeof(rmc));

    snprintf(out_buffer, out_size, "%s%s", gga, rmc);
    return STARFIX_SUCCESS;
}
