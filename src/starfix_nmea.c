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

starfix_status_t starfix_format_nmea_fix(double lat_deg, double lon_deg, time_t time_utc,
                                         int fix_quality, int num_stars, double speed_knots,
                                         double course_deg, char* out_buffer, size_t out_size) {
    if (!out_buffer) return STARFIX_ERR_NULL_POINTER;
    if (!isfinite(lat_deg) || !isfinite(lon_deg) || lat_deg < -90 || lat_deg > 90 ||
        lon_deg < -180 || lon_deg > 180)
        return STARFIX_ERR_INVALID_PARAM;
    if (fix_quality < 0 || fix_quality > 8 || num_stars < 0 || num_stars > 99)
        return STARFIX_ERR_INVALID_PARAM;
    if ((!isnan(speed_knots) && (!isfinite(speed_knots) || speed_knots < 0.0)) ||
        (!isnan(course_deg) && !isfinite(course_deg)))
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

    char speed_str[16] = "";
    if (!isnan(speed_knots)) {
        snprintf(speed_str, sizeof(speed_str), "%.1f", speed_knots);
    }

    char course_str[16] = "";
    if (!isnan(course_deg)) {
        double norm_course = fmod(course_deg, 360.0);
        if (norm_course < 0.0) norm_course += 360.0;
        snprintf(course_str, sizeof(course_str), "%.1f", norm_course);
    }

    char gga[128];
    /* fix_quality: 1 = celestial star fix, 6 = dead reckoning; num_stars = matched inliers */
    snprintf(gga, sizeof(gga), "$GPGGA,%s.00,%s,%c,%s,%c,%d,%02d,1.0,0.0,M,0.0,M,,*", time_str,
             lat_str, lat_dir, lon_str, lon_dir, fix_quality, num_stars);
    append_checksum(gga, sizeof(gga));

    char rmc[128];
    /* A = active, speed_str = knots (when available), course_str = true heading (when available) */
    snprintf(rmc, sizeof(rmc), "$GPRMC,%s.00,A,%s,%c,%s,%c,%s,%s,%s,0.0,E*", time_str, lat_str,
             lat_dir, lon_str, lon_dir, speed_str, course_str, date_str);
    append_checksum(rmc, sizeof(rmc));

    size_t gga_len = strlen(gga), rmc_len = strlen(rmc);
    if (out_size <= gga_len + rmc_len) return STARFIX_ERR_BUFFER_FULL;
    memcpy(out_buffer, gga, gga_len);
    memcpy(out_buffer + gga_len, rmc, rmc_len + 1);
    return STARFIX_SUCCESS;
}

starfix_status_t starfix_format_nmea(double lat_deg, double lon_deg, time_t time_utc,
                                     char* out_buffer, size_t out_size) {
    return starfix_format_nmea_fix(lat_deg, lon_deg, time_utc, 1, 8, 0.0, 0.0, out_buffer,
                                   out_size);
}
