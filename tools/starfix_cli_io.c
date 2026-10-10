#include "starfix_cli_io.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void starfix_rotation_to_euler(const double R[3][3], double* ra_deg, double* dec_deg,
                               double* roll_deg) {
    double z_1 = R[2][0];
    double z_2 = R[2][1];
    double z_3 = R[2][2];

    double dec_rad = asin(z_3 > 1.0 ? 1.0 : (z_3 < -1.0 ? -1.0 : z_3));
    double ra_rad = atan2(z_2, z_1);
    if (ra_rad < 0) ra_rad += 2.0 * M_PI;

    double x_zero[3] = {-sin(ra_rad), cos(ra_rad), 0.0};
    double y_zero[3] = {-sin(dec_rad) * cos(ra_rad), -sin(dec_rad) * sin(ra_rad), cos(dec_rad)};

    double x_cam[3] = {R[0][0], R[0][1], R[0][2]};

    double dot_y = x_cam[0] * y_zero[0] + x_cam[1] * y_zero[1] + x_cam[2] * y_zero[2];
    double dot_x = x_cam[0] * x_zero[0] + x_cam[1] * x_zero[1] + x_cam[2] * x_zero[2];

    double roll_rad = atan2(dot_y, dot_x);
    if (roll_rad < 0) roll_rad += 2.0 * M_PI;

    *ra_deg = ra_rad * 180.0 / M_PI;
    *dec_deg = dec_rad * 180.0 / M_PI;
    *roll_deg = roll_rad * 180.0 / M_PI;
}

unsigned char* read_pgm(const char* path, int* width, int* height) {
    FILE* f = fopen(path, "rb");
    if (f == NULL) return NULL;

    char format[16];
    if (fscanf(f, "%15s", format) != 1 || strcmp(format, "P5") != 0) {
        fclose(f);
        return NULL;
    }

    int w, h, max_val;
    if (fscanf(f, "%d %d %d", &w, &h, &max_val) != 3 || max_val != 255) {
        fclose(f);
        return NULL;
    }

    /* skip single whitespace byte after header */
    fgetc(f);

    unsigned char* data = (unsigned char*)malloc((size_t)w * (size_t)h);
    if (data == NULL) {
        fclose(f);
        return NULL;
    }

    if (fread(data, 1, (size_t)w * (size_t)h, f) != (size_t)w * (size_t)h) {
        free(data);
        fclose(f);
        return NULL;
    }

    fclose(f);
    *width = w;
    *height = h;
    return data;
}

int write_centroids_json(const char* path, const starfix_centroid_t* centroids, int count) {
    FILE* f = fopen(path, "w");
    if (f == NULL) return -1;

    fprintf(f, "[\n");
    int i;
    for (i = 0; i < count; i++) {
        fprintf(f, "  {\n");
        fprintf(f, "    \"u\": %.6f,\n", centroids[i].u);
        fprintf(f, "    \"v\": %.6f,\n", centroids[i].v);
        fprintf(f, "    \"flux\": %.6f\n", centroids[i].flux);
        fprintf(f, "  }%s\n", (i == count - 1) ? "" : ",");
    }
    fprintf(f, "]\n");
    fclose(f);
    return 0;
}

static int read_text(const char* path, char** text) {
    FILE* f = fopen(path, "rb");
    if (!f) return -1;
    int status = -3;
    char* buffer = NULL;
    if (fseek(f, 0, SEEK_END) != 0) goto done;
    long size = ftell(f);
    if (size < 0 || (uintmax_t)size >= SIZE_MAX || fseek(f, 0, SEEK_SET) != 0) goto done;
    buffer = malloc((size_t)size + 1);
    if (!buffer) {
        status = -2;
        goto done;
    }
    if (fread(buffer, 1, (size_t)size, f) != (size_t)size) goto done;
    buffer[size] = '\0';
    *text = buffer;
    status = 0;
done:
    if (status < 0) free(buffer);
    fclose(f);
    return status;
}

/* CLI records are flat numeric objects; terminate each one before looking up its fields */
static char* next_object(char** cursor) {
    char* start = strchr(*cursor, '{');
    if (!start) return NULL;
    char* end = strchr(start, '}');
    if (!end) return NULL;
    *end = '\0';
    *cursor = end + 1;
    return start;
}

static const char* json_value(const char* object, const char* key, const char* fallback) {
    const char* value = strstr(object, key);
    if (value) value = strchr(value, ':');
    return value ? value + 1 : fallback;
}

int read_centroids_json(const char* path, starfix_centroid_t* centroids, int max_count) {
    char* text;
    int status = read_text(path, &text);
    if (status < 0) return status;
    int count = 0;
    char* cursor = text;
    char* object;
    while (count < max_count && (object = next_object(&cursor))) {
        double u = atof(json_value(object, "\"u\"", "-1"));
        double v = atof(json_value(object, "\"v\"", "-1"));
        if (!(u >= 0.0 && v >= 0.0)) continue;
        centroids[count].u = u;
        centroids[count].v = v;
        centroids[count].flux = atof(json_value(object, "\"flux\"", "0"));
        count++;
    }
    free(text);
    return count;
}

int write_identified_json(const char* path, const starfix_match_t* matches, int count,
                          const starfix_catalog_star_t* catalog) {
    FILE* f = fopen(path, "w");
    if (f == NULL) return -1;

    fprintf(f, "[\n");
    int i;
    for (i = 0; i < count; i++) {
        int cat_idx = matches[i].catalog_idx;
        fprintf(f, "  {\n");
        fprintf(f, "    \"centroid_idx\": %d,\n", matches[i].centroid_idx);
        fprintf(f, "    \"catalog_idx\": %d,\n", cat_idx);
        fprintf(f, "    \"hip\": %d,\n", catalog[cat_idx].hip);
        fprintf(f, "    \"ra\": %.8f,\n", catalog[cat_idx].ra);
        fprintf(f, "    \"dec\": %.8f,\n", catalog[cat_idx].dec);
        fprintf(f, "    \"mag\": %.2f\n", catalog[cat_idx].mag);
        fprintf(f, "  }%s\n", (i == count - 1) ? "" : ",");
    }
    fprintf(f, "]\n");
    fclose(f);
    return 0;
}

int read_identified_json(const char* path, starfix_match_t* matches, int max_count) {
    char* text;
    int status = read_text(path, &text);
    if (status < 0) return status;
    int count = 0;
    char* cursor = text;
    char* object;
    while (count < max_count && (object = next_object(&cursor))) {
        int centroid = atoi(json_value(object, "\"centroid_idx\"", "-1"));
        int catalog = atoi(json_value(object, "\"catalog_idx\"", "-1"));
        if (centroid < 0 || catalog < 0) continue;
        matches[count++] = (starfix_match_t){centroid, catalog};
    }
    free(text);
    return count;
}
