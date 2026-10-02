#include "starfix_cli_io.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void starfix_rotation_to_euler(const double R[9], double* ra_deg, double* dec_deg,
                               double* roll_deg) {
    double z_1 = R[6];
    double z_2 = R[7];
    double z_3 = R[8];

    double dec_rad = asin(z_3 > 1.0 ? 1.0 : (z_3 < -1.0 ? -1.0 : z_3));
    double ra_rad = atan2(z_2, z_1);
    if (ra_rad < 0) ra_rad += 2.0 * M_PI;

    double x_zero[3] = {-sin(ra_rad), cos(ra_rad), 0.0};
    double y_zero[3] = {-sin(dec_rad) * cos(ra_rad), -sin(dec_rad) * sin(ra_rad), cos(dec_rad)};

    double x_cam[3] = {R[0], R[1], R[2]};

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

int read_centroids_json(const char* path, starfix_centroid_t* centroids, int max_count) {
    FILE* f = fopen(path, "r");
    if (f == NULL) return -1;

    fseek(f, 0, SEEK_END);
    long sz = ftell(f);
    fseek(f, 0, SEEK_SET);
    char* buf = (char*)malloc((size_t)sz + 1);
    if (buf == NULL) {
        fclose(f);
        return -2;
    }
    if (fread(buf, 1, (size_t)sz, f) != (size_t)sz) {
        free(buf);
        fclose(f);
        return -3;
    }
    buf[sz] = '\0';
    fclose(f);

    int count = 0;
    char* obj_start = buf;
    while (count < max_count) {
        obj_start = strchr(obj_start, '{');
        if (!obj_start) break;
        char* obj_end = strchr(obj_start, '}');
        if (!obj_end) break;

        double u = -1.0, v = -1.0, flux = 0.0;
        char* u_ptr = strstr(obj_start, "\"u\"");
        if (u_ptr && u_ptr < obj_end) {
            char* p = strchr(u_ptr, ':');
            if (p) u = atof(p + 1);
        }
        char* v_ptr = strstr(obj_start, "\"v\"");
        if (v_ptr && v_ptr < obj_end) {
            char* p = strchr(v_ptr, ':');
            if (p) v = atof(p + 1);
        }
        char* flux_ptr = strstr(obj_start, "\"flux\"");
        if (flux_ptr && flux_ptr < obj_end) {
            char* p = strchr(flux_ptr, ':');
            if (p) flux = atof(p + 1);
        }

        if (u >= 0.0 && v >= 0.0) {
            centroids[count].u = u;
            centroids[count].v = v;
            centroids[count].flux = flux;
            count++;
        }
        obj_start = obj_end + 1;
    }
    free(buf);
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
        fprintf(f, "    \"ra\": %.8f,\n", starfix_catalog_ra(&catalog[cat_idx]));
        fprintf(f, "    \"dec\": %.8f,\n", starfix_catalog_dec(&catalog[cat_idx]));
        fprintf(f, "    \"mag\": %.2f\n", starfix_catalog_mag(&catalog[cat_idx]));
        fprintf(f, "  }%s\n", (i == count - 1) ? "" : ",");
    }
    fprintf(f, "]\n");
    fclose(f);
    return 0;
}

int read_identified_json(const char* path, starfix_match_t* matches, int max_count) {
    FILE* f = fopen(path, "r");
    if (f == NULL) return -1;

    fseek(f, 0, SEEK_END);
    long sz = ftell(f);
    fseek(f, 0, SEEK_SET);
    char* buf = (char*)malloc((size_t)sz + 1);
    if (buf == NULL) {
        fclose(f);
        return -2;
    }
    if (fread(buf, 1, (size_t)sz, f) != (size_t)sz) {
        free(buf);
        fclose(f);
        return -3;
    }
    buf[sz] = '\0';
    fclose(f);

    int count = 0;
    char* obj_start = buf;
    while (count < max_count) {
        obj_start = strchr(obj_start, '{');
        if (!obj_start) break;
        char* obj_end = strchr(obj_start, '}');
        if (!obj_end) break;

        int centroid_idx = -1, catalog_idx = -1;
        char* c_ptr = strstr(obj_start, "\"centroid_idx\"");
        if (c_ptr && c_ptr < obj_end) {
            char* p = strchr(c_ptr, ':');
            if (p) centroid_idx = atoi(p + 1);
        }
        char* cat_ptr = strstr(obj_start, "\"catalog_idx\"");
        if (cat_ptr && cat_ptr < obj_end) {
            char* p = strchr(cat_ptr, ':');
            if (p) catalog_idx = atoi(p + 1);
        }

        if (centroid_idx >= 0 && catalog_idx >= 0) {
            matches[count].centroid_idx = centroid_idx;
            matches[count].catalog_idx = catalog_idx;
            count++;
        }
        obj_start = obj_end + 1;
    }
    free(buf);
    return count;
}
