#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "starfix_arena.h"
#include "starfix_attitude.h"
#include "starfix_centroid.h"
#include "starfix_ekf.h"
#include "starfix_graph.h"
#include "starfix_identify.h"
#include "starfix_solver.h"

/* 10 MB statically allocated memory pool for the entire pipeline */
static uint8_t mempool[10 * 1024 * 1024];
static starfix_arena_t arena;

#define RESET_ARENA()                          \
    do {                                       \
        arena.beg = mempool;                   \
        arena.end = mempool + sizeof(mempool); \
    } while (0)

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* helper to extract euler angles from a 3x3 rotation matrix */
static void starfix_rotation_to_euler(const double R[9], double* ra_deg, double* dec_deg,
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

/* helper to read PGM image (P5 binary format) */
static unsigned char* read_pgm(const char* path, int* width, int* height) {
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

/* helper to write detected centroids to JSON */
static int write_centroids_json(const char* path, const starfix_centroid_t* centroids, int count) {
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

/* helper to parse detected centroids JSON */
static int read_centroids_json(const char* path, starfix_centroid_t* centroids, int max_count) {
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

/* helper to write identified stars to JSON */
static int write_identified_json(const char* path, const starfix_match_t* matches, int count,
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

/* helper to parse identified stars JSON */
static int read_identified_json(const char* path, starfix_match_t* matches, int max_count) {
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

/* helper to write attitude to JSON */
static int write_attitude_json(const char* path, double ra, double dec, double roll,
                               const double R[3][3]) {
    FILE* f = fopen(path, "w");
    if (f == NULL) return -1;

    fprintf(f, "{\n");
    fprintf(f, "  \"ra_deg\": %.8f,\n", ra);
    fprintf(f, "  \"dec_deg\": %.8f,\n", dec);
    fprintf(f, "  \"roll_deg\": %.8f,\n", roll);
    fprintf(f, "  \"R\": [\n");
    fprintf(f, "    [%.16f, %.16f, %.16f],\n", R[0][0], R[0][1], R[0][2]);
    fprintf(f, "    [%.16f, %.16f, %.16f],\n", R[1][0], R[1][1], R[1][2]);
    fprintf(f, "    [%.16f, %.16f, %.16f]\n", R[2][0], R[2][1], R[2][2]);
    fprintf(f, "  ]\n");
    fprintf(f, "}\n");
    fclose(f);
    return 0;
}

/* helper to parse attitude JSON R matrix */
static int read_attitude_json_R(const char* path, double R[3][3]) {
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

    char* p = strstr(buf, "\"R\"");
    if (!p) {
        free(buf);
        return -4;
    }

    int count = 0;
    while (count < 9) {
        p = strpbrk(p, "0123456789.-");
        if (!p) break;
        int r = count / 3;
        int c = count % 3;
        R[r][c] = atof(p);
        /* advance past number */
        while (*p && (strchr("0123456789.eE+-", *p) != NULL)) p++;
        count++;
    }

    free(buf);
    return (count == 9) ? 0 : -5;
}

/* actions */
static int do_centroid() {
    printf("Loading image from data/mock_image.pgm...\n");
    int w = 0, h = 0;
    unsigned char* img = read_pgm("data/mock_image.pgm", &w, &h);
    if (img == NULL) {
        printf("Error loading PGM image.\n");
        return -1;
    }

    starfix_centroid_t centroids[200];
    starfix_telemetry_t telem = {0};
    starfix_status_t count_status =
        starfix_find_centroids(img, w, h, 35, 200, centroids, NULL, 1, &telem);
    int count = (int)telem.num_stars_detected;
    free(img);

    if (count_status != STARFIX_SUCCESS) {
        printf("Centroid extraction failed with code %d\n", (int)count_status);
        return -1;
    }

    printf("Detected %d star centroids.\n", count);
    write_centroids_json("data/detected_centroids.json", centroids, count);
    return 0;
}

static int do_identify(int width, int height, double fov_deg) {
    starfix_catalog_star_t* catalog = NULL;
    starfix_hash_entry_t* hash_table = NULL;
    uint32_t num_stars = 0, num_entries = 0, bin_factor = 0;

    printf("Loading database from data/starfix_db.bin...\n");
    if (starfix_load_db("data/starfix_db.bin", &catalog, &num_stars, &hash_table, &num_entries,
                        &bin_factor) != 0) {
        if (starfix_load_db("../data/starfix_db.bin", &catalog, &num_stars, &hash_table,
                            &num_entries, &bin_factor) != 0) {
            printf("Error loading binary database.\n");
            return -1;
        }
    }

    /* calculate current year dynamically and propagate stellar precession */
    time_t t_now = time(NULL);
    struct tm* tm_now = gmtime(&t_now);
    double current_year = 1900.0 + tm_now->tm_year + (tm_now->tm_yday / 365.25);
    printf("Propagating catalog precession to year %.2f...\n", current_year);
    starfix_propagate_precession(catalog, num_stars, current_year, NULL);

    starfix_centroid_t centroids[100];
    int cent_count = read_centroids_json("data/detected_centroids.json", centroids, 100);
    if (cent_count <= 0) {
        printf("Error loading detected centroids.\n");
        free(catalog);
        free(hash_table);
        return -1;
    }

    printf("Using top %d brightest centroids for pattern matching...\n",
           (cent_count < 8) ? cent_count : 8);
    starfix_match_t matches[20];
    starfix_telemetry_t telem = {0};
    RESET_ARENA();
    starfix_status_t match_status = starfix_identify_stars(
        cent_count, centroids, width, height, fov_deg, num_stars, catalog, num_entries, hash_table,
        bin_factor, 20, matches, &arena, &telem, NULL);
    int match_count = (int)telem.identify_matches;

    if (match_status != STARFIX_SUCCESS) {
        printf("Identification failed with code %d\n", (int)match_status);
        free(catalog);
        free(hash_table);
        return -1;
    }

    printf("Identified Stars:\n");
    printf("Centroid Index  | HIP ID     | RA (deg)   | Dec (deg)\n");
    printf("-----------------------------------------------------\n");
    int i;
    for (i = 0; i < match_count; i++) {
        int c_idx = matches[i].centroid_idx;
        int cat_idx = matches[i].catalog_idx;
        printf("%-15d | %-10d | %10.4f | %10.4f\n", c_idx, catalog[cat_idx].hip,
               starfix_catalog_ra(&catalog[cat_idx]) * 180.0 / M_PI,
               starfix_catalog_dec(&catalog[cat_idx]) * 180.0 / M_PI);
    }

    write_identified_json("data/identified_stars.json", matches, match_count, catalog);
    printf("Saved %d identified stars to data/identified_stars.json\n", match_count);

    free(catalog);
    free(hash_table);
    return 0;
}

static int do_estimate_pose() {
    starfix_match_t matches[30];
    int match_count = read_identified_json("data/identified_stars.json", matches, 30);
    if (match_count < 3) {
        printf("Need at least 3 identified stars for attitude estimation.\n");
        return -1;
    }

    /* reload database just to get catalog vectors */
    starfix_catalog_star_t* catalog = NULL;
    starfix_hash_entry_t* hash_table = NULL;
    uint32_t num_stars = 0, num_entries = 0, bin_factor = 0;
    if (starfix_load_db("data/starfix_db.bin", &catalog, &num_stars, &hash_table, &num_entries,
                        &bin_factor) != 0) {
        printf("Error loading binary database.\n");
        return -1;
    }

    /* calculate current year dynamically and propagate stellar precession */
    time_t t_now = time(NULL);
    struct tm* tm_now = gmtime(&t_now);
    double current_year = 1900.0 + tm_now->tm_year + (tm_now->tm_yday / 365.25);
    printf("Propagating catalog precession to year %.2f...\n", current_year);
    starfix_propagate_precession(catalog, num_stars, current_year, NULL);

    /* load detected centroids */
    starfix_centroid_t centroids[100];
    int cent_count = read_centroids_json("data/detected_centroids.json", centroids, 100);
    if (cent_count <= 0) {
        free(catalog);
        free(hash_table);
        return -1;
    }

    double fov_rad = 12.0 * M_PI / 180.0;
    double focal_length = 1024.0 / (2.0 * tan(fov_rad / 2.0));

    starfix_vector3_t* cam_vectors =
        (starfix_vector3_t*)malloc((size_t)match_count * sizeof(starfix_vector3_t));
    starfix_vector3_t* cat_vectors =
        (starfix_vector3_t*)malloc((size_t)match_count * sizeof(starfix_vector3_t));

    int i;
    for (i = 0; i < match_count; i++) {
        int c_idx = matches[i].centroid_idx;
        int cat_idx = matches[i].catalog_idx;

        /* camera vector */
        double cx = (centroids[c_idx].u - 512.0) / focal_length;
        double cy = -(centroids[c_idx].v - 512.0) / focal_length;
        double cz = 1.0;
        double norm = sqrt(cx * cx + cy * cy + cz * cz);
        cam_vectors[i].x = cx / norm;
        cam_vectors[i].y = cy / norm;
        cam_vectors[i].z = cz / norm;

        /* catalog vector */
        double dec_s = starfix_catalog_dec(&catalog[cat_idx]);
        double ra_s = starfix_catalog_ra(&catalog[cat_idx]);
        cat_vectors[i].x = cos(dec_s) * cos(ra_s);
        cat_vectors[i].y = cos(dec_s) * sin(ra_s);
        cat_vectors[i].z = sin(dec_s);
    }

    /* Apply atmospheric refraction correction to observed camera vectors using standard Zenith
     * (0,0,1) */
    starfix_vector3_t z_cam = {0.0, 0.0, 1.0};
    printf("Applying atmospheric refraction correction (temp: 15.0C, press: 1013.25hPa)...\n");
    starfix_correct_refraction(match_count, cam_vectors, &z_cam, 15.0, 1013.25);

    starfix_quaternion_t q_est;
    double R[3][3];
    starfix_telemetry_t telem = {0};
    /* use RANSAC solver with a 0.1 degree threshold to filter out any outliers */
    RESET_ARENA();
    starfix_status_t inlier_status = starfix_solve_attitude_ransac(
        match_count, cam_vectors, cat_vectors, 0.1, &q_est, R, &arena, &telem);
    int num_inliers = (int)telem.ransac_inliers;
    if (inlier_status != STARFIX_SUCCESS || num_inliers < 2) {
        printf("Attitude solver failed (RANSAC found too few inliers: %d).\n", num_inliers);
        free(catalog);
        free(hash_table);
        free(cam_vectors);
        free(cat_vectors);
        return -1;
    }

    /* extract pointing angles */
    double R_flat[9] = {R[0][0], R[0][1], R[0][2], R[1][0], R[1][1],
                        R[1][2], R[2][0], R[2][1], R[2][2]};
    double ra, dec, roll;
    starfix_rotation_to_euler(R_flat, &ra, &dec, &roll);

    printf("\nAttitude Estimation Results:\n");
    printf("-------------------------------------------------------\n");
    printf("Parameter       | Estimated (deg)  | True (deg)\n");
    printf("-------------------------------------------------------\n");
    printf("Right Ascension |       %10.6f |    80.000000\n", ra);
    printf("Declination     |       %10.6f |    37.774900\n", dec);
    printf("Roll            |       %10.6f |     0.000000\n", roll);
    printf("-------------------------------------------------------\n");

    write_attitude_json("data/attitude.json", ra, dec, roll, (const double (*)[3])R);

    free(catalog);
    free(hash_table);
    free(cam_vectors);
    free(cat_vectors);
    return 0;
}

static int do_solve_fix() {
    double R[3][3];
    if (read_attitude_json_R("data/attitude.json", R) != 0) {
        printf("Error reading estimated attitude matrix R.\n");
        return -1;
    }

    /* true position (San Francisco) */
    double true_lat = 37.7749;
    double true_lon = -122.4194;

    /* assumed position (AP) */
    double ap_lat = 37.5;
    double ap_lon = -122.0;
    double gha_aries = 202.4194;

    /* generate 50 mock Zenith vectors in camera frame using true position */
    double true_lat_rad = true_lat * M_PI / 180.0;
    double true_lst_rad = (gha_aries + true_lon) * M_PI / 180.0;
    double true_zenith_cel[3] = {cos(true_lat_rad) * cos(true_lst_rad),
                                 cos(true_lat_rad) * sin(true_lst_rad), sin(true_lat_rad)};

    /* project true zenith to camera frame: zenith_cam = R_true * zenith_cel */
    double zenith_cam[3];
    zenith_cam[0] =
        R[0][0] * true_zenith_cel[0] + R[0][1] * true_zenith_cel[1] + R[0][2] * true_zenith_cel[2];
    zenith_cam[1] =
        R[1][0] * true_zenith_cel[0] + R[1][1] * true_zenith_cel[1] + R[1][2] * true_zenith_cel[2];
    zenith_cam[2] =
        R[2][0] * true_zenith_cel[0] + R[2][1] * true_zenith_cel[1] + R[2][2] * true_zenith_cel[2];

    starfix_vector3_t* zenith_cam_meas = (starfix_vector3_t*)malloc(50 * sizeof(starfix_vector3_t));
    double gha_aries_arr[50];

    srand(42);
    int i;
    for (i = 0; i < 50; i++) {
        /* add small noise (0.0001 rad corresponds to ~20 arcseconds of tilt error) */
        double n_x = (rand() / (double)RAND_MAX - 0.5) * 0.0001;
        double n_y = (rand() / (double)RAND_MAX - 0.5) * 0.0001;
        double n_z = (rand() / (double)RAND_MAX - 0.5) * 0.0001;

        double vx = zenith_cam[0] + n_x;
        double vy = zenith_cam[1] + n_y;
        double vz = zenith_cam[2] + n_z;
        double norm = sqrt(vx * vx + vy * vy + vz * vz);

        zenith_cam_meas[i].x = vx / norm;
        zenith_cam_meas[i].y = vy / norm;
        zenith_cam_meas[i].z = vz / norm;

        gha_aries_arr[i] = gha_aries;
    }

    double solved_lat = 0.0, solved_lon = 0.0;
    starfix_telemetry_t telem = {0};
    RESET_ARENA();
    starfix_status_t status =
        starfix_solve_position(50, zenith_cam_meas, (const double (*)[3])R, gha_aries_arr, ap_lat,
                               ap_lon, &solved_lat, &solved_lon, &arena, &telem);

    if (status != STARFIX_SUCCESS) {
        printf("Position fix solver failed.\n");
        free(zenith_cam_meas);
        return -1;
    }

    /* position error in nautical miles */
    double err_lat = solved_lat - true_lat;
    double err_lon = (solved_lon - true_lon) * cos(true_lat * M_PI / 180.0);
    double dist_err_nm = sqrt(err_lat * err_lat + err_lon * err_lon) * 60.0;

    /* save solved position to JSON */
    FILE* sf = fopen("data/solved_position.json", "w");
    if (sf != NULL) {
        fprintf(sf, "{\n");
        fprintf(sf, "  \"lat\": %.6f,\n", solved_lat);
        fprintf(sf, "  \"lon\": %.6f\n", solved_lon);
        fprintf(sf, "}\n");
        fclose(sf);
    }

    printf("\nCelestial Position Fix Results:\n");
    printf("-----------------------------------------------------------------\n");
    printf("Parameter          | AP (deg)   | Solved (deg) | True (deg)\n");
    printf("-----------------------------------------------------------------\n");
    printf("Latitude           |   %8.4f |     %8.4f |   %8.4f\n", ap_lat, solved_lat, true_lat);
    printf("Longitude          |  %9.4f |    %9.4f |  %9.4f\n", ap_lon, solved_lon, true_lon);
    printf("-----------------------------------------------------------------\n");
    printf("Position Solver Error: %.4f Nautical Miles (nm)\n", dist_err_nm);
    printf("-----------------------------------------------------------------\n");

    free(zenith_cam_meas);
    return 0;
}

static int do_solve_photo(double gha_aries, double g_x, double g_y, double g_z, double ap_lat,
                          double ap_lon, double fov_deg) {
    starfix_catalog_star_t* catalog = NULL;
    starfix_hash_entry_t* hash_table = NULL;
    uint32_t num_stars = 0, num_entries = 0, bin_factor = 0;

    if (starfix_load_db("data/starfix_db.bin", &catalog, &num_stars, &hash_table, &num_entries,
                        &bin_factor) != 0) {
        if (starfix_load_db("../data/starfix_db.bin", &catalog, &num_stars, &hash_table,
                            &num_entries, &bin_factor) != 0) {
            printf("Failed to load db (tried data/starfix_db.bin and ../data/starfix_db.bin).\n");
            return 1;
        }
    }
    starfix_propagate_precession(catalog, num_stars, 2026.56, NULL);

    /* load detected centroids */
    starfix_centroid_t centroids[100];
    int cent_count = read_centroids_json("data/detected_centroids.json", centroids, 100);
    if (cent_count <= 0) {
        cent_count = read_centroids_json("../data/detected_centroids.json", centroids, 100);
    }
    if (cent_count <= 0) {
        printf("Error: Could not load centroids from data/detected_centroids.json\n");
        free(catalog);
        free(hash_table);
        return 1;
    }

    /* load identified stars */
    starfix_match_t matches[100];
    int match_count = read_identified_json("data/identified_stars.json", matches, 100);
    if (match_count <= 0) {
        match_count = read_identified_json("../data/identified_stars.json", matches, 100);
    }
    if (match_count < 3) {
        printf(
            "Error: Need at least 3 identified stars in data/identified_stars.json (found: %d)\n",
            match_count);
        free(catalog);
        free(hash_table);
        return 1;
    }

    double fov_rad = fov_deg * M_PI / 180.0;
    double crop_size = 560.0;
    double f_len = crop_size / (2.0 * tan(fov_rad / 2.0));

    starfix_vector3_t* cam_vecs =
        (starfix_vector3_t*)malloc((size_t)match_count * sizeof(starfix_vector3_t));
    starfix_vector3_t* cel_vecs =
        (starfix_vector3_t*)malloc((size_t)match_count * sizeof(starfix_vector3_t));

    int i;
    for (i = 0; i < match_count; i++) {
        int c_idx = matches[i].centroid_idx;
        int cat_idx = matches[i].catalog_idx;

        /* convert centroids to camera frame unit vectors */
        double cx = (centroids[c_idx].u - crop_size / 2.0) / f_len;
        double cy = (centroids[c_idx].v - crop_size / 2.0) / f_len;
        double cz = 1.0;
        double norm = sqrt(cx * cx + cy * cy + cz * cz);
        cam_vecs[i].x = cx / norm;
        cam_vecs[i].y = cy / norm;
        cam_vecs[i].z = cz / norm;

        /* convert catalog coordinates to celestial frame unit vectors */
        double dec = starfix_catalog_dec(&catalog[cat_idx]);
        double ra = starfix_catalog_ra(&catalog[cat_idx]);
        cel_vecs[i].x = cos(dec) * cos(ra);
        cel_vecs[i].y = cos(dec) * sin(ra);
        cel_vecs[i].z = sin(dec);
    }

    /* estimate camera attitude matrix R */
    double R[3][3];
    starfix_quaternion_t q_est;
    starfix_status_t status = starfix_solve_attitude(match_count, cam_vecs, cel_vecs, &q_est, R);

    /* calculate rms projection error */
    double sum_sq_err = 0.0;
    if (status == STARFIX_SUCCESS) {
        for (i = 0; i < match_count; i++) {
            int c_idx = matches[i].centroid_idx;
            double vx = R[0][0] * cel_vecs[i].x + R[0][1] * cel_vecs[i].y + R[0][2] * cel_vecs[i].z;
            double vy = R[1][0] * cel_vecs[i].x + R[1][1] * cel_vecs[i].y + R[1][2] * cel_vecs[i].z;
            double vz = R[2][0] * cel_vecs[i].x + R[2][1] * cel_vecs[i].y + R[2][2] * cel_vecs[i].z;
            if (vz > 0.0) {
                double u_proj = crop_size / 2.0 + (vx / vz) * f_len;
                double v_proj = crop_size / 2.0 + (vy / vz) * f_len;
                double du = u_proj - centroids[c_idx].u;
                double dv = v_proj - centroids[c_idx].v;
                sum_sq_err += du * du + dv * dv;
            }
        }
        double rms_err = sqrt(sum_sq_err / (double)match_count);
        printf("RMS_ERROR: %f pixels\n", rms_err);
    }

    free(cam_vecs);
    free(cel_vecs);

    if (status != STARFIX_SUCCESS) {
        printf("Attitude estimation failed with status %d\n", status);
        free(catalog);
        free(hash_table);
        return 1;
    }

    /* Convert R to Euler Angles */
    double R_flat[9] = {R[0][0], R[0][1], R[0][2], R[1][0], R[1][1],
                        R[1][2], R[2][0], R[2][1], R[2][2]};
    double ra_deg, dec_deg, roll_deg;
    starfix_rotation_to_euler(R_flat, &ra_deg, &dec_deg, &roll_deg);

    printf("\nEstimated Attitude:\n");
    printf("  Right Ascension : %.4f degrees\n", ra_deg);
    printf("  Declination     : %.4f degrees\n", dec_deg);
    printf("  Roll            : %.4f degrees\n", roll_deg);
    printf("ROTATION_MATRIX: %f %f %f %f %f %f %f %f %f\n", R[0][0], R[0][1], R[0][2], R[1][0],
           R[1][1], R[1][2], R[2][0], R[2][1], R[2][2]);

    /* Solve absolute position using command line parameters */
    double solved_lat = 0.0;
    double solved_lon = 0.0;

    starfix_vector3_t zenith_meas[1];
    zenith_meas[0].x = g_x;
    zenith_meas[0].y = g_y;
    zenith_meas[0].z = g_z;

    RESET_ARENA();
    starfix_telemetry_t telem = {0};
    starfix_status_t solve_status =
        starfix_solve_position(1, zenith_meas, (const double (*)[3])R, &gha_aries, ap_lat, ap_lon,
                               &solved_lat, &solved_lon, &arena, &telem);

    if (solve_status == STARFIX_SUCCESS) {
        printf("\n>>> CELESTIAL POSITION SOLVED! <<<\n");
        printf("  Solved Latitude  : %f degrees North\n", solved_lat);
        printf("  Solved Longitude : %f degrees East\n", solved_lon);
    } else {
        printf("Position fix solver failed with status %d\n", solve_status);
    }

    free(catalog);
    free(hash_table);
    return (solve_status == STARFIX_SUCCESS) ? 0 : 1;
}

static int do_fuse_ekf() {
    printf("Starting EKF Sensor Fusion Simulation...\n");
    starfix_ekf_t ekf;
    starfix_ekf_init(&ekf, 34.0522, -118.2437);

    double true_speed = 12.0;
    double true_heading = 45.0;
    double compass_bias = 4.0;
    double speed_scale = 0.93;
    double dt = 10.0;

    double cur_true_lat = 34.0522;
    double cur_true_lon = -118.2437;

    int step;
    srand(42);
    for (step = 0; step < 4320; step++) {
        double current_heading = ((step * dt) < 10800) ? true_heading : 135.0;

        /* update true location */
        double dist_nm = true_speed * dt / 3600.0;
        cur_true_lat += dist_nm * cos(current_heading * M_PI / 180.0) / 60.0;
        cur_true_lon += dist_nm * sin(current_heading * M_PI / 180.0) /
                        (60.0 * cos(cur_true_lat * M_PI / 180.0));

        /* simulate noisy measurements */
        double m_speed = true_speed * speed_scale + (rand() / (double)RAND_MAX - 0.5) * 0.6;
        double m_heading = current_heading + compass_bias + (rand() / (double)RAND_MAX - 0.5) * 2.0;

        starfix_ekf_predict(&ekf, m_speed, m_heading, dt);

        /* celestial fixes every 20 minutes (120 steps) */
        if ((step + 1) % 120 == 0) {
            double n_lat = (rand() / (double)RAND_MAX - 0.5) * 0.08; /* ~3.5 nm error */
            double n_lon = (rand() / (double)RAND_MAX - 0.5) * 0.08;
            starfix_telemetry_t telem = {0};
            starfix_ekf_correct(&ekf, cur_true_lat + n_lat, cur_true_lon + n_lon, &telem);
        }
    }

    double err_lat = ekf.lat - cur_true_lat;
    double err_lon = (ekf.lon - cur_true_lon) * cos(cur_true_lat * M_PI / 180.0);
    double dist_err_nm = sqrt(err_lat * err_lat + err_lon * err_lon) * 60.0;

    printf("Simulation Finished!\n");
    printf("=================================================================\n");
    printf("Total distance traveled: 144.0 nautical miles\n");
    printf("-----------------------------------------------------------------\n");
    printf("EKF Sensor Fusion Error: %.4f nm\n", dist_err_nm);
    printf("=================================================================\n\n");
    printf("Online Calibration of Sensor Biases:\n");
    printf("-----------------------------------------------------------------\n");
    printf("Parameter            | Estimated Value      | True Value     \n");
    printf("-----------------------------------------------------------------\n");
    printf("Compass Bias (deg)   |           %10.4f |          4.0000\n", ekf.bias * 180.0 / M_PI);
    printf("Speed Scale Factor   |           %10.4f |          0.9300\n", ekf.scale);
    printf("-----------------------------------------------------------------\n");

    return 0;
}

static int do_fuse_graph() {
    printf("Generating Sailboat Path & Sensor Readings...\n");
    int num_nodes = 50;
    int steps_per_node = 50;
    double dt_sec = 10.0;
    double init_lat = 34.0522;
    double init_lon = -118.2437;
    int total_steps = num_nodes * steps_per_node;

    starfix_odometry_t* odom =
        (starfix_odometry_t*)malloc((size_t)total_steps * sizeof(starfix_odometry_t));
    starfix_pose_t* poses =
        (starfix_pose_t*)malloc((size_t)(num_nodes + 1) * sizeof(starfix_pose_t));

    if (odom == NULL || poses == NULL) {
        printf("Memory allocation failed.\n");
        return -1;
    }

    double true_speed = 12.0;
    double true_scale = 0.93;

    srand(42);
    int i;
    for (i = 0; i < total_steps; i++) {
        double current_heading = (i < total_steps / 2) ? 45.0 : 135.0;
        odom[i].speed = true_speed * true_scale + (rand() / (double)RAND_MAX - 0.5) * 0.1;
        odom[i].heading = current_heading + 4.0 + (rand() / (double)RAND_MAX - 0.5) * 0.2;
    }

    double* true_lats = (double*)malloc((size_t)(num_nodes + 1) * sizeof(double));
    double* true_lons = (double*)malloc((size_t)(num_nodes + 1) * sizeof(double));
    true_lats[0] = init_lat;
    true_lons[0] = init_lon;

    double cur_lat = init_lat;
    double cur_lon = init_lon;
    for (i = 1; i <= num_nodes; i++) {
        int j;
        for (j = 0; j < steps_per_node; j++) {
            int s_idx = (i - 1) * steps_per_node + j;
            double current_heading = (s_idx < total_steps / 2) ? 45.0 : 135.0;
            double dist_nm = true_speed * dt_sec / 3600.0;
            double lat_rad = cur_lat * M_PI / 180.0;
            cur_lat += dist_nm * cos(current_heading * M_PI / 180.0) / 60.0;
            cur_lon += dist_nm * sin(current_heading * M_PI / 180.0) / (60.0 * cos(lat_rad));
        }
        true_lats[i] = cur_lat;
        true_lons[i] = cur_lon;
    }

    int num_fixes = 25;
    starfix_fix_t* fixes = (starfix_fix_t*)malloc((size_t)num_fixes * sizeof(starfix_fix_t));
    for (i = 0; i < num_fixes; i++) {
        fixes[i].node_idx = (i + 1) * 2;
        fixes[i].lat = true_lats[(i + 1) * 2] + (rand() / (double)RAND_MAX - 0.5) * 0.001;
        fixes[i].lon = true_lons[(i + 1) * 2] + (rand() / (double)RAND_MAX - 0.5) * 0.001;
    }

    double est_scale = 1.0;
    double est_bias = 0.0;
    starfix_telemetry_t telem = {0};
    RESET_ARENA();
    starfix_status_t status =
        starfix_solve_graph(num_nodes, steps_per_node, dt_sec, init_lat, init_lon, odom, num_fixes,
                            fixes, poses, &est_scale, &est_bias, &arena, &telem);

    if (status == STARFIX_SUCCESS) {
        double err_lat = poses[num_nodes].lat - true_lats[num_nodes];
        double err_lon = (poses[num_nodes].lon - true_lons[num_nodes]) *
                         cos(true_lats[num_nodes] * M_PI / 180.0);
        double dist_err_nm = sqrt(err_lat * err_lat + err_lon * err_lon) * 60.0;

        printf("Factor Graph Optimization Finished!\n");
        printf("=================================================================\n");
        printf("Total distance traveled: 166.7 nautical miles\n");
        printf("-----------------------------------------------------------------\n");
        printf("Factor Graph Smoother Error: %.4f nm\n", dist_err_nm);
        printf("=================================================================\n\n");
        printf("Calibrated Sensor Biases:\n");
        printf("-----------------------------------------------------------------\n");
        printf("Parameter                 | Estimated Value      | True Value     \n");
        printf("-----------------------------------------------------------------\n");
        printf("Compass Bias (deg)        |           %10.4f |          4.0000\n",
               est_bias * 180.0 / M_PI);
        printf("Speed Scale Factor        |           %10.4f |          0.9300\n", est_scale);
        printf("-----------------------------------------------------------------\n");
    } else {
        printf("Factor Graph Optimization failed with code %d\n", status);
    }

    free(odom);
    free(poses);
    free(fixes);
    free(true_lats);
    free(true_lons);
    return 0;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        printf("StarFix Native C Command-Line Interface\n");
        printf("Usage: %s [option]\n", argv[0]);
        printf("Options:\n");
        printf("  --centroid       Extract sub-pixel star centroids from raw PGM image\n");
        printf("  --identify       Match detected centroids to database catalog quads\n");
        printf("  --estimate-pose  Estimate camera attitude pointing matrix R\n");
        printf("  --solve-fix      Solve absolute absolute Latitude/Longitude position fix\n");
        printf(
            "  --solve-photo    Solve terrestrial position from a real photo (uses: <gha_aries> "
            "<g_x> <g_y> <g_z> <ap_lat> <ap_lon>)\n");
        printf("  --fuse-ekf       Run real-time EKF trajectory fusion simulation\n");
        printf("  --fuse-graph     Run Factor Graph batch trajectory smoother simulation\n");
        printf("  --pipeline       Run the entire end-to-end star tracker pipeline natively\n");
        return 1;
    }

    if (strcmp(argv[1], "--centroid") == 0) {
        return do_centroid();
    } else if (strcmp(argv[1], "--identify") == 0) {
        int w = (argc >= 3) ? atoi(argv[2]) : 1024;
        int h = (argc >= 4) ? atoi(argv[3]) : 1024;
        double fov = (argc >= 5) ? atof(argv[4]) : 12.0;
        return do_identify(w, h, fov);
    } else if (strcmp(argv[1], "--estimate-pose") == 0) {
        return do_estimate_pose();
    } else if (strcmp(argv[1], "--solve-fix") == 0) {
        return do_solve_fix();
    } else if (strcmp(argv[1], "--solve-photo") == 0) {
        if (argc < 8) {
            printf(
                "Usage: %s --solve-photo <gha_aries> <g_x> <g_y> <g_z> <ap_lat> <ap_lon> "
                "[fov_deg]\n",
                argv[0]);
            return 1;
        }
        double gha = atof(argv[2]);
        double gx = atof(argv[3]);
        double gy = atof(argv[4]);
        double gz = atof(argv[5]);
        double lat = atof(argv[6]);
        double lon = atof(argv[7]);
        double fov_deg = (argc >= 9) ? atof(argv[8]) : 11.28;
        return do_solve_photo(gha, gx, gy, gz, lat, lon, fov_deg);
    } else if (strcmp(argv[1], "--fuse-ekf") == 0) {
        return do_fuse_ekf();
    } else if (strcmp(argv[1], "--fuse-graph") == 0) {
        return do_fuse_graph();
    } else if (strcmp(argv[1], "--pipeline") == 0) {
        printf("Starting Native C End-to-End StarFix Pipeline...\n");
        if (do_centroid() != 0) return 1;
        if (do_identify(1024, 1024, 12.0) != 0) return 1;
        if (do_estimate_pose() != 0) return 1;
        if (do_solve_fix() != 0) return 1;
        printf("\nPipeline execution completed successfully!\n");
        return 0;
    } else {
        printf("Invalid option: %s\n", argv[1]);
        return 1;
    }
}
