#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "starfix_arena.h"

static uint8_t mempool[10 * 1024 * 1024];
static starfix_arena_t arena;
#define RESET_ARENA() (arena.beg = mempool, arena.end = mempool + sizeof(mempool))

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "starfix_attitude.h"
#include "starfix_cli_io.h"
#include "starfix_solver.h"

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
    starfix_propagate_astrometry(catalog, num_stars, current_year, NULL);

    /* load detected centroids */
    starfix_centroid_t centroids[1000];
    int cent_count = read_centroids_json("data/detected_centroids.json", centroids, 1000);
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
    if (cam_vectors == NULL || cat_vectors == NULL) {
        free(catalog);
        free(hash_table);
        free(cam_vectors);
        free(cat_vectors);
        return -1;
    }

    int i;
    for (i = 0; i < match_count; i++) {
        int c_idx = matches[i].centroid_idx;
        int cat_idx = matches[i].catalog_idx;
        if (c_idx < 0 || c_idx >= cent_count || cat_idx < 0 || (uint32_t)cat_idx >= num_stars) {
            free(catalog);
            free(hash_table);
            free(cam_vectors);
            free(cat_vectors);
            return -1;
        }

        /* camera vector */
        double cx = (centroids[c_idx].u - 512.0) / focal_length;
        double cy = (centroids[c_idx].v - 512.0) / focal_length;
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
        starfix_solve_position(50, zenith_cam_meas, (const double (*)[3])R, gha_aries_arr,
                               &solved_lat, &solved_lon, &telem);

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
    printf("Parameter          | Solved (deg) | True (deg)\n");
    printf("-----------------------------------------------------------------\n");
    printf("Latitude           |     %8.4f |   %8.4f\n", solved_lat, true_lat);
    printf("Longitude          |    %9.4f |  %9.4f\n", solved_lon, true_lon);
    printf("-----------------------------------------------------------------\n");
    printf("Position Solver Error: %.4f Nautical Miles (nm)\n", dist_err_nm);
    printf("-----------------------------------------------------------------\n");

    free(zenith_cam_meas);
    return 0;
}

int main(int argc, char** argv) {
    if (argc == 2 && strcmp(argv[1], "--estimate-pose") == 0) return do_estimate_pose();
    if (argc == 2 && strcmp(argv[1], "--solve-fix") == 0) return do_solve_fix();
    fprintf(stderr, "Usage: %s --estimate-pose | --solve-fix\n", argv[0]);
    return 1;
}
