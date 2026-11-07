#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "starfix_arena.h"

static uint8_t mempool[10 * 1024 * 1024];
static starfix_arena_t arena = {mempool, mempool + sizeof(mempool)};

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "starfix_cli_io.h"
#include "starfix_photo.h"
#include "starfix_solver.h"

static int do_estimate_pose(void) {
    starfix_match_t matches[30];
    int match_count = read_identified_json("data/identified_stars.json", matches, 30);
    if (match_count < 3) {
        printf("Need at least 3 identified stars for attitude estimation.\n");
        return -1;
    }

    /* reload database just to get catalog vectors */
    starfix_db_t db = {0};
    if (starfix_load_db("data/starfix_db.bin", &db) != 0) {
        printf("Error loading binary database.\n");
        return -1;
    }

    /* calculate current year dynamically and propagate stellar precession */
    double current_year = current_year_utc();
    printf("Propagating catalog precession to year %.2f...\n", current_year);
    starfix_propagate_astrometry(db.stars, db.num_stars, current_year);

    /* load detected centroids */
    starfix_centroid_t centroids[1000];
    int cent_count = read_centroids_json("data/detected_centroids.json", centroids, 1000);
    if (cent_count <= 0) {
        starfix_free_db(&db);
        return -1;
    }

    starfix_photo_options_t options = {
        .width = 1024,
        .height = 1024,
        .fov_deg = 12,
        .gha_aries = 202.4194,
        .zenith_camera = {0, 0, 1},
    };
    starfix_photo_result_t result;

    starfix_status_t status = starfix_solve_photo(cent_count, centroids, db.num_stars, db.stars,
                                                  match_count, matches, &options, &result, arena);
    starfix_free_db(&db);
    if (status != STARFIX_SUCCESS) {
        fprintf(stderr, "Simulated photo fit failed (status %d)\n", status);
        return 1;
    }
    const double (*R)[3] = (const double (*)[3])result.rotation;

    double ra, dec, roll;
    starfix_rotation_to_euler(R, &ra, &dec, &roll);

    printf("\nAttitude Estimation Results:\n");
    printf("-------------------------------------------------------\n");
    printf("Parameter       | Estimated (deg)  | True (deg)\n");
    printf("-------------------------------------------------------\n");
    printf("Right Ascension |       %10.6f |    80.000000\n", ra);
    printf("Declination     |       %10.6f |    37.774900\n", dec);
    printf("Roll            |       %10.6f |     0.000000\n", roll);
    printf("-------------------------------------------------------\n");

    if (write_attitude_json("data/attitude.json", ra, dec, roll, (const double (*)[3])R) < 0)
        return -1;

    return 0;
}

static int do_solve_fix(void) {
    double R[3][3];
    if (read_attitude_json("data/attitude.json", R) != 0) {
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

    starfix_status_t status = starfix_solve_position(50, zenith_cam_meas, (const double (*)[3])R,
                                                     gha_aries_arr, &solved_lat, &solved_lon);

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
