#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "starfix_cli_io.h"

static void write_text(const char* path, const char* text) {
    FILE* f = fopen(path, "w");
    assert(f);
    assert(fputs(text, f) >= 0);
    assert(fclose(f) == 0);
}

int main(int argc, char** argv) {
    assert(argc == 2);
    const char* path = argv[1];
    starfix_centroid_t input[2] = {{1.25, 2.5, 0, 30, 1}, {4.5, 6.75, 0, 40, 1}};
    starfix_centroid_t output[2] = {{0}, {.u = -99}};
    assert(write_centroids_json(path, input, 2) == 0);
    assert(read_centroids_json(path, output, 1) == 1);
    assert(fabs(output[0].u - 1.25) < 1e-12 && fabs(output[0].v - 2.5) < 1e-12);
    assert(fabs(output[0].flux - 30) < 1e-12 && fabs(output[1].u + 99) < 1e-12);
    assert(read_centroids_json(path, output, 2) == 2);
    assert(fabs(output[1].u - 4.5) < 1e-12 && fabs(output[1].v - 6.75) < 1e-12);

    /* missing fields must not borrow values from the following record */
    write_text(path, "[{\"u\"}, {\"u\": 7, \"v\": 8}, {\"u\": -1, \"v\": 9}, {\"u\": 3}]");
    assert(read_centroids_json(path, output, 2) == 1);
    assert(fabs(output[0].u - 7) < 1e-12 && fabs(output[0].v - 8) < 1e-12);
    assert(fabs(output[0].flux) < 1e-12);
    write_text(path, "[{\"u\": 1, \"v\": 2}, {\"u\": 3");
    assert(read_centroids_json(path, output, 2) == 1);
    assert(read_centroids_json(path, output, 0) == 0);

    starfix_catalog_star_t catalog[2] = {{.hip = 10}, {.hip = 20}};
    starfix_match_t matches[2] = {{3, 1}, {5, 0}}, readback[2] = {{0}, {-99, -99}};
    assert(write_identified_json(path, matches, 2, catalog) == 0);
    assert(read_identified_json(path, readback, 1) == 1);
    assert(readback[0].centroid_idx == 3 && readback[0].catalog_idx == 1);
    assert(readback[1].centroid_idx == -99 && readback[1].catalog_idx == -99);
    assert(read_identified_json(path, readback, 2) == 2);
    assert(memcmp(matches, readback, sizeof(matches)) == 0);
    write_text(path, "[{\"centroid_idx\"}, {\"centroid_idx\": 4, \"catalog_idx\": 0}]");
    assert(read_identified_json(path, readback, 2) == 1);
    assert(readback[0].centroid_idx == 4 && readback[0].catalog_idx == 0);
    write_text(path, "");
    assert(read_centroids_json(path, output, 2) == 0);
    assert(read_identified_json(path, readback, 2) == 0);
    assert(remove(path) == 0);
    assert(read_centroids_json(path, output, 2) < 0);
    assert(read_identified_json(path, readback, 2) < 0);

    const double rotation[3][3] = {{0, 1, 0}, {0, 0, 1}, {1, 0, 0}};
    double ra, dec, roll;
    starfix_rotation_to_euler(rotation, &ra, &dec, &roll);
    assert(fabs(ra) < 1e-12 && fabs(dec) < 1e-12 && fabs(roll) < 1e-12);
    puts("CLI I/O checks passed");
    return 0;
}
