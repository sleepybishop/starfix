#include <assert.h>
#include <limits.h>
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

    const char* invalid[] = {"\"garbage\"", "null",    "true", "",   "1junk",
                             "1e9999",      "1e-9999", "nan",  "inf"};
    char record[256];
    for (size_t i = 0; i < sizeof(invalid) / sizeof(invalid[0]); i++) {
        snprintf(record, sizeof(record), "[{\"u\": %s, \"v\": 2}]", invalid[i]);
        write_text(path, record);
        assert(read_centroids_json(path, output, 2) == STARFIX_ERR_DATA_CORRUPTION);
        snprintf(record, sizeof(record), "[{\"u\": 1, \"v\": 2, \"flux\": %s}]", invalid[i]);
        write_text(path, record);
        assert(read_centroids_json(path, output, 2) == STARFIX_ERR_DATA_CORRUPTION);
        snprintf(record, sizeof(record), "[{\"centroid_idx\": %s, \"catalog_idx\": 0}]",
                 invalid[i]);
        write_text(path, record);
        assert(read_identified_json(path, readback, 2) == STARFIX_ERR_DATA_CORRUPTION);
    }
    write_text(path, "[{\"centroid_idx\": 1.5, \"catalog_idx\": 0}]");
    assert(read_identified_json(path, readback, 2) == STARFIX_ERR_DATA_CORRUPTION);
    snprintf(record, sizeof(record), "[{\"centroid_idx\": %llu, \"catalog_idx\": 0}]",
             (unsigned long long)INT_MAX + 1);
    write_text(path, record);
    assert(read_identified_json(path, readback, 2) == STARFIX_ERR_DATA_CORRUPTION);
    write_text(path, "[{\"u\": 1.25e1, \"v\": 2.5e-1, \"flux\": 3e2 }]");
    assert(read_centroids_json(path, output, 2) == 1);
    assert(fabs(output[0].u - 12.5) < 1e-12 && fabs(output[0].v - 0.25) < 1e-12);
    assert(fabs(output[0].flux - 300) < 1e-12);
#ifdef __linux__
    assert(write_centroids_json("/dev/full", input, 2) < 0);
    assert(write_identified_json("/dev/full", matches, 2, catalog) < 0);
    starfix_centroid_t many[128] = {{0}};
    assert(write_centroids_json("/dev/full", many, 128) < 0);
#endif

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
