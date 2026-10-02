#ifndef STARFIX_CLI_IO_H
#define STARFIX_CLI_IO_H

#include "starfix_identify.h"

void starfix_rotation_to_euler(const double R[9], double* ra_deg, double* dec_deg,
                               double* roll_deg);
unsigned char* read_pgm(const char* path, int* width, int* height);
int write_centroids_json(const char* path, const starfix_centroid_t* centroids, int count);
int read_centroids_json(const char* path, starfix_centroid_t* centroids, int max_count);
int write_identified_json(const char* path, const starfix_match_t* matches, int count,
                          const starfix_catalog_star_t* catalog);
int read_identified_json(const char* path, starfix_match_t* matches, int max_count);

#endif
