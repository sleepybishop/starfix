#ifndef STARFIX_CENTROID_H
#define STARFIX_CENTROID_H

#include "starfix_status.h"

typedef struct {
    double u;       /* sub-pixel horizontal coordinate */
    double v;       /* sub-pixel vertical coordinate */
    double peak;    /* peak intensity value */
    double flux;    /* integrated intensity (brightness) */
    int num_pixels; /* component size in pixels */
} starfix_centroid_t;

/* find flux-weighted centroids in a row-major byte image, subtracting threshold from intensity.
 * centroids holds max_centroids; optional visited_mask has width * height generation entries.
 * frame_number zero clears/restarts the mask.
 * returns the count (up to max_centroids) or a negative error */
int starfix_find_centroids(const unsigned char* image, int width, int height,
                           unsigned char threshold, int max_centroids,
                           starfix_centroid_t* centroids, uint32_t* visited_mask,
                           uint32_t frame_number);

/* find hot pixels; hot_u/hot_v hold max_hot entries; returns count or negative error */
int starfix_calibrate_hot_pixels(const unsigned char* dark_image, int width, int height,
                                 unsigned char threshold, int max_hot, int* hot_u, int* hot_v);

/* zero hot_u/hot_v pixels in the row-major image; returns zero or a negative error */
int starfix_apply_hot_pixel_mask(unsigned char* image, int width, int height, int num_hot,
                                 const int* hot_u, const int* hot_v);

#endif /* STARFIX_CENTROID_H */
