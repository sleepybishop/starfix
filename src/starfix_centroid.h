#ifndef STARFIX_CENTROID_H
#define STARFIX_CENTROID_H

#include <stdint.h>

#include "starfix_status.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double u;       /* sub-pixel horizontal coordinate */
    double v;       /* sub-pixel vertical coordinate */
    double peak;    /* peak intensity above threshold */
    double flux;    /* summed intensity above threshold */
    int num_pixels; /* component size in pixels */
} starfix_centroid_t;

/* find flux-weighted centroids in a row-major byte image, subtracting threshold from intensity.
 * keep the first max_centroids components, sorted by flux; visited_mask has width * height entries.
 * zero-initialize optional visited_mask; increment frame_number each call (zero resets it).
 * returns the count (up to max_centroids) or a negative error */
int starfix_find_centroids(const unsigned char* image, int width, int height,
                           unsigned char threshold, int max_centroids,
                           starfix_centroid_t* centroids, uint32_t* visited_mask,
                           uint32_t frame_number);

/* find hot pixels; hot_u/hot_v hold max_hot entries; returns count or negative error */
int starfix_calibrate_hot_pixels(const unsigned char* dark_image, int width, int height,
                                 unsigned char threshold, int max_hot, int* hot_u, int* hot_v);

/* zero hot_u/hot_v pixels in the row-major image; returns zero or a negative error */
starfix_status_t starfix_apply_hot_pixel_mask(unsigned char* image, int width, int height,
                                              int num_hot, const int* hot_u, const int* hot_v);

#ifdef __cplusplus
}
#endif
#endif /* STARFIX_CENTROID_H */
