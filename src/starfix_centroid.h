#include "starfix_status.h"

#ifndef STARFIX_CENTROID_H
#define STARFIX_CENTROID_H

/* starfix_centroid: C API for sub-pixel star centroiding. detects star blobs
   in a raw grayscale image buffer using a threshold and 8-connected flood fill,
   computing their sub-pixel weighted Centers of Mass. */

typedef struct {
    double u;       /* sub-pixel horizontal coordinate */
    double v;       /* sub-pixel vertical coordinate */
    double peak;    /* peak intensity value */
    double flux;    /* integrated intensity (brightness) */
    int num_pixels; /* component size in pixels */
} starfix_centroid_t;

/* finds centroids in a raw grayscale image buffer:
   - image: flat 1D array of pixel bytes (row-major order)
   - width: image width in pixels
   - height: image height in pixels
   - threshold: background subtraction intensity threshold
   - max_centroids: capacity of the output centroids array
   - centroids: output array to be populated with detected stars
   - returns: number of centroids detected on success, or a negative value on error */
starfix_status_t starfix_find_centroids(const unsigned char* image, int width, int height,
                           unsigned char threshold, int max_centroids,
                           starfix_centroid_t* centroids, starfix_telemetry_t* telem);

/* scans a dark image (with lens cap on) to identify hot pixels above threshold,
   populating the output mask arrays:
   - dark_image: raw dark frame pixels
   - width, height: image dimensions
   - threshold: intensity above which a pixel is considered "hot"
   - max_hot: capacity of the output coordinate arrays
   - hot_u, hot_v: output arrays for the detected hot pixel coordinates
   - returns: number of hot pixels detected, or negative error code */
int starfix_calibrate_hot_pixels(const unsigned char* dark_image, int width, int height,
                                 unsigned char threshold, int max_hot, int* hot_u, int* hot_v);

/* zeroes out the masked hot pixels in the target image buffer:
   - image: target raw image buffer (modified in-place)
   - width, height: image dimensions
   - num_hot: number of hot pixels in the mask
   - hot_u, hot_v: coordinates of the hot pixels to mask
   - returns: 0 on success, negative error code */
int starfix_apply_hot_pixel_mask(unsigned char* image, int width, int height, int num_hot,
                                 const int* hot_u, const int* hot_v);

#endif /* STARFIX_CENTROID_H */
