#include "starfix_centroid.h"
#include "starfix_status.h"

#include <stdlib.h>
#include <string.h>

#define MAX_BLOB_PIXELS 1000
#define MAX_STAR_PIXELS 300
#define MAX_BFS_ITERS 1500


/* helper to sort centroids in place (flux descending) */
static int compare_centroids(const void* a, const void* b) {
    const starfix_centroid_t* ca = (const starfix_centroid_t*)a;
    const starfix_centroid_t* cb = (const starfix_centroid_t*)b;
    if (ca->flux < cb->flux) return 1;
    if (ca->flux > cb->flux) return -1;
    return 0;
}

static void sort_centroids(starfix_centroid_t* arr, int len) {
    qsort(arr, (size_t)len, sizeof(starfix_centroid_t), compare_centroids);
}

starfix_status_t starfix_find_centroids(const unsigned char* image, int width, int height,
                           unsigned char threshold, int max_centroids,
                           starfix_centroid_t* centroids, starfix_telemetry_t* telem) {
    if (image == NULL || centroids == NULL) {
        return STARFIX_ERR_CENTROID_BUFFER_FULL;
    }
    if (width <= 0 || height <= 0 || max_centroids <= 0) {
        return -3;
    }
    int count = 0;
    int x, y, dx, dy;

    /* allocate visited mask (0 for unvisited, 1 for visited) */
    unsigned char* visited =
        (unsigned char*)calloc((size_t)width * (size_t)height, sizeof(unsigned char));
    if (visited == NULL) {
        return STARFIX_ERR_NULL_POINTER;
    }

    /* static buffers for the flood fill bfs queue to avoid dynamic allocation in loop */
    int queue_x[MAX_BLOB_PIXELS];
    int queue_y[MAX_BLOB_PIXELS];

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            int idx = y * width + x;

            /* start a new connected component if we find a bright unvisited pixel */
            if (image[idx] > threshold && !visited[idx]) {
                int head = 0;
                int tail = 1;

                queue_x[0] = x;
                queue_y[0] = y;
                visited[idx] = 1;

                double sum_intensity = 0.0;
                double sum_u = 0.0;
                double sum_v = 0.0;
                double peak = 0.0;

                int bfs_iters = 0;
                /* run bfs flood fill */
                while (head < tail && bfs_iters < MAX_BFS_ITERS) {
                    bfs_iters++;
                    int cx = queue_x[head];
                    int cy = queue_y[head];
                    head++;

                    double val = (double)image[cy * width + cx] - (double)threshold;
                    if (val < 0.0) val = 0.0;
                    sum_intensity += val;
                    sum_u += val * cx;
                    sum_v += val * cy;
                    if (val > peak) {
                        peak = val;
                    }

                    /* check 8-connected neighbors */
                    for (dy = -1; dy <= 1; dy++) {
                        for (dx = -1; dx <= 1; dx++) {
                            int nx = cx + dx;
                            int ny = cy + dy;

                            if (nx >= 0 && nx < width && ny >= 0 && ny < height) {
                                int n_idx = ny * width + nx;
                                if (image[n_idx] > threshold && !visited[n_idx]) {
                                    /* check queue capacity limit to avoid buffer overflow */
                                    if (tail < MAX_BLOB_PIXELS) {
                                        visited[n_idx] = 1;
                                        queue_x[tail] = nx;
                                        queue_y[tail] = ny;
                                        tail++;
                                    }
                                }
                            }
                        }
                    }
                }

                /* validate star blob size */
                if (sum_intensity > 0.0 && tail <= MAX_STAR_PIXELS) {
                    if (count < max_centroids) {
                        centroids[count].u = sum_u / sum_intensity;
                        centroids[count].v = sum_v / sum_intensity;
                        centroids[count].peak = peak;
                        centroids[count].flux = sum_intensity;
                        centroids[count].num_pixels = tail;
                        count++;
                    }
                }
            }
        }
    }

    /* sort detected centroids by integrated flux (brightest first) */
    sort_centroids(centroids, count);

    if (telem) { telem->num_stars_detected = (uint32_t)count; }
    return STARFIX_SUCCESS;
}

int starfix_calibrate_hot_pixels(const unsigned char* dark_image, int width, int height,
                                 unsigned char threshold, int max_hot, int* hot_u, int* hot_v) {
    if (dark_image == NULL || hot_u == NULL || hot_v == NULL || width <= 0 || height <= 0 ||
        max_hot <= 0) {
        return STARFIX_ERR_NULL_POINTER;
    }

    int count = 0;
    int u, v;
    for (v = 0; v < height; v++) {
        for (u = 0; u < width; u++) {
            if (dark_image[v * width + u] > threshold) {
                if (count < max_hot) {
                    hot_u[count] = u;
                    hot_v[count] = v;
                    count++;
                } else {
                    return count; /* reached capacity limit */
                }
            }
        }
    }
    return count;
}

int starfix_apply_hot_pixel_mask(unsigned char* image, int width, int height, int num_hot,
                                 const int* hot_u, const int* hot_v) {
    if (image == NULL || hot_u == NULL || hot_v == NULL || width <= 0 || height <= 0) {
        return STARFIX_ERR_NULL_POINTER;
    }
    if (num_hot < 0) {
        return STARFIX_ERR_CENTROID_BUFFER_FULL;
    }

    int i;
    for (i = 0; i < num_hot; i++) {
        int u = hot_u[i];
        int v = hot_v[i];
        if (u >= 0 && u < width && v >= 0 && v < height) {
            image[v * width + u] = 0;
        }
    }
    return 0;
}
