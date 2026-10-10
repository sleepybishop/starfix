#include "starfix_centroid.h"

#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "starfix_status.h"

#define MAX_BLOB_PIXELS 1000
#define MAX_STAR_PIXELS 300

static int compare_centroids(const void* a, const void* b) {
    const starfix_centroid_t* ca = (const starfix_centroid_t*)a;
    const starfix_centroid_t* cb = (const starfix_centroid_t*)b;
    if (ca->flux < cb->flux) return 1;
    if (ca->flux > cb->flux) return -1;
    return 0;
}

/* blend the center of mass with a parabolic fit around an interior peak */
static void refine_centroid(const unsigned char* image, int width, int height, int x, int y,
                            starfix_centroid_t* centroid) {
    if (x <= 0 || x >= width - 1 || y <= 0 || y >= height - 1) return;
    double center = image[y * width + x];
    double left = image[y * width + x - 1], right = image[y * width + x + 1];
    double top = image[(y - 1) * width + x], bottom = image[(y + 1) * width + x];
    double dx = left - 2.0 * center + right, dy = top - 2.0 * center + bottom;
    if (!(fabs(dx) > 1e-5 && fabs(dy) > 1e-5)) return;
    double u = 0.5 * (left - right) / dx, v = 0.5 * (top - bottom) / dy;
    if (fabs(u) > 0.5 || fabs(v) > 0.5) return;
    centroid->u = 0.5 * centroid->u + 0.5 * ((double)x + u);
    centroid->v = 0.5 * centroid->v + 0.5 * ((double)y + v);
}

int starfix_find_centroids(const unsigned char* image, int width, int height,
                           unsigned char threshold, int max_centroids,
                           starfix_centroid_t* centroids, uint32_t* visited_mask,
                           uint32_t frame_number) {
    if (image == NULL || centroids == NULL) {
        return STARFIX_ERR_NULL_POINTER;
    }
    if (width <= 0 || height <= 0 || max_centroids <= 0) {
        return STARFIX_ERR_INVALID_PARAM;
    }
    if ((size_t)width > SIZE_MAX / (size_t)height || width > INT_MAX / height) {
        return STARFIX_ERR_INVALID_PARAM;
    }
    int count = 0;

    uint32_t* visited = visited_mask;
    if (visited == NULL) {
        visited = (uint32_t*)calloc((size_t)width * (size_t)height, sizeof(uint32_t));
        if (visited == NULL) return STARFIX_ERR_MEMORY;
        frame_number = 1; /* first frame */
    } else if (frame_number == 0) {
        /* clear to a different generation, preserving zero for this frame */
        memset(visited, 0xff, (size_t)width * (size_t)height * sizeof(uint32_t));
    }

    int queue_x[MAX_BLOB_PIXELS];
    int queue_y[MAX_BLOB_PIXELS];

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int idx = y * width + x;

            if (image[idx] <= threshold || visited[idx] == frame_number) continue;
            int head = 0;
            int tail = 1;

            queue_x[0] = x;
            queue_y[0] = y;
            visited[idx] = frame_number;

            double sum_intensity = 0.0;
            double sum_u = 0.0;
            double sum_v = 0.0;
            double peak = 0.0;
            int peak_x = x;
            int peak_y = y;

            while (head < tail) {
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
                    peak_x = cx;
                    peak_y = cy;
                }

                /* check 8-connected neighbors */
                for (int dy = -1; dy <= 1; dy++) {
                    for (int dx = -1; dx <= 1; dx++) {
                        int nx = cx + dx;
                        int ny = cy + dy;

                        if (nx < 0 || nx >= width || ny < 0 || ny >= height) continue;
                        int n_idx = ny * width + nx;
                        if (image[n_idx] <= threshold || visited[n_idx] == frame_number ||
                            tail >= MAX_BLOB_PIXELS)
                            continue;
                        visited[n_idx] = frame_number;
                        queue_x[tail] = nx;
                        queue_y[tail++] = ny;
                    }
                }
            }

            if (sum_intensity <= 0.0 || tail > MAX_STAR_PIXELS || count >= max_centroids) continue;
            starfix_centroid_t* centroid = &centroids[count++];
            *centroid = (starfix_centroid_t){sum_u / sum_intensity, sum_v / sum_intensity, peak,
                                             sum_intensity, tail};
            refine_centroid(image, width, height, peak_x, peak_y, centroid);
        }
    }

    qsort(centroids, (size_t)count, sizeof(*centroids), compare_centroids);

    if (!visited_mask) {
        free(visited);
    } else if (frame_number == 0) {
        /* remove the temporary sentinel before another generation can use it */
        memset(visited, 0, (size_t)width * (size_t)height * sizeof(uint32_t));
    }
    return count;
}

int starfix_calibrate_hot_pixels(const unsigned char* dark_image, int width, int height,
                                 unsigned char threshold, int max_hot, int* hot_u, int* hot_v) {
    if (!dark_image || !hot_u || !hot_v) return STARFIX_ERR_NULL_POINTER;
    if (width <= 0 || height <= 0 || width > INT_MAX / height || max_hot <= 0)
        return STARFIX_ERR_INVALID_PARAM;

    int count = 0;
    for (int v = 0; v < height && count < max_hot; v++) {
        for (int u = 0; u < width && count < max_hot; u++) {
            if (dark_image[v * width + u] <= threshold) continue;
            hot_u[count] = u;
            hot_v[count++] = v;
        }
    }
    return count;
}

starfix_status_t starfix_apply_hot_pixel_mask(unsigned char* image, int width, int height,
                                              int num_hot, const int* hot_u, const int* hot_v) {
    if (!image || !hot_u || !hot_v) return STARFIX_ERR_NULL_POINTER;
    if (width <= 0 || height <= 0 || width > INT_MAX / height || num_hot < 0)
        return STARFIX_ERR_INVALID_PARAM;

    for (int i = 0; i < num_hot; i++) {
        int u = hot_u[i];
        int v = hot_v[i];
        if (u >= 0 && u < width && v >= 0 && v < height) {
            image[v * width + u] = 0;
        }
    }
    return 0;
}
