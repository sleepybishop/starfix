#include "starfix_centroid.h"

#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "starfix_status.h"

#define INITIAL_QUEUE_PIXELS 1000
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
        /* reset the mask; any non-zero frame_number on the next call is fresh */
        memset(visited, 0, (size_t)width * (size_t)height * sizeof(uint32_t));
        frame_number = UINT32_MAX;
    }

    int queue_storage[INITIAL_QUEUE_PIXELS];
    int* queue = queue_storage;
    int queue_capacity = INITIAL_QUEUE_PIXELS;

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int idx = y * width + x;

            if (image[idx] <= threshold || visited[idx] == frame_number) continue;
            int head = 0;
            int tail = 1;
            int pending = 1;
            int total_pixels = 1;

            queue[0] = idx;
            visited[idx] = frame_number;

            double sum_intensity = 0.0;
            double sum_u = 0.0;
            double sum_v = 0.0;
            double peak = 0.0;
            int peak_x = x;
            int peak_y = y;

            while (pending > 0) {
                int current = queue[head];
                int cx = current % width;
                int cy = current / width;
                head = (head + 1) % queue_capacity;
                pending--;

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
                        if (image[n_idx] <= threshold || visited[n_idx] == frame_number) continue;
                        if (pending == queue_capacity) {
                            /* keep every frontier pixel so large blobs do not fragment */
                            int pixels = width * height;
                            int capacity =
                                queue_capacity > pixels / 2 ? pixels : 2 * queue_capacity;
                            if ((size_t)capacity > SIZE_MAX / sizeof(*queue)) goto no_memory;
                            int* grown = malloc((size_t)capacity * sizeof(*queue));
                            if (!grown) goto no_memory;
                            size_t first = (size_t)(queue_capacity - head);
                            memcpy(grown, queue + head, first * sizeof(*queue));
                            memcpy(grown + first, queue, (size_t)head * sizeof(*queue));
                            if (queue != queue_storage) free(queue);
                            queue = grown;
                            queue_capacity = capacity;
                            head = 0;
                            tail = pending;
                        }
                        visited[n_idx] = frame_number;
                        total_pixels++;
                        queue[tail] = n_idx;
                        tail = (tail + 1) % queue_capacity;
                        pending++;
                    }
                }
            }

            if (sum_intensity <= 0.0 || total_pixels > MAX_STAR_PIXELS) continue;
            starfix_centroid_t candidate = {sum_u / sum_intensity, sum_v / sum_intensity, peak,
                                            sum_intensity, total_pixels};
            refine_centroid(image, width, height, peak_x, peak_y, &candidate);
            if (count < max_centroids) {
                centroids[count++] = candidate;
            } else {
                int min_idx = 0;
                for (int i = 1; i < count; i++) {
                    if (centroids[i].flux < centroids[min_idx].flux) min_idx = i;
                }
                if (candidate.flux > centroids[min_idx].flux) centroids[min_idx] = candidate;
            }
        }
    }

    qsort(centroids, (size_t)count, sizeof(*centroids), compare_centroids);

    if (queue != queue_storage) free(queue);
    if (!visited_mask) free(visited);
    return count;

no_memory:
    if (queue != queue_storage) free(queue);
    if (!visited_mask) free(visited);
    return STARFIX_ERR_MEMORY;
}
