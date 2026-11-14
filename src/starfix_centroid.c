#include "starfix_centroid.h"

#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define MAX_STAR_PIXELS 1000
#define INITIAL_SPAN_CAPACITY 2048

typedef struct {
    int xl;
    int xr;
    int y;
} span_t;

static int compare_centroids(const void* a, const void* b) {
    const starfix_centroid_t* ca = (const starfix_centroid_t*)a;
    const starfix_centroid_t* cb = (const starfix_centroid_t*)b;
    if (ca->flux < cb->flux) return 1;
    if (ca->flux > cb->flux) return -1;
    return 0;
}

/* blend the center of mass with a parabolic fit around an unsaturated interior peak */
static void refine_centroid(const unsigned char* image, int width, int height, int x, int y,
                            starfix_centroid_t* centroid) {
    if (x <= 0 || x >= width - 1 || y <= 0 || y >= height - 1) return;
    unsigned char c_px = image[y * width + x];
    if (c_px == 255U) return;
    unsigned char l_px = image[y * width + x - 1], r_px = image[y * width + x + 1];
    unsigned char t_px = image[(y - 1) * width + x], b_px = image[(y + 1) * width + x];
    if (c_px == l_px || c_px == r_px || c_px == t_px || c_px == b_px) return;

    double center = (double)c_px;
    double left = (double)l_px, right = (double)r_px;
    double top = (double)t_px, bottom = (double)b_px;
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

    size_t total_px = (size_t)width * (size_t)height;
    uint32_t* visited = visited_mask;
    uint32_t stamp = frame_number;
    int clear_mask_on_exit = 0;

    if (visited == NULL) {
        visited = (uint32_t*)calloc(total_px, sizeof(uint32_t));
        if (visited == NULL) return STARFIX_ERR_MEMORY;
        stamp = 1U;
    } else if (frame_number == 0U) {
        /* reset mask and clear again on exit so every non-zero generation is fresh */
        memset(visited, 0, total_px * sizeof(uint32_t));
        stamp = 1U;
        clear_mask_on_exit = 1;
    }

    span_t span_storage[INITIAL_SPAN_CAPACITY];
    span_t* stack = span_storage;
    int span_cap = INITIAL_SPAN_CAPACITY;
    int count = 0;
    int min_idx = -1;

    for (int y = 0; y < height; y++) {
        const unsigned char* row_img = image + (size_t)y * (size_t)width;
        uint32_t* row_vis = visited + (size_t)y * (size_t)width;
        for (int x = 0; x < width; x++) {
            if (row_img[x] <= threshold || row_vis[x] == stamp) continue;

            double sum_intensity = 0.0;
            double sum_u = 0.0;
            double sum_v = 0.0;
            double peak = 0.0;
            int peak_x = x;
            int peak_y = y;
            int total_pixels = 0;

            /* expand initial horizontal span on row y */
            int xl = x;
            int xr = x;
            while (xr + 1 < width && row_img[xr + 1] > threshold && row_vis[xr + 1] != stamp) xr++;
            for (int cx = xl; cx <= xr; cx++) {
                row_vis[cx] = stamp;
                double val = (double)row_img[cx] - (double)threshold;
                sum_intensity += val;
                sum_u += val * (double)cx;
                sum_v += val * (double)y;
                if (val > peak) {
                    peak = val;
                    peak_x = cx;
                    peak_y = y;
                }
                total_pixels++;
            }

            int sp = 0;
            stack[sp++] = (span_t){xl, xr, y};
            x = xr;

            while (sp > 0) {
                span_t cur = stack[--sp];
                int scan_l = cur.xl > 0 ? cur.xl - 1 : 0;
                int scan_r = cur.xr + 1 < width ? cur.xr + 1 : width - 1;

                for (int dy = -1; dy <= 1; dy += 2) {
                    int ny = cur.y + dy;
                    if (ny < 0 || ny >= height) continue;
                    const unsigned char* nrow_img = image + (size_t)ny * (size_t)width;
                    uint32_t* nrow_vis = visited + (size_t)ny * (size_t)width;

                    for (int nx = scan_l; nx <= scan_r; nx++) {
                        if (nrow_img[nx] <= threshold || nrow_vis[nx] == stamp) continue;
                        int nxl = nx;
                        while (nxl > 0 && nrow_img[nxl - 1] > threshold &&
                               nrow_vis[nxl - 1] != stamp) {
                            nxl--;
                        }
                        int nxr = nx;
                        while (nxr + 1 < width && nrow_img[nxr + 1] > threshold &&
                               nrow_vis[nxr + 1] != stamp) {
                            nxr++;
                        }
                        for (int cx = nxl; cx <= nxr; cx++) {
                            nrow_vis[cx] = stamp;
                            double val = (double)nrow_img[cx] - (double)threshold;
                            sum_intensity += val;
                            sum_u += val * (double)cx;
                            sum_v += val * (double)ny;
                            if (val > peak) {
                                peak = val;
                                peak_x = cx;
                                peak_y = ny;
                            }
                            total_pixels++;
                        }
                        if (sp == span_cap) {
                            if ((size_t)span_cap > SIZE_MAX / (2U * sizeof(*stack)) ||
                                span_cap > INT_MAX / 2)
                                goto no_memory;
                            int next_cap = span_cap * 2;
                            span_t* grown = (span_t*)malloc((size_t)next_cap * sizeof(*stack));
                            if (!grown) goto no_memory;
                            memcpy(grown, stack, (size_t)sp * sizeof(*stack));
                            if (stack != span_storage) free(stack);
                            stack = grown;
                            span_cap = next_cap;
                        }
                        stack[sp++] = (span_t){nxl, nxr, ny};
                        nx = nxr + 1;
                    }
                }
            }

            if (sum_intensity <= 0.0 || total_pixels > MAX_STAR_PIXELS) continue;
            starfix_centroid_t candidate = {sum_u / sum_intensity, sum_v / sum_intensity, peak,
                                            sum_intensity, total_pixels};
            refine_centroid(image, width, height, peak_x, peak_y, &candidate);
            if (count < max_centroids) {
                centroids[count++] = candidate;
                min_idx = -1;
            } else {
                if (min_idx < 0) {
                    min_idx = 0;
                    for (int i = 1; i < count; i++) {
                        if (centroids[i].flux < centroids[min_idx].flux) min_idx = i;
                    }
                }
                if (candidate.flux > centroids[min_idx].flux) {
                    centroids[min_idx] = candidate;
                    min_idx = 0;
                    for (int i = 1; i < count; i++) {
                        if (centroids[i].flux < centroids[min_idx].flux) min_idx = i;
                    }
                }
            }
        }
    }

    qsort(centroids, (size_t)count, sizeof(*centroids), compare_centroids);

    if (clear_mask_on_exit) memset(visited, 0, total_px * sizeof(uint32_t));
    if (stack != span_storage) free(stack);
    if (!visited_mask) free(visited);
    return count;

no_memory:
    if (clear_mask_on_exit) memset(visited, 0, total_px * sizeof(uint32_t));
    if (stack != span_storage) free(stack);
    if (!visited_mask) free(visited);
    return STARFIX_ERR_MEMORY;
}
