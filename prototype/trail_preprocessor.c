#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define MAX_BLOB_PIXELS 50000
#define MAX_TRAILS 1000

typedef struct {
    int x;
    int y;
} Point;

typedef struct {
    Point* pixels;
    int num_pixels;
} StarTrail;

/* Simple PGM Image Reader */
static unsigned char* read_pgm(const char* filename, int* width, int* height) {
    FILE* f = fopen(filename, "rb");
    if (!f) return NULL;

    char magic[3];
    if (fscanf(f, "%2s", magic) != 1 || strcmp(magic, "P5") != 0) {
        fclose(f);
        return NULL;
    }

    int w, h, maxval;
    if (fscanf(f, "%d %d %d", &w, &h, &maxval) != 3) {
        fclose(f);
        return NULL;
    }
    fgetc(f); /* consume newline */

    *width = w;
    *height = h;

    unsigned char* img = malloc((size_t)w * (size_t)h);
    if (img) {
        if (fread(img, 1, (size_t)w * (size_t)h, f) != (size_t)w * (size_t)h) {
            free(img);
            img = NULL;
        }
    }
    fclose(f);
    return img;
}

/* Bresenham's line algorithm to increment accumulator */
static void draw_line_accumulator(int* accum, int width, int height, int x0, int y0, int x1, int y1) {
    int dx = abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
    int dy = -abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
    int err = dx + dy, e2; /* error value e_xy */

    for (;;) {
        if (x0 >= 0 && x0 < width && y0 >= 0 && y0 < height) {
            accum[y0 * width + x0]++;
        }
        if (x0 == x1 && y0 == y1) break;
        e2 = 2 * err;
        if (e2 >= dy) { err += dy; x0 += sx; }
        if (e2 <= dx) { err += dx; y0 += sy; }
    }
}

int main(int argc, char** argv) {
    if (argc != 3) {
        printf("Usage: %s <input.pgm> <output.json>\n", argv[0]);
        return 1;
    }

    int width, height;
    unsigned char* image = read_pgm(argv[1], &width, &height);
    if (!image) {
        printf("Failed to read PGM image %s\n", argv[1]);
        return 1;
    }
    printf("Loaded image %dx%d\n", width, height);

    /* Allocate memory */
    unsigned char* visited = calloc((size_t)width * height, 1);
    int* accum = calloc((size_t)width * height, sizeof(int));
    StarTrail* trails = calloc(MAX_TRAILS, sizeof(StarTrail));
    int* queue_x = malloc(MAX_BLOB_PIXELS * sizeof(int));
    int* queue_y = malloc(MAX_BLOB_PIXELS * sizeof(int));
    
    if (!visited || !accum || !trails || !queue_x || !queue_y) {
        printf("Memory allocation failed\n");
        return 1;
    }

    int num_trails = 0;
    unsigned char threshold = 50;

    /* 1. Extract Trails using BFS Flood-fill */
    printf("Extracting star trails...\n");
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int idx = y * width + x;
            if (image[idx] > threshold && !visited[idx]) {
                int head = 0, tail = 1;
                queue_x[0] = x;
                queue_y[0] = y;
                visited[idx] = 1;

                while (head < tail && tail < MAX_BLOB_PIXELS) {
                    int cx = queue_x[head];
                    int cy = queue_y[head];
                    head++;

                    /* 8-connected neighbors */
                    for (int dy = -1; dy <= 1; dy++) {
                        for (int dx = -1; dx <= 1; dx++) {
                            int nx = cx + dx;
                            int ny = cy + dy;
                            if (nx >= 0 && nx < width && ny >= 0 && ny < height) {
                                int n_idx = ny * width + nx;
                                if (image[n_idx] > threshold && !visited[n_idx]) {
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

                /* If blob is long enough to be a trail (e.g. > 50 pixels) */
                if (tail > 50 && num_trails < MAX_TRAILS) {
                    trails[num_trails].num_pixels = tail;
                    trails[num_trails].pixels = malloc(tail * sizeof(Point));
                    for (int i = 0; i < tail; i++) {
                        trails[num_trails].pixels[i].x = queue_x[i];
                        trails[num_trails].pixels[i].y = queue_y[i];
                    }
                    num_trails++;
                }
            }
        }
    }
    printf("Found %d valid star trails.\n", num_trails);

    /* 2. Vote for Celestial Pole in Accumulator */
    printf("Computing Celestial Pole via perpendicular bisectors...\n");
    for (int t = 0; t < num_trails; t++) {
        StarTrail* trail = &trails[t];
        /* Pick points at 10%, 50%, 90% of the BFS queue (which correlates roughly to spatial spread) */
        int p1_idx = trail->num_pixels / 10;
        int p2_idx = trail->num_pixels - (trail->num_pixels / 10);
        
        Point p1 = trail->pixels[p1_idx];
        Point p2 = trail->pixels[p2_idx];

        double mx = (p1.x + p2.x) / 2.0;
        double my = (p1.y + p2.y) / 2.0;
        
        double dx = p2.x - p1.x;
        double dy = p2.y - p1.y;
        
        /* Perpendicular vector */
        double px = -dy;
        double py = dx;
        
        /* Normalize and extend to cover image */
        double norm = sqrt(px*px + py*py);
        if (norm > 0) {
            px /= norm;
            py /= norm;
            
            int ext = width > height ? width : height;
            int start_x = (int)(mx - px * ext);
            int start_y = (int)(my - py * ext);
            int end_x = (int)(mx + px * ext);
            int end_y = (int)(my + py * ext);
            
            draw_line_accumulator(accum, width, height, start_x, start_y, end_x, end_y);
        }
    }

    /* Find max vote in accumulator */
    int max_votes = 0;
    int pole_x = width / 2;
    int pole_y = height / 2;
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            if (accum[y * width + x] > max_votes) {
                max_votes = accum[y * width + x];
                pole_x = x;
                pole_y = y;
            }
        }
    }
    printf("Celestial Pole detected at (%d, %d) with %d votes.\n", pole_x, pole_y, max_votes);

    /* 3. Compute Synthetic Centroids (Angular Midpoints) */
    printf("Extracting synthetic centroids...\n");
    FILE* out_f = fopen(argv[2], "w");
    fprintf(out_f, "[\n");
    
    int valid_centroids = 0;
    for (int t = 0; t < num_trails; t++) {
        StarTrail* trail = &trails[t];
        double sum_r = 0;
        double sum_vx = 0, sum_vy = 0;
        
        for (int i = 0; i < trail->num_pixels; i++) {
            double rx = trail->pixels[i].x - pole_x;
            double ry = trail->pixels[i].y - pole_y;
            double r = sqrt(rx*rx + ry*ry);
            
            sum_r += r;
            if (r > 0) {
                sum_vx += rx / r;
                sum_vy += ry / r;
            }
        }
        
        double avg_r = sum_r / trail->num_pixels;
        double avg_theta = atan2(sum_vy, sum_vx);
        
        double synth_x = pole_x + avg_r * cos(avg_theta);
        double synth_y = pole_y + avg_r * sin(avg_theta);
        
        if (valid_centroids > 0) fprintf(out_f, ",\n");
        fprintf(out_f, "  {\"u\": %.2f, \"v\": %.2f}", synth_x, synth_y);
        valid_centroids++;
    }
    fprintf(out_f, "\n]\n");
    fclose(out_f);
    printf("Saved %d synthetic centroids to %s\n", valid_centroids, argv[2]);

    /* Cleanup */
    for(int t=0; t<num_trails; t++) {
        free(trails[t].pixels);
    }
    free(trails);
    free(image);
    free(visited);
    free(accum);
    free(queue_x);
    free(queue_y);
    
    return 0;
}
