#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "starfix_centroid.h"
#include "starfix_status.h"

static int large_component_checks(void) {
    const int width = 1024, height = 1024;
    unsigned char* image = calloc((size_t)width * (size_t)height, sizeof(*image));
    uint32_t* visited = calloc((size_t)width * (size_t)height, sizeof(*visited));
    if (!image || !visited) {
        free(image);
        free(visited);
        return 0;
    }
    starfix_centroid_t centroids[4];
    /* 3301-pixel ellipse must not leave a false star after queue wrap */
    for (int y = 480; y <= 544; y++) {
        for (int x = 476; x <= 548; x++) {
            double dx = (x - 512) / 34.27, dy = (y - 512) / 30.69;
            if (dx * dx + dy * dy <= 1.0) image[y * width + x] = 255;
        }
    }
    image[20 * width + 20] = 200;
    int count = starfix_find_centroids(image, width, height, 35, 4, centroids, visited, 1);
    int passed =
        count == 1 && fabs(centroids[0].u - 20.0) < 1e-6 && fabs(centroids[0].v - 20.0) < 1e-6;

    /* broad flood frontier exceeds 1000 queued pixels */
    memset(image, 0, (size_t)width * (size_t)height);
    for (int y = 100; y < 900; y++) memset(image + y * width + 100, 255, 800);
    for (int y = 900; y < 1000; y++) image[y * width + 899] = 255;
    image[20 * width + 20] = 200;
    image[1000 * width + 1000] = 150;
    count = starfix_find_centroids(image, width, height, 35, 4, centroids, visited, 2);
    passed &= count == 2 && fabs(centroids[0].u - 20.0) < 1e-6 &&
              fabs(centroids[0].v - 20.0) < 1e-6 && fabs(centroids[1].u - 1000.0) < 1e-6 &&
              fabs(centroids[1].v - 1000.0) < 1e-6;
    free(image);
    free(visited);
    return passed;
}

int main(void) {
    int failures = 0;
    printf("--- Running Centroid Comprehensive Tests ---\n");

    int width = 128;
    int height = 128;
    unsigned char* image =
        (unsigned char*)malloc((size_t)width * (size_t)height * sizeof(unsigned char));
    starfix_centroid_t centroids[20];

    if (image == NULL) {
        printf("Memory allocation failed.\n");
        return -1;
    }

    /* 1. NULL Pointer & Invalid Size Checks */
    int ret1 = starfix_find_centroids(NULL, width, height, 35, 10, centroids, NULL, 1);
    int ret2 = starfix_find_centroids(image, 0, height, 35, 10, centroids, NULL, 1);
    int ret3 = starfix_find_centroids(image, width, height, 35, 0, centroids, NULL, 1);

    if (ret1 < 0 && ret2 < 0 && ret3 < 0) {
        printf("Test 1: Boundary check... Passed\n");
    } else {
        failures++;
        printf("Test 1: Boundary check... FAILED (got %d, %d, %d)\n", ret1, ret2, ret3);
    }

    /* 2. Completely Black Image (Thresholding check) */
    int i;
    for (i = 0; i < width * height; i++) {
        image[i] = 10;
    }
    int count_black = starfix_find_centroids(image, width, height, 35, 10, centroids, NULL, 1);
    if (count_black == 0) {
        printf("Test 2: Black image (no stars)... Passed\n");
    } else {
        failures++;
        printf("Test 2: Black image (no stars)... FAILED (detected %d centroids)\n", count_black);
    }

    /* 3. Normal Happy Path Star Detection */
    for (i = 0; i < width * height; i++) {
        image[i] = 15;
    }
    /* draw star 0: center (30.4, 50.8) */
    double cx0 = 30.4, cy0 = 50.8;
    int x, y;
    for (y = 47; y <= 54; y++) {
        for (x = 27; x <= 34; x++) {
            double dist_sq = (x - cx0) * (x - cx0) + (y - cy0) * (y - cy0);
            if (dist_sq <= 9.0) {
                double val = 200.0 * exp(-dist_sq / 2.0);
                unsigned char pix = (unsigned char)(15.0 + val);
                if (pix > image[y * width + x]) {
                    image[y * width + x] = pix;
                }
            }
        }
    }
    /* draw star 1: center (85.2, 20.3) */
    double cx1 = 85.2, cy1 = 20.3;
    for (y = 18; y <= 23; y++) {
        for (x = 83; x <= 88; x++) {
            double dist_sq = (x - cx1) * (x - cx1) + (y - cy1) * (y - cy1);
            if (dist_sq <= 4.0) {
                double val = 150.0 * exp(-dist_sq / 1.0);
                unsigned char pix = (unsigned char)(15.0 + val);
                if (pix > image[y * width + x]) {
                    image[y * width + x] = pix;
                }
            }
        }
    }

    int count_happy = starfix_find_centroids(image, width, height, 35, 10, centroids, NULL, 1);
    if (count_happy == 2) {
        double err0 = sqrt(pow(centroids[0].u - 30.4, 2) + pow(centroids[0].v - 50.8, 2));
        double err1 = sqrt(pow(centroids[1].u - 85.2, 2) + pow(centroids[1].v - 20.3, 2));
        if (err0 < 0.1 && err1 < 0.1) {
            printf("Test 3: Happy path multi-star detection... Passed\n");
        } else {
            failures++;
            printf(
                "Test 3: Happy path multi-star detection... FAILED (large subpixel error: %.4f, "
                "%.4f)\n",
                err0, err1);
        }
    } else {
        failures++;
        printf(
            "Test 3: Happy path multi-star detection... FAILED (detected %d centroids, expected "
            "2)\n",
            count_happy);
    }

    /* 4. Giant/Over-saturated Blob (should be filtered out without fragmenting) */
    /* draw a giant filled square of size 40x40 = 1600 pixels (exceeds queue size 1000) */
    for (y = 60; y < 100; y++) {
        for (x = 60; x < 100; x++) {
            image[y * width + x] = 255;
        }
    }
    int count_giant = starfix_find_centroids(image, width, height, 35, 10, centroids, NULL, 1);
    /* should detect star 0 and star 1, and discard the entire 1600-pixel block */
    if (count_giant == 2) {
        printf("Test 4: Giant saturated blob rejection... Passed\n");
    } else {
        failures++;
        printf("Test 4: Giant saturated blob rejection... FAILED (count = %d, expected 2)\n",
               count_giant);
    }

    /* 5. Top-N Brightest Selection When Count Exceeds Capacity */
    printf("Test 5: Top-N brightest selection... ");
    memset(image, 0, (size_t)width * (size_t)height);
    image[10 * width + 10] = 60;  /* dim star in top row */
    image[20 * width + 20] = 80;  /* medium star in middle row */
    image[90 * width + 90] = 250; /* brightest star in bottom row */
    int top_count = starfix_find_centroids(image, width, height, 35, 2, centroids, NULL, 1);
    if (top_count == 2 && fabs(centroids[0].u - 90.0) < 1e-6 &&
        fabs(centroids[1].u - 20.0) < 1e-6) {
        printf("Passed\n");
    } else {
        failures++;
        printf("FAILED (top_count=%d, u0=%.2f, u1=%.2f)\n", top_count, centroids[0].u,
               centroids[1].u);
    }

    /* reset and wrap must leave the next generation usable */
    printf("Test 6: Generation-zero mask handling... ");
    memset(image, 0, (size_t)width * (size_t)height);
    image[60 * width + 60] = 255;
    uint32_t* visited = (uint32_t*)calloc((size_t)width * (size_t)height, sizeof(uint32_t));
    if (!visited) return 1;
    uint32_t generations[] = {0, 1, 2, UINT32_MAX, 0, 1};
    int generations_ok = 1;
    for (size_t generation = 0; generation < sizeof(generations) / sizeof(generations[0]);
         generation++) {
        memset(image, 0, (size_t)width * (size_t)height);
        image[60 * width + (generation >= 3 ? 61 : 60)] = 255;
        int count = starfix_find_centroids(image, width, height, 35, 10, centroids, visited,
                                           generations[generation]);
        generations_ok &= count == 1;
    }
    if (generations_ok) {
        printf("Passed\n");
    } else {
        failures++;
        printf("FAILED (star lost across reset or wrap)\n");
    }
    free(visited);

    printf("Test 7: Count limited to output capacity... ");
    image[80 * width + 80] = 255;
    starfix_centroid_t limited[2] = {{0}, {.num_pixels = -1}};
    int limited_count = starfix_find_centroids(image, width, height, 35, 1, limited, NULL, 1);
    if (limited_count == 1 && limited[1].num_pixels == -1) {
        printf("Passed\n");
    } else {
        failures++;
        printf("FAILED (count=%d)\n", limited_count);
    }

    printf("Test 8: Complete draining of wrapped and full blob queues... ");
    if (large_component_checks()) {
        printf("Passed\n");
    } else {
        failures++;
        printf("FAILED (large component fragmented or isolated stars lost)\n");
    }

    free(image);
    return failures ? 1 : 0;
}
