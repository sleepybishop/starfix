#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "starfix_centroid.h"
#include "starfix_status.h"

/* test runner to verify the C centroiding implementation under multiple conditions */

int main() {
    starfix_telemetry_t telem = {0};
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
    int ret1 = starfix_find_centroids(NULL, width, height, 35, 10, centroids, NULL, 1, &telem);
    int ret2 = starfix_find_centroids(image, 0, height, 35, 10, centroids, NULL, 1, &telem);
    int ret3 = starfix_find_centroids(image, width, height, 35, 0, centroids, NULL, 1, &telem);

    if (ret1 < 0 && ret2 < 0 && ret3 < 0) {
        printf("Test 1: Boundary check... Passed\n");
    } else {
        printf("Test 1: Boundary check... FAILED (got %d, %d, %d)\n", ret1, ret2, ret3);
    }

    /* 2. Completely Black Image (Thresholding check) */
    int i;
    for (i = 0; i < width * height; i++) {
        image[i] = 10;
    }
    starfix_status_t status_black =
        starfix_find_centroids(image, width, height, 35, 10, centroids, NULL, 1, &telem);
    int count_black = (int)telem.num_stars_detected;
    if (status_black == STARFIX_SUCCESS && count_black == 0) {
        printf("Test 2: Black image (no stars)... Passed\n");
    } else {
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

    starfix_status_t status_happy =
        starfix_find_centroids(image, width, height, 35, 10, centroids, NULL, 1, &telem);
    int count_happy = (int)telem.num_stars_detected;
    if (status_happy == STARFIX_SUCCESS && count_happy == 2) {
        double err0 = sqrt(pow(centroids[0].u - 30.4, 2) + pow(centroids[0].v - 50.8, 2));
        double err1 = sqrt(pow(centroids[1].u - 85.2, 2) + pow(centroids[1].v - 20.3, 2));
        if (err0 < 0.1 && err1 < 0.1) {
            printf("Test 3: Happy path multi-star detection... Passed\n");
        } else {
            printf(
                "Test 3: Happy path multi-star detection... FAILED (large subpixel error: %.4f, "
                "%.4f)\n",
                err0, err1);
        }
    } else {
        printf(
            "Test 3: Happy path multi-star detection... FAILED (detected %d centroids, expected "
            "2)\n",
            count_happy);
    }

    /* 4. Giant/Over-saturated Blob (should be filtered out by size constraint) */
    /* draw a giant filled square of size 25x25 = 625 pixels (limit is 300) */
    for (y = 10; y < 35; y++) {
        for (x = 10; x < 35; x++) {
            image[y * width + x] = 255;
        }
    }
    starfix_status_t status_giant =
        starfix_find_centroids(image, width, height, 35, 10, centroids, NULL, 1, &telem);
    int count_giant = (int)telem.num_stars_detected;
    /* should detect star 0 and star 1, but discard the giant block */
    if (status_giant == STARFIX_SUCCESS && count_giant == 2) {
        printf("Test 4: Giant saturated blob rejection... Passed\n");
    } else {
        printf("Test 4: Giant saturated blob rejection... FAILED (count = %d, expected 2)\n",
               count_giant);
    }

    /* 5. Dynamic Hot Pixel Masking Checks */
    printf("Test 5: Hot pixel calibration and masking... ");
    /* create a mock dark frame (all zeros except 3 hot pixels) */
    unsigned char* dark_frame =
        (unsigned char*)calloc((size_t)width * (size_t)height, sizeof(unsigned char));
    dark_frame[12 * width + 15] = 255; /* hot pixel 1 */
    dark_frame[45 * width + 88] = 240; /* hot pixel 2 */
    dark_frame[77 * width + 33] = 220; /* hot pixel 3 */

    int mask_u[10];
    int mask_v[10];
    int num_hot = starfix_calibrate_hot_pixels(dark_frame, width, height, 100, 10, mask_u, mask_v);

    if (num_hot == 3 && mask_u[0] == 15 && mask_v[0] == 12 && mask_u[1] == 88 && mask_v[1] == 45 &&
        mask_u[2] == 33 && mask_v[2] == 77) {
        /* clear previous image, draw 1 true star and inject same 3 hot pixels */
        memset(image, 0, (size_t)width * (size_t)height * sizeof(unsigned char));
        /* true star at (50, 50) */
        image[50 * width + 50] = 250;
        image[50 * width + 51] = 180;
        image[51 * width + 50] = 180;

        /* inject hot pixels */
        image[12 * width + 15] = 255;
        image[45 * width + 88] = 240;
        image[77 * width + 33] = 220;

        /* apply mask */
        int mask_status =
            starfix_apply_hot_pixel_mask(image, width, height, num_hot, mask_u, mask_v);

        /* detect centroids */
        starfix_status_t status_masked =
            starfix_find_centroids(image, width, height, 35, 10, centroids, NULL, 1, &telem);
        int count_masked = (int)telem.num_stars_detected;

        if (mask_status == 0 && status_masked == STARFIX_SUCCESS && count_masked == 1 &&
            fabs(centroids[0].u - 50.295) < 0.1 && fabs(centroids[0].v - 50.295) < 0.1) {
            printf("Passed\n");
        } else {
            printf("FAILED (status=%d, count_masked=%d, u=%.2f, v=%.2f)\n", mask_status,
                   count_masked, centroids[0].u, centroids[0].v);
        }
    } else {
        printf("FAILED (calibration found %d hot pixels, expected 3)\n", num_hot);
    }

    free(dark_frame);

    free(image);
    return 0;
}
