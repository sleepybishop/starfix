#ifndef PROTOTYPE_CAMERA_H
#define PROTOTYPE_CAMERA_H

#include <stdint.h>

/* Initializes the V4L2 camera (e.g. /dev/video0)
 * Returns 0 on success, < 0 on error. */
int camera_init(const char* dev_name, int width, int height);

/* Captures a single 8-bit grayscale raw frame.
 * The buffer is managed internally and should not be freed by the caller.
 * Returns 0 on success, < 0 on error. */
int camera_capture_frame(uint8_t** buffer);

/* Closes the camera and frees resources */
void camera_close(void);

#endif /* PROTOTYPE_CAMERA_H */
