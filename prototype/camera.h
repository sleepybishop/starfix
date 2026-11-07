#ifndef PROTOTYPE_CAMERA_H
#define PROTOTYPE_CAMERA_H

#include <stdint.h>

/* initialize V4L2 camera (e.g. /dev/video0); returns 0 on success, < 0 on error */
int camera_init(const char* dev_name, int width, int height);

/* capture a single 8-bit grayscale frame; buffer is owned by driver.
 * optional captured_utc receives fractional unix seconds from driver timestamp,
 * or capture completion time when driver clock is unknown. */
int camera_capture_frame(uint8_t** buffer, double* captured_utc);

/* close camera and release resources */
void camera_close(void);

#endif /* PROTOTYPE_CAMERA_H */
