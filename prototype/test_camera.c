#define _POSIX_C_SOURCE 200809L

#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <sys/ioctl.h>
#include <time.h>

static int test_ioctl(int fd, unsigned long request, ...);
static int test_clock_gettime(clockid_t clock, struct timespec* value);

/* exercise the driver capture path without requiring a V4L2 device */
#define ioctl test_ioctl
#define clock_gettime test_clock_gettime
#include "camera.c"
#undef ioctl
#undef clock_gettime

static uint32_t capture_flags = V4L2_BUF_FLAG_TIMESTAMP_MONOTONIC;
static struct timeval capture_time = {495, 375000};
static int fail_realtime = 0;
static int fail_monotonic = 0;

static int test_ioctl(int fd, unsigned long request, ...) {
    assert(fd == 42);
    va_list args;
    va_start(args, request);
    struct v4l2_buffer* buffer = va_arg(args, struct v4l2_buffer*);
    va_end(args);
    assert(request == VIDIOC_QBUF || request == VIDIOC_DQBUF);
    if (request == VIDIOC_DQBUF) {
        buffer->flags = capture_flags;
        buffer->timestamp = capture_time;
    }
    return 0;
}

static int test_clock_gettime(clockid_t clock, struct timespec* value) {
    if (clock == CLOCK_REALTIME) {
        if (fail_realtime) return -1;
        *value = (struct timespec){1790834406, 500000000};
    } else {
        assert(clock == CLOCK_MONOTONIC);
        if (fail_monotonic) return -1;
        *value = (struct timespec){500, 750000000};
    }
    return 0;
}

int main(void) {
    uint8_t pixel = 0;
    cam_fd = 42;
    cam_buffer = &pixel;
    uint8_t* frame = NULL;
    double utc = 0.0;

    /* driver frame is 5.375 seconds old; preserve its fractional capture time */
    assert(camera_capture_frame(&frame, &utc) == 0);
    assert(frame == &pixel && fabs(utc - 1790834401.125) < 1e-6);

    /* unknown clock domains and missing/invalid metadata use capture completion */
    capture_flags = V4L2_BUF_FLAG_TIMESTAMP_UNKNOWN;
    assert(camera_capture_frame(&frame, &utc) == 0);
    assert(fabs(utc - 1790834406.5) < 1e-6);
    capture_flags = V4L2_BUF_FLAG_TIMESTAMP_COPY;
    assert(camera_capture_frame(&frame, &utc) == 0);
    assert(fabs(utc - 1790834406.5) < 1e-6);
    capture_flags = V4L2_BUF_FLAG_TIMESTAMP_MONOTONIC;
    capture_time = (struct timeval){0, 0};
    assert(camera_capture_frame(&frame, &utc) == 0);
    assert(fabs(utc - 1790834406.5) < 1e-6);
    capture_time = (struct timeval){501, 0};
    assert(camera_capture_frame(&frame, &utc) == 0);
    assert(fabs(utc - 1790834406.5) < 1e-6);
    capture_time = (struct timeval){495, 250000};
    fail_monotonic = 1;
    assert(camera_capture_frame(&frame, &utc) == 0);
    assert(fabs(utc - 1790834406.5) < 1e-6);

    fail_realtime = 1;
    utc = -1.0;
    assert(camera_capture_frame(&frame, &utc) < 0 && utc < 0.0);
    assert(camera_capture_frame(&frame, NULL) == 0);
    puts("Camera timestamp checks passed");
    return 0;
}
