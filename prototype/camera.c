#include "camera.h"

#include <fcntl.h>
#include <linux/videodev2.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/ioctl.h>
#include <sys/mman.h>
#include <unistd.h>

static int cam_fd = -1;
static uint8_t* cam_buffer = NULL;
static size_t cam_buf_length = 0;

int camera_init(const char* dev_name, int width, int height) {
    cam_fd = open(dev_name, O_RDWR);
    if (cam_fd < 0) {
        perror("Failed to open camera device");
        return -1;
    }

    struct v4l2_format fmt = {0};
    fmt.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
    fmt.fmt.pix.width = width;
    fmt.fmt.pix.height = height;
    /* Request 8-bit grayscale format */
    fmt.fmt.pix.pixelformat = V4L2_PIX_FMT_GREY; 
    fmt.fmt.pix.field = V4L2_FIELD_NONE;

    if (ioctl(cam_fd, VIDIOC_S_FMT, &fmt) < 0) {
        perror("Setting Pixel Format");
        close(cam_fd);
        cam_fd = -1;
        return -1;
    }

    struct v4l2_requestbuffers req = {0};
    req.count = 1;
    req.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
    req.memory = V4L2_MEMORY_MMAP;

    if (ioctl(cam_fd, VIDIOC_REQBUFS, &req) < 0) {
        perror("Requesting Buffer");
        close(cam_fd);
        cam_fd = -1;
        return -1;
    }

    struct v4l2_buffer buf = {0};
    buf.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
    buf.memory = V4L2_MEMORY_MMAP;
    buf.index = 0;

    if (ioctl(cam_fd, VIDIOC_QUERYBUF, &buf) < 0) {
        perror("Querying Buffer");
        close(cam_fd);
        cam_fd = -1;
        return -1;
    }

    cam_buf_length = buf.length;
    cam_buffer = mmap(NULL, buf.length, PROT_READ | PROT_WRITE, MAP_SHARED, cam_fd, buf.m.offset);
    if (cam_buffer == MAP_FAILED) {
        perror("mmap");
        close(cam_fd);
        cam_fd = -1;
        return -1;
    }

    /* Start streaming */
    int type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
    if (ioctl(cam_fd, VIDIOC_STREAMON, &type) < 0) {
        perror("VIDIOC_STREAMON");
        return -1;
    }

    return 0;
}

int camera_capture_frame(uint8_t** buffer) {
    if (cam_fd < 0 || cam_buffer == NULL) return -1;

    struct v4l2_buffer buf = {0};
    buf.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
    buf.memory = V4L2_MEMORY_MMAP;

    /* Enqueue buffer */
    if (ioctl(cam_fd, VIDIOC_QBUF, &buf) < 0) {
        perror("VIDIOC_QBUF");
        return -1;
    }

    /* Dequeue buffer (blocks until frame is ready) */
    if (ioctl(cam_fd, VIDIOC_DQBUF, &buf) < 0) {
        perror("VIDIOC_DQBUF");
        return -1;
    }

    *buffer = cam_buffer;
    return 0;
}

void camera_close(void) {
    if (cam_fd >= 0) {
        int type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
        ioctl(cam_fd, VIDIOC_STREAMOFF, &type);
        if (cam_buffer != NULL && cam_buffer != MAP_FAILED) {
            munmap(cam_buffer, cam_buf_length);
        }
        close(cam_fd);
        cam_fd = -1;
        cam_buffer = NULL;
    }
}
