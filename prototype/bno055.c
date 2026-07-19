#include "bno055.h"

#include <fcntl.h>
#include <linux/i2c-dev.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/ioctl.h>
#include <unistd.h>

static int i2c_fd = -1;

/* Register Addresses */
#define BNO055_OPR_MODE_ADDR 0x3D
#define BNO055_GRAVITY_DATA_X_LSB_ADDR 0x2E
#define BNO055_EULER_H_LSB_ADDR 0x1A

/* Operation Modes */
#define OPERATION_MODE_NDOF 0x0C

static int i2c_read_bytes(uint8_t reg, uint8_t* data, int length) {
    if (i2c_fd < 0) return -1;
    if (write(i2c_fd, &reg, 1) != 1) return -1;
    if (read(i2c_fd, data, length) != length) return -1;
    return 0;
}

static int i2c_write_byte(uint8_t reg, uint8_t data) {
    if (i2c_fd < 0) return -1;
    uint8_t buf[2] = {reg, data};
    if (write(i2c_fd, buf, 2) != 2) return -1;
    return 0;
}

int bno055_init(const char* i2c_dev, uint8_t addr) {
    i2c_fd = open(i2c_dev, O_RDWR);
    if (i2c_fd < 0) {
        perror("Failed to open I2C bus");
        return -1;
    }

    if (ioctl(i2c_fd, I2C_SLAVE, addr) < 0) {
        perror("Failed to acquire bus access and/or talk to slave");
        close(i2c_fd);
        i2c_fd = -1;
        return -1;
    }

    /* Set operation mode to NDOF (Nine Degrees of Freedom Sensor Fusion) */
    if (i2c_write_byte(BNO055_OPR_MODE_ADDR, OPERATION_MODE_NDOF) < 0) {
        return -1;
    }

    /* Wait for mode switch to complete (takes up to 7ms) */
    usleep(20000);

    return 0;
}

int bno055_read_gravity(double* gx, double* gy, double* gz) {
    uint8_t buf[6];
    if (i2c_read_bytes(BNO055_GRAVITY_DATA_X_LSB_ADDR, buf, 6) < 0) {
        return -1;
    }

    int16_t raw_x = (int16_t)((buf[1] << 8) | buf[0]);
    int16_t raw_y = (int16_t)((buf[3] << 8) | buf[2]);
    int16_t raw_z = (int16_t)((buf[5] << 8) | buf[4]);

    /* BNO055 Gravity output is 100 LSB per m/s^2 */
    double x = (double)raw_x / 100.0;
    double y = (double)raw_y / 100.0;
    double z = (double)raw_z / 100.0;

    /* Normalize the gravity vector to give Zenith pointing vector in IMU frame */
    double norm = sqrt(x * x + y * y + z * z);
    if (norm > 1e-6) {
        *gx = x / norm;
        *gy = y / norm;
        *gz = z / norm;
    } else {
        *gx = 0.0;
        *gy = 0.0;
        *gz = 1.0;
    }

    return 0;
}

int bno055_read_heading(double* heading_deg) {
    uint8_t buf[2];
    if (i2c_read_bytes(BNO055_EULER_H_LSB_ADDR, buf, 2) < 0) {
        return -1;
    }

    int16_t raw_h = (int16_t)((buf[1] << 8) | buf[0]);
    
    /* BNO055 Euler output is 16 LSB per degree */
    *heading_deg = (double)raw_h / 16.0;
    return 0;
}

int bno055_read_pitch_roll(double* pitch_deg, double* roll_deg) {
    uint8_t buf[4];
    /* Address 0x1C is Euler Pitch LSB, 0x1E is Euler Roll LSB */
    if (i2c_read_bytes(0x1C, buf, 4) < 0) {
        return -1;
    }

    int16_t raw_p = (int16_t)((buf[1] << 8) | buf[0]);
    int16_t raw_r = (int16_t)((buf[3] << 8) | buf[2]);

    *pitch_deg = (double)raw_p / 16.0;
    *roll_deg = (double)raw_r / 16.0;
    return 0;
}
