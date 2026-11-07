#ifndef PROTOTYPE_BNO055_H
#define PROTOTYPE_BNO055_H

#include <stdint.h>

/* I2C addresses */
#define BNO055_I2C_ADDR1 0x28
#define BNO055_I2C_ADDR2 0x29

/* initialize BNO055 sensor; returns 0 on success, < 0 on error */
int bno055_init(const char* i2c_dev, uint8_t addr);

/* read normalized gravity vector [gx, gy, gz] */
int bno055_read_gravity(double* gx, double* gy, double* gz);

/* read euler heading (yaw) in degrees [0, 360] */
int bno055_read_heading(double* heading_deg);

/* read euler pitch and roll in degrees */
int bno055_read_pitch_roll(double* pitch_deg, double* roll_deg);

/* close BNO055 sensor and release resources */
void bno055_close(void);

#endif /* PROTOTYPE_BNO055_H */
