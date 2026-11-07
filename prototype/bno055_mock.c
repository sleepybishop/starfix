#include <math.h>
#include <stdlib.h>

#include "bno055.h"
#include "sim_state.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static prototype_sim_state_t current_sim_state;

void prototype_sim_set_state(const prototype_sim_state_t* state) {
    if (state) current_sim_state = *state;
}

const prototype_sim_state_t* prototype_sim_get_state(void) { return &current_sim_state; }

int bno055_init(const char* i2c_dev, uint8_t addr) {
    (void)i2c_dev;
    (void)addr;
    srand(42);
    return 0;
}

static double rand_symmetric(void) {
    /* sum of 3 uniforms for zero-mean bell-shaped jitter in [-1, 1] */
    double u1 = (double)rand() / (double)RAND_MAX - 0.5;
    double u2 = (double)rand() / (double)RAND_MAX - 0.5;
    double u3 = (double)rand() / (double)RAND_MAX - 0.5;
    return (u1 + u2 + u3) * (2.0 / 3.0);
}

int bno055_read_gravity(double* gx, double* gy, double* gz) {
    if (!gx || !gy || !gz) return -1;
    double tx = current_sim_state.tilt_x;
    double ty = current_sim_state.tilt_y;
    double tz = sqrt(fmax(0.0, 1.0 - tx * tx - ty * ty));

    /* residual orbital wave surge/sway acceleration + IMU-camera sync lag + MEMS vibration */
    double dt_sync = rand_symmetric() * 0.0020;
    double nx = current_sim_state.wave_accel_x + current_sim_state.roll_rate_rad_s * dt_sync +
                rand_symmetric() * 2.5e-4;
    double ny = current_sim_state.wave_accel_y + current_sim_state.pitch_rate_rad_s * dt_sync +
                rand_symmetric() * 2.5e-4;
    double nz = rand_symmetric() * 1.2e-4;

    /* apply perturbation in the tangent plane orthogonal to true zenith (tx, ty, tz) */
    double x = tx + nx * tz;
    double y = ty + ny * tz;
    double z = tz - (tx * nx + ty * ny) + nz;
    double norm = hypot(hypot(x, y), z);
    if (norm <= 1e-12) return -1;

    *gx = x / norm;
    *gy = y / norm;
    *gz = z / norm;
    return 0;
}

int bno055_read_heading(double* heading_deg) {
    if (!heading_deg) return -1;
    /* simulated +4.0 deg magnetic compass bias plus sea-state wave yaw jitter */
    double yaw_amp =
        current_sim_state.yaw_jitter_deg > 0.0 ? current_sim_state.yaw_jitter_deg : 2.5;
    double noise = rand_symmetric() * yaw_amp;
    double h = fmod(current_sim_state.true_heading_deg + 4.0 + noise, 360.0);
    if (h < 0.0) h += 360.0;
    *heading_deg = h;
    return 0;
}

int bno055_read_pitch_roll(double* pitch_deg, double* roll_deg) {
    if (!pitch_deg || !roll_deg) return -1;
    *roll_deg = asin(current_sim_state.tilt_x) * 180.0 / M_PI;
    *pitch_deg = asin(current_sim_state.tilt_y) * 180.0 / M_PI;
    return 0;
}

void bno055_close(void) {}
