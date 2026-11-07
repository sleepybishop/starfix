#ifndef PROTOTYPE_SIM_H
#define PROTOTYPE_SIM_H

#include "../src/starfix_identify.h"

typedef struct {
    const starfix_db_t* db;
    double gha_aries;
    double true_lat;
    double true_lon;
    double true_heading_deg;
    double true_speed_knots;
    double tilt_x;
    double tilt_y;
    double wave_accel_x;
    double wave_accel_y;
    double roll_rate_rad_s;
    double pitch_rate_rad_s;
    double yaw_jitter_deg;
    int is_daytime;
    double cloud_cover;
} prototype_sim_state_t;

void prototype_sim_set_state(const prototype_sim_state_t* state);
const prototype_sim_state_t* prototype_sim_get_state(void);

#endif /* PROTOTYPE_SIM_H */
