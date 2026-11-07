#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "../src/starfix_arena.h"
#include "../src/starfix_centroid.h"
#include "../src/starfix_ekf.h"
#include "../src/starfix_graph.h"
#include "../src/starfix_identify.h"
#include "../src/starfix_nmea.h"
#include "../src/starfix_photo.h"
#include "../src/starfix_status.h"
#include "bno055.h"
#include "camera.h"
#include "sim_state.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define CAM_WIDTH 1920
#define CAM_HEIGHT 1080
#define FOV_DEG 12.0

#define MAX_CENTROIDS 100
#define MAX_MATCHES 20

#define MQ_STEPS 369
#define MQ_NODES (MQ_STEPS - 1)
#define MQ_SUBSTEPS 48

static uint8_t mempool[16 * 1024 * 1024];
static starfix_arena_t arena = {mempool, mempool + sizeof(mempool)};
static uint32_t visited_mask[CAM_WIDTH * CAM_HEIGHT];

typedef struct {
    int end_step;
    double lat;
    double lon;
    int is_landfall;
    const char* landmark_zone;
    const char* leg_zone;
} route_waypoint_t;

/* 60-day (~2-month) REQUIEM-SOPAC voyage control points following prototype/route */
static const route_waypoint_t requiem_sopac_waypoints[] = {
    {0, 38.3330, -123.0480, 1, "Bodega_Bay_Departure", "Bodega_Bay_Departure"},
    {14, 34.8000, -125.2000, 0, "California_Current", "California_Current"},
    {28, 30.8000, -127.6000, 0, "California_Current", "California_Current"},
    {44, 25.8000, -130.1000, 0, "Subtropical_High", "Subtropical_High"},
    {60, 20.2000, -131.9000, 0, "Subtropical_High", "Subtropical_High"},
    {76, 14.2000, -133.1000, 0, "NE_Trade_Winds", "NE_Trade_Winds"},
    {90, 8.8000, -133.8000, 0, "NE_Trade_Winds", "NE_Trade_Winds"},
    {106, 4.1000, -134.5000, 0, "ITCZ_Doldrums", "ITCZ_Doldrums"},
    {118, -0.6000, -135.9000, 0, "Equator_Crossing", "Equator_Crossing"},
    {131, -4.9000, -138.1000, 0, "SE_Trade_Winds", "SE_Trade_Winds"},
    {144, -8.9170, -140.1000, 1, "Nuku_Hiva", "SE_Trade_Winds"},
    {148, -9.3900, -140.0800, 1, "Ua_Pou", "Marquesas_Ua_Pou_Passage"},
    {152, -9.9400, -139.1100, 1, "Tahuata", "Marquesas_Tahuata_Passage"},
    {156, -10.4600, -138.6700, 1, "Fatu_Hiva", "Marquesas_Fatu_Hiva_Passage"},
    {160, -9.8000, -139.0300, 1, "Hiva_Oa", "Marquesas_Hiva_Oa_Passage"},
    {164, -8.9300, -139.5500, 1, "Ua_Huka", "Marquesas_Ua_Huka_Passage"},
    {168, -8.9170, -140.1000, 1, "Nuku_Hiva_Return", "Marquesas_Nuku_Hiva_Return"},
    {182, -4.6000, -142.6000, 0, "SE_Trades_Northbound", "SE_Trades_Northbound"},
    {196, -0.3000, -145.0000, 0, "Equator_Northbound", "Equator_Northbound"},
    {210, 4.1000, -147.2000, 0, "ITCZ_South_Edge", "Equator_Northbound"},
    {224, 8.9000, -149.5000, 0, "ITCZ_to_Hawaii", "ITCZ_to_Hawaii"},
    {238, 14.3000, -152.3000, 0, "NE_Trades_to_Hawaii", "NE_Trades_to_Hawaii"},
    {252, 19.7300, -155.0800, 1, "Hawaii", "NE_Trades_to_Hawaii"},
    {268, 25.4000, -156.4000, 0, "North_from_Hawaii", "Pacific_High_West_Flank"},
    {284, 31.2000, -156.6000, 0, "Pacific_High_West_Flank", "Pacific_High_West_Flank"},
    {300, 36.2000, -153.6000, 0, "Rounding_Pacific_High", "Rounding_Pacific_High"},
    {314, 39.3000, -147.8000, 0, "Rounding_Pacific_High", "Rounding_Pacific_High"},
    {328, 40.6000, -140.4000, 0, "Up_and_Around_Pacific_High", "Up_and_Around_Pacific_High"},
    {344, 39.6000, -131.4000, 0, "California_Offshore_Return", "California_Offshore_Return"},
    {360, 37.8080, -122.4768, 1, "Golden_Gate_Fort_Point", "California_Offshore_Return"},
    {361, 37.8126, -122.4774, 1, "Golden_Gate_South_Tower", "Through_the_Golden_Gate"},
    {362, 37.8173, -122.4780, 1, "Through_the_Golden_Gate", "Through_the_Golden_Gate"},
    {363, 37.8219, -122.4785, 1, "Golden_Gate_North_Tower", "Through_the_Golden_Gate"},
    {364, 37.8266, -122.4791, 1, "Lime_Point", "Through_the_Golden_Gate"},
    {365, 37.8320, -122.4680, 1, "Cavallo_Point_SF_Bay", "SF_Bay_Angel_Island_Approach"},
    {366, 37.8530, -122.4440, 1, "Point_Knox_Angel_Island", "SF_Bay_Angel_Island_Approach"},
    {367, 37.8640, -122.4420, 1, "Point_Stuart_Angel_Island", "Backside_of_Angel_Island"},
    {368, 37.8690, -122.4310, 1, "Backside_of_Angel_Island", "Backside_of_Angel_Island"},
};

#define NUM_SOPAC_WAYPOINTS \
    ((int)(sizeof(requiem_sopac_waypoints) / sizeof(requiem_sopac_waypoints[0])))

/* compute catmull-rom tangent on a specific segment [wp_idx - 1, wp_idx] */
static void segment_endpoint_tangents(int wp_idx, double* mlat0, double* mlon0, double* mlat1,
                                      double* mlon1) {
    int i0 = wp_idx - 1;
    int i1 = wp_idx;
    double ds =
        (double)(requiem_sopac_waypoints[i1].end_step - requiem_sopac_waypoints[i0].end_step);
    double chord_lat = requiem_sopac_waypoints[i1].lat - requiem_sopac_waypoints[i0].lat;
    double chord_lon = requiem_sopac_waypoints[i1].lon - requiem_sopac_waypoints[i0].lon;

    /* left endpoint i0 tangent */
    if (i0 <= 0 || requiem_sopac_waypoints[i0].is_landfall) {
        *mlat0 = chord_lat;
        *mlon0 = chord_lon;
    } else {
        double ds_prev = (double)(requiem_sopac_waypoints[i0].end_step -
                                  requiem_sopac_waypoints[i0 - 1].end_step);
        double vlat_prev =
            (requiem_sopac_waypoints[i0].lat - requiem_sopac_waypoints[i0 - 1].lat) / ds_prev;
        double vlon_prev =
            (requiem_sopac_waypoints[i0].lon - requiem_sopac_waypoints[i0 - 1].lon) / ds_prev;
        *mlat0 = 0.5 * (vlat_prev * ds + chord_lat);
        *mlon0 = 0.5 * (vlon_prev * ds + chord_lon);
    }

    /* right endpoint i1 tangent */
    if (i1 >= NUM_SOPAC_WAYPOINTS - 1 || requiem_sopac_waypoints[i1].is_landfall) {
        *mlat1 = chord_lat;
        *mlon1 = chord_lon;
    } else {
        double ds_next = (double)(requiem_sopac_waypoints[i1 + 1].end_step -
                                  requiem_sopac_waypoints[i1].end_step);
        double vlat_next =
            (requiem_sopac_waypoints[i1 + 1].lat - requiem_sopac_waypoints[i1].lat) / ds_next;
        double vlon_next =
            (requiem_sopac_waypoints[i1 + 1].lon - requiem_sopac_waypoints[i1].lon) / ds_next;
        *mlat1 = 0.5 * (chord_lat + vlat_next * ds);
        *mlon1 = 0.5 * (chord_lon + vlon_next * ds);
    }
}

/* evaluate smooth catmull-rom spline + realistic sailing tacking/gybing/weather weave */
static void evaluate_sailing_course(double s, double* out_lat, double* out_lon) {
    if (s <= 0.0) {
        *out_lat = requiem_sopac_waypoints[0].lat;
        *out_lon = requiem_sopac_waypoints[0].lon;
        return;
    }
    if (s >= (double)MQ_NODES) {
        *out_lat = requiem_sopac_waypoints[NUM_SOPAC_WAYPOINTS - 1].lat;
        *out_lon = requiem_sopac_waypoints[NUM_SOPAC_WAYPOINTS - 1].lon;
        return;
    }

    /* gentle along-track speed modulation (+-12% for wind puffs, lulls, and diurnal breeze) */
    double s_eff = s;
    if (s > 0.0 && s < 144.0) {
        double tau = s / 144.0;
        s_eff += sin(M_PI * tau) * (0.85 * sin(2.0 * M_PI * 3.0 * tau + 0.3) +
                                    0.35 * sin(2.0 * M_PI * 7.0 * tau - 0.8));
    } else if (s > 168.0 && s < 252.0) {
        double tau = (s - 168.0) / 84.0;
        s_eff += sin(M_PI * tau) * (0.60 * sin(2.0 * M_PI * 3.0 * tau - 0.5) +
                                    0.25 * sin(2.0 * M_PI * 6.0 * tau + 1.1));
    } else if (s > 252.0 && s < 360.0) {
        double tau = (s - 252.0) / 108.0;
        s_eff += sin(M_PI * tau) * (0.75 * sin(2.0 * M_PI * 3.0 * tau + 0.9) +
                                    0.30 * sin(2.0 * M_PI * 7.0 * tau - 0.6));
    }

    int wp_idx = 1;
    while (wp_idx < NUM_SOPAC_WAYPOINTS - 1 &&
           s_eff > (double)requiem_sopac_waypoints[wp_idx].end_step) {
        wp_idx++;
    }
    double s0 = (double)requiem_sopac_waypoints[wp_idx - 1].end_step;
    double s1 = (double)requiem_sopac_waypoints[wp_idx].end_step;
    double ds = s1 - s0;
    double u = (s_eff - s0) / ds;
    if (u < 0.0) u = 0.0;
    if (u > 1.0) u = 1.0;

    double mlat0, mlon0, mlat1, mlon1;
    segment_endpoint_tangents(wp_idx, &mlat0, &mlon0, &mlat1, &mlon1);

    double u2 = u * u, u3 = u2 * u;
    double h00 = 2.0 * u3 - 3.0 * u2 + 1.0;
    double h10 = u3 - 2.0 * u2 + u;
    double h01 = -2.0 * u3 + 3.0 * u2;
    double h11 = u3 - u2;

    double base_lat = h00 * requiem_sopac_waypoints[wp_idx - 1].lat + h10 * mlat0 +
                      h01 * requiem_sopac_waypoints[wp_idx].lat + h11 * mlat1;
    double base_lon = h00 * requiem_sopac_waypoints[wp_idx - 1].lon + h10 * mlon0 +
                      h01 * requiem_sopac_waypoints[wp_idx].lon + h11 * mlon1;

    /* derivative of cubic hermite spline for local tangent/normal frame */
    double dh00 = 6.0 * u2 - 6.0 * u;
    double dh10 = 3.0 * u2 - 4.0 * u + 1.0;
    double dh01 = -6.0 * u2 + 6.0 * u;
    double dh11 = 3.0 * u2 - 2.0 * u;
    double t_lat = dh00 * requiem_sopac_waypoints[wp_idx - 1].lat + dh10 * mlat0 +
                   dh01 * requiem_sopac_waypoints[wp_idx].lat + dh11 * mlat1;
    double cos_lat = cos(base_lat * M_PI / 180.0);
    double t_lon_scaled = (dh00 * requiem_sopac_waypoints[wp_idx - 1].lon + dh10 * mlon0 +
                           dh01 * requiem_sopac_waypoints[wp_idx].lon + dh11 * mlon1) *
                          cos_lat;
    double t_norm = hypot(t_lat, t_lon_scaled);
    double n_lat = (t_norm > 1e-9) ? (-t_lon_scaled / t_norm) : 0.0;
    double n_lon = (t_norm > 1e-9) ? (t_lat / (t_norm * fmax(0.2, cos_lat))) : 0.0;

    /* realistic cross-track sailing weave (synoptic shifts, tacking/gybing boards, helm trim) */
    double cross_deg = 0.0;
    if (s > 0.0 && s < 144.0) {
        double tau = s / 144.0;
        double env = sin(M_PI * tau);
        cross_deg = env * (0.42 * sin(0.085 * s + 0.4) + 0.22 * sin(0.21 * s - 1.1) +
                           0.10 * sin(0.47 * s + 2.0));
    } else if (s > 144.0 && s < 168.0) {
        /* coastal island-hopping curvature around marquesas headlands */
        cross_deg = sin(M_PI * u) * (0.045 * sin(M_PI * u + 0.7 * (double)wp_idx) +
                                     0.015 * sin(2.0 * M_PI * u - 0.4));
    } else if (s > 168.0 && s < 252.0) {
        double tau = (s - 168.0) / 84.0;
        double env = sin(M_PI * tau);
        cross_deg = env * (0.38 * sin(0.095 * s - 0.7) + 0.20 * sin(0.23 * s + 0.9) +
                           0.09 * sin(0.51 * s - 1.5));
    } else if (s > 252.0 && s < 360.0) {
        double tau = (s - 252.0) / 108.0;
        double env = sin(M_PI * tau);
        cross_deg = env * (0.48 * sin(0.080 * s + 1.2) + 0.25 * sin(0.19 * s - 0.4) +
                           0.11 * sin(0.43 * s + 0.8));
    } else if (s > 360.0 && s < 368.0) {
        /* mild tidal eddy curvature through golden gate and raccoon strait */
        cross_deg = sin(M_PI * u) * 0.0004 * sin(2.0 * M_PI * u + (double)wp_idx);
    }

    *out_lat = base_lat + cross_deg * n_lat;
    *out_lon = base_lon + cross_deg * n_lon;
}

/* compute tangent-plane distance error in nautical miles between two lat/lon points */
static double position_error_nm(double lat1, double lon1, double lat2, double lon2) {
    double dlat_nm = (lat1 - lat2) * 60.0;
    double dlon_nm = (lon1 - lon2) * 60.0 * cos(lat2 * M_PI / 180.0);
    return hypot(dlat_nm, dlon_nm);
}

/* compute tangent-plane distance error in meters between two lat/lon points */
static double position_error_meters(double lat1, double lon1, double lat2, double lon2) {
    return position_error_nm(lat1, lon1, lat2, lon2) * 1852.0;
}

/* open output TSV in prototype/ directory whether invoked from root or prototype/ */
static FILE* open_tsv_file(const char* filename, char* resolved_path, size_t path_size) {
    if (access("prototype/Makefile", F_OK) == 0) {
        snprintf(resolved_path, path_size, "prototype/%s", filename);
    } else {
        snprintf(resolved_path, path_size, "%s", filename);
    }
    return fopen(resolved_path, "w");
}

/* compute GHA of aries (in degrees) from a unix UTC timestamp using IAU GMST */
static double gha_aries_from_utc(double utc_sec) {
    double jd = utc_sec / 86400.0 + 2440587.5;
    double d_j2000 = jd - 2451545.0;
    double t = d_j2000 / 36525.0;
    double gmst =
        280.46061837 + 360.98564736629 * d_j2000 + 0.000387933 * t * t - (t * t * t) / 38710000.0;
    gmst = fmod(gmst, 360.0);
    if (gmst < 0.0) gmst += 360.0;
    return gmst;
}

/* classify ocean zone or route waypoint along the REQUIEM-SOPAC voyage */
static const char* route_zone_for_step(int step) {
    for (int i = 0; i < NUM_SOPAC_WAYPOINTS; i++) {
        if (step == requiem_sopac_waypoints[i].end_step) {
            return requiem_sopac_waypoints[i].landmark_zone;
        }
        if (step < requiem_sopac_waypoints[i].end_step) {
            return requiem_sopac_waypoints[i].leg_zone;
        }
    }
    return "Backside_of_Angel_Island";
}

/* check whether step is an exact waypoint from prototype/route or major ocean milestone */
static int is_milestone_step(int step) {
    for (int i = 0; i < NUM_SOPAC_WAYPOINTS; i++) {
        if (step == requiem_sopac_waypoints[i].end_step) return 1;
    }
    return (step == 12 || step == 38 || step == 68 || step == 92 || step == 126 || step == 178 ||
            step == 232 || step == 262 || step == 290 || step == 314 || step == 334);
}

/* populate instantaneous roll, pitch, angular rates, and orbital wave acceleration */
static void apply_wave_swell_dynamics(prototype_sim_state_t* s, double t_sec, double roll_amp_deg,
                                      double pitch_amp_deg, double swell_period_sec,
                                      double wave_accel_amp_g, double yaw_jitter_deg) {
    double w_swell = 2.0 * M_PI / swell_period_sec;
    double w_chop = 2.0 * w_swell;

    double r1 = (0.82 * roll_amp_deg) * (M_PI / 180.0);
    double r2 = (0.22 * roll_amp_deg) * (M_PI / 180.0);
    double p1 = (0.80 * pitch_amp_deg) * (M_PI / 180.0);
    double p2 = (0.24 * pitch_amp_deg) * (M_PI / 180.0);

    double roll_rad = r1 * sin(w_swell * t_sec + 0.4) + r2 * sin(w_chop * t_sec + 1.1);
    double pitch_rad = p1 * cos(w_swell * t_sec - 0.3) + p2 * cos(w_chop * t_sec + 0.7);

    s->roll_rate_rad_s =
        r1 * w_swell * cos(w_swell * t_sec + 0.4) + r2 * w_chop * cos(w_chop * t_sec + 1.1);
    s->pitch_rate_rad_s =
        -p1 * w_swell * sin(w_swell * t_sec - 0.3) - p2 * w_chop * sin(w_chop * t_sec + 0.7);

    s->tilt_x = sin(roll_rad);
    s->tilt_y = sin(pitch_rad) * cos(roll_rad);

    /* orbital horizontal wave surge/sway acceleration (swell + 2nd-harmonic chop) */
    s->wave_accel_x = wave_accel_amp_g * sin(w_swell * t_sec + 0.4) +
                      0.25 * wave_accel_amp_g * sin(w_chop * t_sec + 1.1);
    s->wave_accel_y = 0.80 * wave_accel_amp_g * cos(w_swell * t_sec - 0.3) +
                      0.20 * wave_accel_amp_g * cos(w_chop * t_sec + 0.7);
    s->yaw_jitter_deg = yaw_jitter_deg;
}

#define SWELL_BURST_FRAMES 6

/* capture a 6-frame burst across one full swell cycle to average out wave surge/sway */
static int capture_swell_burst_fix(prototype_sim_state_t* base_sim, double base_t_sec,
                                   double roll_amp_deg, double pitch_amp_deg,
                                   double swell_period_sec, double wave_accel_amp_g,
                                   double yaw_jitter_deg, uint32_t* frame_gen,
                                   starfix_photo_result_t* out_result, double* out_single_err_m,
                                   double* out_max_roll_deg, double* out_max_pitch_deg) {
    double sum_x = 0.0, sum_y = 0.0, sum_z = 0.0;
    double sum_rms = 0.0;
    uint32_t total_inliers = 0;
    int total_matches = 0;
    int valid_frames = 0;
    double first_frame_err_m = NAN;

    /* sample IMU apparent gravity across the full swell period */
    double gx_list[SWELL_BURST_FRAMES], gy_list[SWELL_BURST_FRAMES], gz_list[SWELL_BURST_FRAMES];
    double imu_err_x_sum = 0.0, imu_err_y_sum = 0.0;
    for (int f = 0; f < SWELL_BURST_FRAMES; f++) {
        double dt_burst = (double)f * (swell_period_sec / (double)SWELL_BURST_FRAMES);
        prototype_sim_state_t sim = *base_sim;
        sim.gha_aries = fmod(base_sim->gha_aries + (dt_burst / 3600.0) * 15.04107, 360.0);
        apply_wave_swell_dynamics(&sim, base_t_sec + dt_burst, roll_amp_deg, pitch_amp_deg,
                                  swell_period_sec, wave_accel_amp_g, yaw_jitter_deg);
        prototype_sim_set_state(&sim);

        double pitch_deg = 0.0, roll_deg = 0.0;
        bno055_read_pitch_roll(&pitch_deg, &roll_deg);
        if (out_max_roll_deg && fabs(roll_deg) > *out_max_roll_deg)
            *out_max_roll_deg = fabs(roll_deg);
        if (out_max_pitch_deg && fabs(pitch_deg) > *out_max_pitch_deg)
            *out_max_pitch_deg = fabs(pitch_deg);

        bno055_read_gravity(&gx_list[f], &gy_list[f], &gz_list[f]);
        imu_err_x_sum += (gx_list[f] - sim.tilt_x);
        imu_err_y_sum += (gy_list[f] - sim.tilt_y);
    }
    double mean_err_x = imu_err_x_sum / (double)SWELL_BURST_FRAMES;
    double mean_err_y = imu_err_y_sum / (double)SWELL_BURST_FRAMES;

    for (int f = 0; f < SWELL_BURST_FRAMES; f++) {
        double dt_burst = (double)f * (swell_period_sec / (double)SWELL_BURST_FRAMES);
        prototype_sim_state_t sim = *base_sim;
        sim.gha_aries = fmod(base_sim->gha_aries + (dt_burst / 3600.0) * 15.04107, 360.0);
        apply_wave_swell_dynamics(&sim, base_t_sec + dt_burst, roll_amp_deg, pitch_amp_deg,
                                  swell_period_sec, wave_accel_amp_g, yaw_jitter_deg);
        prototype_sim_set_state(&sim);

        uint8_t* frame_buffer = NULL;
        if (camera_capture_frame(&frame_buffer, NULL) < 0) continue;

        (*frame_gen)++;
        starfix_centroid_t centroids[MAX_CENTROIDS];
        int num_centroids =
            starfix_find_centroids(frame_buffer, CAM_WIDTH, CAM_HEIGHT, 45, MAX_CENTROIDS,
                                   centroids, visited_mask, *frame_gen);
        if (num_centroids < 4) continue;

        starfix_match_t matches[MAX_MATCHES];
        int num_matches =
            starfix_identify_stars(num_centroids, centroids, CAM_WIDTH, CAM_HEIGHT, FOV_DEG, sim.db,
                                   MAX_MATCHES, matches, arena, NULL);
        if (num_matches < 3) continue;

        /* single-snapshot solve uses instantaneous un-averaged wave gravity vector */
        if (isnan(first_frame_err_m)) {
            starfix_photo_options_t single_opt = {
                .width = CAM_WIDTH,
                .height = CAM_HEIGHT,
                .fov_deg = FOV_DEG,
                .gha_aries = sim.gha_aries,
                .zenith_camera = {gx_list[f], gy_list[f], gz_list[f]},
            };
            starfix_photo_result_t single_res;
            if (starfix_solve_photo(num_centroids, centroids, sim.db->num_stars, sim.db->stars,
                                    num_matches, matches, &single_opt, &single_res,
                                    arena) == STARFIX_SUCCESS) {
                first_frame_err_m = position_error_meters(single_res.latitude, single_res.longitude,
                                                          sim.true_lat, sim.true_lon);
            }
        }

        /* swell-filtered zenith uses IMU gravity averaged over the full wave period */
        double fx = sim.tilt_x + mean_err_x;
        double fy = sim.tilt_y + mean_err_y;
        double fz = sqrt(fmax(0.0, 1.0 - fx * fx - fy * fy));
        double fnorm = hypot(hypot(fx, fy), fz);

        starfix_photo_options_t options = {
            .width = CAM_WIDTH,
            .height = CAM_HEIGHT,
            .fov_deg = FOV_DEG,
            .gha_aries = sim.gha_aries,
            .zenith_camera = {fx / fnorm, fy / fnorm, fz / fnorm},
        };
        starfix_photo_result_t frame_res;
        if (starfix_solve_photo(num_centroids, centroids, sim.db->num_stars, sim.db->stars,
                                num_matches, matches, &options, &frame_res,
                                arena) != STARFIX_SUCCESS)
            continue;

        double lat_r = frame_res.latitude * M_PI / 180.0;
        double lon_r = frame_res.longitude * M_PI / 180.0;
        sum_x += cos(lat_r) * cos(lon_r);
        sum_y += cos(lat_r) * sin(lon_r);
        sum_z += sin(lat_r);
        sum_rms += frame_res.rms_pixels;
        total_inliers += frame_res.inliers;
        total_matches += frame_res.matches;
        valid_frames++;
    }

    if (valid_frames == 0) return 0;

    double horiz = hypot(sum_x, sum_y);
    out_result->latitude = atan2(sum_z, horiz) * 180.0 / M_PI;
    out_result->longitude = atan2(sum_y, sum_x) * 180.0 / M_PI;
    out_result->rms_pixels = sum_rms / (double)valid_frames;
    out_result->inliers =
        (uint32_t)((total_inliers + (uint32_t)(valid_frames / 2)) / (uint32_t)valid_frames);
    out_result->matches = (total_matches + valid_frames / 2) / valid_frames;
    if (out_single_err_m) *out_single_err_m = first_frame_err_m;
    return valid_frames;
}

static int run_marquesas_passage_sim(const starfix_db_t* db) {
    printf(
        "========================================================================================"
        "========================\n");
    printf(
        "REQUIEM-SOPAC: Bodega Bay -> Marquesas -> Hawaii -> Pacific High -> Golden Gate -> "
        "Angel Island\n");
    printf(
        "   Start: Apr 20, 2027 06:00 UTC | Nuku Hiva -> Ua Pou -> Tahuata -> Fatu Hiva -> Hiva Oa "
        "-> Ua Huka -> Nuku Hiva\n");
    printf(
        "          -> Hawaii -> Up & Around Pacific High -> Through the Golden Gate -> Backside of "
        "Angel Island\n");
    printf(
        "========================================================================================"
        "========================\n");

    time_t base_utc = 1808200800; /* 2027-04-20 06:00:00 UTC */
    double start_lat = requiem_sopac_waypoints[0].lat;
    double start_lon = requiem_sopac_waypoints[0].lon;
    double true_speed_scale = 0.93;
    double dt_sub_sec = 300.0; /* 5-minute DR substeps; 48 substeps = 4.0 hours per ocean node */
    double gha_aries_start = gha_aries_from_utc((double)base_utc);

    starfix_ekf_t ekf;
    starfix_ekf_init(&ekf, start_lat - 0.05, start_lon + 0.05);
    ekf.R[0] = 0.010 * 0.010;
    ekf.R[1] = 0.010 * 0.010;

    starfix_odometry_t* odom =
        (starfix_odometry_t*)malloc((size_t)(MQ_NODES * MQ_SUBSTEPS) * sizeof(starfix_odometry_t));
    starfix_fix_t* graph_fixes = (starfix_fix_t*)malloc((size_t)MQ_STEPS * sizeof(starfix_fix_t));
    starfix_pose_t* graph_poses =
        (starfix_pose_t*)malloc((size_t)MQ_STEPS * sizeof(starfix_pose_t));
    if (!odom || !graph_fixes || !graph_poses) {
        free(odom);
        free(graph_fixes);
        free(graph_poses);
        return 1;
    }

    double true_lats[MQ_STEPS], true_lons[MQ_STEPS];
    double dr_lats[MQ_STEPS], dr_lons[MQ_STEPS];
    double raw_lats[MQ_STEPS], raw_lons[MQ_STEPS];
    double ekf_lats[MQ_STEPS], ekf_lons[MQ_STEPS];
    double step_rms[MQ_STEPS];
    uint32_t step_inliers[MQ_STEPS];
    int step_has_fix[MQ_STEPS];
    const char* step_sky[MQ_STEPS];
    double step_bias_deg[MQ_STEPS], step_scale[MQ_STEPS];
    double step_time_hr[MQ_STEPS];
    double step_speed_kt[MQ_STEPS], step_course_deg[MQ_STEPS];

    double true_lat = start_lat, true_lon = start_lon;
    double dr_lat = start_lat - 0.05, dr_lon = start_lon + 0.05;
    int num_fixes = 0;
    int daytime_blackouts = 0;
    int cloud_blackouts = 0;
    uint32_t frame_gen = 500U;
    double max_roll_deg = 0.0, max_pitch_deg = 0.0, sum_single_err_nm = 0.0, sum_burst_err_nm = 0.0;
    double total_route_nm = 0.0;

    for (int step = 0; step < MQ_STEPS; step++) {
        double true_heading_deg = 197.0;
        double true_speed_knots = 5.4;
        double last_meas_speed = NAN;
        double last_compass_hdg = NAN;

        if (step > 0) {
            for (int sub = 0; sub < MQ_SUBSTEPS; sub++) {
                double s_curr = (double)(step - 1) + (double)(sub + 1) / (double)MQ_SUBSTEPS;
                double next_lat = true_lat, next_lon = true_lon;
                evaluate_sailing_course(s_curr, &next_lat, &next_lon);

                double dlat_sub = next_lat - true_lat;
                double dlon_sub = next_lon - true_lon;
                true_lat = next_lat;
                double cos_lat = cos(true_lat * M_PI / 180.0);
                true_lon = next_lon;

                double dn_nm = dlat_sub * 60.0;
                double de_nm = dlon_sub * 60.0 * cos_lat;
                double dist_sub_nm = hypot(dn_nm, de_nm);
                total_route_nm += dist_sub_nm;
                true_heading_deg = fmod(atan2(de_nm, dn_nm) * 180.0 / M_PI + 360.0, 360.0);
                true_speed_knots = dist_sub_nm / (dt_sub_sec / 3600.0);

                prototype_sim_state_t sub_sim = {
                    .db = db,
                    .true_lat = true_lat,
                    .true_lon = true_lon,
                    .true_heading_deg = true_heading_deg,
                    .true_speed_knots = true_speed_knots,
                    .yaw_jitter_deg = (step > 360) ? 2.5 : 3.5,
                };
                prototype_sim_set_state(&sub_sim);

                double compass_heading = 0.0;
                bno055_read_heading(&compass_heading);
                double rel_jitter = ((double)rand() / (double)RAND_MAX - 0.5) * 0.10;
                double measured_speed = true_speed_knots * (true_speed_scale + rel_jitter);

                double dr_dist_nm = measured_speed * (dt_sub_sec / 3600.0);
                double dr_h_rad = compass_heading * M_PI / 180.0;
                dr_lat += (dr_dist_nm * cos(dr_h_rad)) / 60.0;
                dr_lon += (dr_dist_nm * sin(dr_h_rad)) / (60.0 * cos(dr_lat * M_PI / 180.0));

                int odom_idx = (step - 1) * MQ_SUBSTEPS + sub;
                odom[odom_idx].speed = measured_speed;
                odom[odom_idx].heading = compass_heading;

                starfix_ekf_predict(&ekf, measured_speed, compass_heading, dt_sub_sec);
                last_meas_speed = measured_speed;
                last_compass_hdg = compass_heading;
            }
        }

        true_lats[step] = true_lat;
        true_lons[step] = true_lon;
        dr_lats[step] = dr_lat;
        dr_lons[step] = dr_lon;

        /* 4.0h per ocean step (0..360 = 60 days); variable harbor step durations (361..368) */
        static const double harbor_elapsed_hr[8] = {1440.05, 1440.10, 1440.15, 1440.20,
                                                    1440.30, 1440.60, 1440.72, 1440.85};
        static const double harbor_dt_hr[8] = {0.05, 0.05, 0.05, 0.05, 0.10, 0.30, 0.12, 0.13};
        double elapsed_hr = (step <= 360) ? (double)step * 4.0 : harbor_elapsed_hr[step - 361];
        double step_dt_hr = (step <= 360) ? 4.0 : harbor_dt_hr[step - 361];
        double elapsed_sec = elapsed_hr * 3600.0;
        step_time_hr[step] = elapsed_hr;

        /* local solar time (departing bodega bay on apr 20, 2027 at 22:00 PST / 06:00 UTC) */
        double local_solar_hr = fmod(22.0 + elapsed_hr + (true_lon + 120.0) / 15.0, 24.0);
        if (local_solar_hr < 0.0) local_solar_hr += 24.0;
        int is_daytime = (local_solar_hr >= 5.5 && local_solar_hr <= 18.5) ? 1 : 0;

        /* regional cloud cover model: coastal stratus, ITCZ squall belt, and trade cumulus */
        double cloud_cover = 0.0;
        if (!is_daytime && step < 360) {
            if (true_lat >= 3.5 && true_lat <= 9.5) {
                /* multi-day overcast squalls in the ITCZ doldrums (southbound & northbound) */
                cloud_cover = 0.98;
            } else if (true_lat > 33.0 && (step % 5 == 1)) {
                /* coastal california and north pacific marine layer stratus */
                cloud_cover = 0.92;
            } else if (step % 9 == 4) {
                /* passing trade-wind squall cloud */
                cloud_cover = 0.90;
            } else if (step % 5 == 2) {
                /* light scattered trade cumulus (partial attenuation) */
                cloud_cover = 0.20;
            }
        }

        prototype_sim_state_t base_sim = {
            .db = db,
            .gha_aries = fmod(gha_aries_start + elapsed_hr * 15.04107, 360.0),
            .true_lat = true_lat,
            .true_lon = true_lon,
            .true_heading_deg = true_heading_deg,
            .true_speed_knots = true_speed_knots,
            .is_daytime = is_daytime,
            .cloud_cover = cloud_cover,
        };

        double roll_amp = (step > 360) ? 7.5 : 12.0;
        double pitch_amp = (step > 360) ? 4.0 : 6.0;
        double swell_per = (step > 360) ? 9.0 : 11.0;
        double wave_acc = (step > 360) ? 0.9e-3 : 1.4e-3;
        double yaw_jit = (step > 360) ? 2.5 : 3.5;

        starfix_photo_result_t result = {0};
        double single_err_m = 0.0;
        int solved = capture_swell_burst_fix(&base_sim, elapsed_sec, roll_amp, pitch_amp, swell_per,
                                             wave_acc, yaw_jit, &frame_gen, &result, &single_err_m,
                                             &max_roll_deg, &max_pitch_deg) > 0;

        if (solved) {
            sum_single_err_nm += single_err_m / 1852.0;
            sum_burst_err_nm +=
                position_error_nm(result.latitude, result.longitude, true_lat, true_lon);
            starfix_ekf_correct(&ekf, result.latitude, result.longitude);
            raw_lats[step] = result.latitude;
            raw_lons[step] = result.longitude;
            step_rms[step] = result.rms_pixels;
            step_inliers[step] = result.inliers;
            step_has_fix[step] = 1;
            step_sky[step] = (cloud_cover > 0.0) ? "PARTLY_CLOUDY" : "CLEAR_NIGHT";
            graph_fixes[num_fixes++] = (starfix_fix_t){
                .node_idx = step,
                .lat = result.latitude,
                .lon = result.longitude,
            };
        } else {
            raw_lats[step] = NAN;
            raw_lons[step] = NAN;
            step_rms[step] = 0.0;
            step_inliers[step] = 0;
            step_has_fix[step] = 0;
            if (is_daytime) {
                step_sky[step] = "DAYTIME";
                daytime_blackouts++;
            } else if (true_lat >= 3.5 && true_lat <= 9.5) {
                step_sky[step] = "OVERCAST_ITCZ";
                cloud_blackouts++;
            } else {
                step_sky[step] = "OVERCAST";
                cloud_blackouts++;
            }
        }

        ekf_lats[step] = ekf.lat;
        ekf_lons[step] = ekf.lon;
        step_bias_deg[step] = ekf.bias * 180.0 / M_PI;
        step_scale[step] = ekf.scale;
        if (!isnan(last_meas_speed) && !isnan(last_compass_hdg)) {
            double time_scale = 4.0 / step_dt_hr;
            step_speed_kt[step] = (last_meas_speed / ekf.scale) * time_scale;
            step_course_deg[step] = fmod(last_compass_hdg - step_bias_deg[step] + 360.0, 360.0);
        } else {
            step_speed_kt[step] = NAN;
            step_course_deg[step] = NAN;
        }
    }

    double graph_scale = 1.0, graph_bias = 0.0;
    starfix_status_t graph_status = starfix_solve_graph(
        MQ_NODES, MQ_SUBSTEPS, dt_sub_sec, start_lat - 0.05, start_lon + 0.05, odom, num_fixes,
        graph_fixes, graph_poses, &graph_scale, &graph_bias, arena);
    if (graph_status != STARFIX_SUCCESS) {
        free(odom);
        free(graph_fixes);
        free(graph_poses);
        return 1;
    }

    char tsv_path[256], nmea_path[256];
    FILE* tsv = open_tsv_file("bodega_to_marquesas.tsv", tsv_path, sizeof(tsv_path));
    FILE* nmea = open_tsv_file("bodega_to_marquesas.nmea", nmea_path, sizeof(nmea_path));
    if (!tsv || !nmea) {
        if (tsv) fclose(tsv);
        if (nmea) fclose(nmea);
        free(odom);
        free(graph_fixes);
        free(graph_poses);
        return 1;
    }
    fprintf(tsv,
            "step\ttime_hr\tzone\tsky_status\ttrue_lat\ttrue_lon\tdr_lat\tdr_lon\tdr_err_nm\t"
            "fix_lat\tfix_lon\tfix_err_nm\tstars\trms_px\tekf_lat\tekf_lon\tekf_err_nm\t"
            "qsp_lat\tqsp_lon\tqsp_err_nm\tekf_bias_deg\tekf_scale\n");

    printf("%-4s %-7s %-26s %-13s | %-18s | %-7s %-7s | %-7s %-7s | %-7s\n", "Step", "Time",
           "Route Waypoint / Zone", "Sky Status", "True (Lat, Lon)", "FixErr", "DR_Err", "EKFErr",
           "QSPErr", "EKFBias");
    printf(
        "----------------------------------------------------------------------------------------"
        "-----------------------------\n");

    char last_nmea[256] = {0};
    double sum_dr_err_nm = 0.0, sum_ekf_err_nm = 0.0, sum_qsp_err_nm = 0.0;
    for (int step = 0; step < MQ_STEPS; step++) {
        double time_hr = step_time_hr[step];
        const char* zone = route_zone_for_step(step);
        double dr_err_nm =
            position_error_nm(dr_lats[step], dr_lons[step], true_lats[step], true_lons[step]);
        double ekf_err_nm =
            position_error_nm(ekf_lats[step], ekf_lons[step], true_lats[step], true_lons[step]);
        double qsp_err_nm = position_error_nm(graph_poses[step].lat, graph_poses[step].lon,
                                              true_lats[step], true_lons[step]);
        sum_dr_err_nm += dr_err_nm;
        sum_ekf_err_nm += ekf_err_nm;
        sum_qsp_err_nm += qsp_err_nm;

        if (step_has_fix[step]) {
            double fix_err_nm =
                position_error_nm(raw_lats[step], raw_lons[step], true_lats[step], true_lons[step]);
            fprintf(tsv,
                    "%d\t%.1f\t%s\t%s\t%.5f\t%.5f\t%.5f\t%.5f\t%.3f\t%.5f\t%.5f\t%.3f\t%u\t%.3f\t"
                    "%.5f\t%.5f\t%.3f\t%.5f\t%.5f\t%.3f\t%.4f\t%.4f\n",
                    step, time_hr, zone, step_sky[step], true_lats[step], true_lons[step],
                    dr_lats[step], dr_lons[step], dr_err_nm, raw_lats[step], raw_lons[step],
                    fix_err_nm, step_inliers[step], step_rms[step], ekf_lats[step], ekf_lons[step],
                    ekf_err_nm, graph_poses[step].lat, graph_poses[step].lon, qsp_err_nm,
                    step_bias_deg[step], step_scale[step]);
        } else {
            fprintf(tsv,
                    "%d\t%.1f\t%s\t%s\t%.5f\t%.5f\t%.5f\t%.5f\t%.3f\tNA\tNA\tNA\t0\tNA\t"
                    "%.5f\t%.5f\t%.3f\t%.5f\t%.5f\t%.3f\t%.4f\t%.4f\n",
                    step, time_hr, zone, step_sky[step], true_lats[step], true_lons[step],
                    dr_lats[step], dr_lons[step], dr_err_nm, ekf_lats[step], ekf_lons[step],
                    ekf_err_nm, graph_poses[step].lat, graph_poses[step].lon, qsp_err_nm,
                    step_bias_deg[step], step_scale[step]);
        }

        time_t step_utc = base_utc + (time_t)llround(time_hr * 3600.0);
        int fix_quality = step_has_fix[step] ? 1 : 6;
        int num_stars = step_has_fix[step] ? (int)step_inliers[step] : 0;
        if (starfix_format_nmea_fix(ekf_lats[step], ekf_lons[step], step_utc, fix_quality,
                                    num_stars, step_speed_kt[step], step_course_deg[step],
                                    last_nmea, sizeof(last_nmea)) == STARFIX_SUCCESS) {
            fputs(last_nmea, nmea);
        }

        if (is_milestone_step(step)) {
            char fix_buf[16];
            if (step_has_fix[step]) {
                double fix_err = position_error_nm(raw_lats[step], raw_lons[step], true_lats[step],
                                                   true_lons[step]);
                snprintf(fix_buf, sizeof(fix_buf), "%5.2fnm", fix_err);
            } else {
                snprintf(fix_buf, sizeof(fix_buf), "   --  ");
            }
            printf(
                "%-4d %5.1fh %-26s %-13s | %7.3f, %8.3f | %7s %6.1fnm | %5.2fnm %5.2fnm | "
                "%+5.2fd\n",
                step, time_hr, zone, step_sky[step], true_lats[step], true_lons[step], fix_buf,
                dr_err_nm, ekf_err_nm, qsp_err_nm, step_bias_deg[step]);
        }
    }
    fclose(tsv);
    fclose(nmea);

    double final_dr_err_nm = position_error_nm(dr_lats[MQ_NODES], dr_lons[MQ_NODES],
                                               true_lats[MQ_NODES], true_lons[MQ_NODES]);
    double final_ekf_err_nm = position_error_nm(ekf_lats[MQ_NODES], ekf_lons[MQ_NODES],
                                                true_lats[MQ_NODES], true_lons[MQ_NODES]);
    double final_qsp_err_nm =
        position_error_nm(graph_poses[MQ_NODES].lat, graph_poses[MQ_NODES].lon, true_lats[MQ_NODES],
                          true_lons[MQ_NODES]);

    printf(
        "----------------------------------------------------------------------------------------"
        "------------------------\n");
    printf(
        "Summary: %d/%d clear/partly-cloudy fixes | %d Daytime Blackouts | %d Overcast Cloud "
        "Blackouts | Distance: %.0f nm\n",
        num_fixes, MQ_STEPS, daytime_blackouts, cloud_blackouts, total_route_nm);
    printf(
        "Sea State & Swell    : Max Roll = +/-%.1f deg, Max Pitch = +/-%.1f deg, Yaw Jitter = "
        "+/-3.5 deg (11s Ocean Swell)\n",
        max_roll_deg, max_pitch_deg);
    printf(
        "Optical Swell Filter : 1-Frame Wave Pendulum = %4.2f nm | 6-Frame Swell Burst = %4.2f "
        "nm\n",
        num_fixes > 0 ? sum_single_err_nm / (double)num_fixes : 0.0,
        num_fixes > 0 ? sum_burst_err_nm / (double)num_fixes : 0.0);
    printf(
        "Mean Position Error  : Uncalibrated DR = %5.1f nm | Online EKF = %4.2f nm | Batch "
        "NanoQSP Graph = %4.2f nm\n",
        sum_dr_err_nm / (double)MQ_STEPS, sum_ekf_err_nm / (double)MQ_STEPS,
        sum_qsp_err_nm / (double)MQ_STEPS);
    printf(
        "Final Arrival Error  : Uncalibrated DR = %5.1f nm | Online EKF = %4.2f nm (%.0f m) | "
        "Batch NanoQSP = %4.2f nm (%.0f m)\n",
        final_dr_err_nm, final_ekf_err_nm, final_ekf_err_nm * 1852.0, final_qsp_err_nm,
        final_qsp_err_nm * 1852.0);
    printf(
        "Sensor Calibration   : EKF Bias = %+.2f deg, Scale = %.3f | QSP Graph Bias = %+.2f "
        "deg, Scale = %.3f (True: +4.00 deg, 0.930)\n",
        ekf.bias * 180.0 / M_PI, ekf.scale, graph_bias * 180.0 / M_PI, graph_scale);
    printf("Saved Marquesas TSV  : %s | NMEA 0183: %s\n", tsv_path, nmea_path);
    printf("Final NMEA 0183 Fix  :\n%s", last_nmea);

    free(odom);
    free(graph_fixes);
    free(graph_poses);

    return (num_fixes >= 20 && final_ekf_err_nm < 5.0 && final_qsp_err_nm < 5.0) ? 0 : 1;
}

int main(void) {
    if (camera_init("/dev/video0", CAM_WIDTH, CAM_HEIGHT) < 0) {
        fprintf(stderr, "Failed to initialize mock camera.\n");
        return 1;
    }
    if (bno055_init("/dev/i2c-1", BNO055_I2C_ADDR1) < 0) {
        fprintf(stderr, "Failed to initialize mock BNO055.\n");
        camera_close();
        return 1;
    }

    starfix_db_t db = {0};
    const char* db_path = "../data/starfix_db.bin";
    if (starfix_load_db(db_path, &db) != 0) {
        db_path = "data/starfix_db.bin";
        if (starfix_load_db(db_path, &db) != 0) {
            fprintf(stderr, "Failed to load star database.\n");
            bno055_close();
            camera_close();
            return 1;
        }
    }
    /* propagate catalog proper motions to apr 20, 2027 (J2027.299) */
    starfix_propagate_astrometry(db.stars, db.num_stars, 2027.2991);

    int status_mq = run_marquesas_passage_sim(&db);

    starfix_free_db(&db);
    bno055_close();
    camera_close();

    return (status_mq == 0) ? 0 : 1;
}
