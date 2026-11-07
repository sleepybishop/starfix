#define _POSIX_C_SOURCE 200809L

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../src/starfix_astrometry.h"
#include "../src/starfix_attitude.h"
#include "../src/starfix_geometry.h"
#include "camera.h"
#include "sim_state.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static uint8_t* mock_buffer = NULL;
static int mock_width = 0;
static int mock_height = 0;

/* build celestial-to-camera rotation matrix so true zenith maps to camera tilt vector */
static void build_camera_attitude(const prototype_sim_state_t* s, double R[3][3],
                                  starfix_vector3_t* g_cam) {
    double lat = s->true_lat * M_PI / 180.0;
    double lst = (s->gha_aries + s->true_lon) * M_PI / 180.0;
    double hdg = s->true_heading_deg * M_PI / 180.0;

    /* local tangent frame on celestial sphere: east, north, zenith */
    starfix_vector3_t east = {-sin(lst), cos(lst), 0.0};
    starfix_vector3_t north = {-sin(lat) * cos(lst), -sin(lat) * sin(lst), cos(lat)};
    starfix_vector3_t zenith_cel = {cos(lat) * cos(lst), cos(lat) * sin(lst), sin(lat)};

    /* rotate horizontal tangent vectors by vessel heading */
    double ch = cos(hdg), sh = sin(hdg);
    starfix_vector3_t body_x = {ch * east.x - sh * north.x, ch * east.y - sh * north.y,
                                ch * east.z - sh * north.z};
    starfix_vector3_t body_y = {sh * east.x + ch * north.x, sh * east.y + ch * north.y,
                                sh * east.z + ch * north.z};

    /* camera zenith unit vector in camera coordinates */
    double gx = s->tilt_x;
    double gy = s->tilt_y;
    double gz = sqrt(fmax(0.0, 1.0 - gx * gx - gy * gy));
    double g_norm = hypot(hypot(gx, gy), gz);
    gx /= g_norm;
    gy /= g_norm;
    gz /= g_norm;
    *g_cam = (starfix_vector3_t){gx, gy, gz};

    /* orthonormal basis in camera frame with column 2 = g_cam */
    double ux = 1.0, uy = 0.0, uz = 0.0;
    double dot = ux * gx + uy * gy + uz * gz;
    ux -= dot * gx;
    uy -= dot * gy;
    uz -= dot * gz;
    double u_norm = hypot(hypot(ux, uy), uz);
    ux /= u_norm;
    uy /= u_norm;
    uz /= u_norm;

    double vx = gy * uz - gz * uy;
    double vy = gz * ux - gx * uz;
    double vz = gx * uy - gy * ux;

    /* R = [u_cam, v_cam, g_cam] * [body_x, body_y, zenith_cel]^T */
    R[0][0] = ux * body_x.x + vx * body_y.x + gx * zenith_cel.x;
    R[0][1] = ux * body_x.y + vx * body_y.y + gx * zenith_cel.y;
    R[0][2] = ux * body_x.z + vx * body_y.z + gx * zenith_cel.z;

    R[1][0] = uy * body_x.x + vy * body_y.x + gy * zenith_cel.x;
    R[1][1] = uy * body_x.y + vy * body_y.y + gy * zenith_cel.y;
    R[1][2] = uy * body_x.z + vy * body_y.z + gy * zenith_cel.z;

    R[2][0] = uz * body_x.x + vz * body_y.x + gz * zenith_cel.x;
    R[2][1] = uz * body_x.y + vz * body_y.y + gz * zenith_cel.y;
    R[2][2] = uz * body_x.z + vz * body_y.z + gz * zenith_cel.z;
}

/* shift true camera ray toward zenith by atmospheric refraction (10 C, 1010 hPa) */
static starfix_vector3_t refract_camera_ray(starfix_vector3_t cam, starfix_vector3_t g_cam) {
    double cos_theta = starfix_dot(cam, g_cam);
    if (cos_theta > 1.0) cos_theta = 1.0;
    if (cos_theta < -1.0) cos_theta = -1.0;

    double theta_t = acos(cos_theta);
    double alt_t = M_PI / 2.0 - theta_t;
    if (theta_t <= 1e-6 || alt_t <= 0.0) return cam;

    /* invert starfix_unrefract via 2-step fixed-point iteration */
    double alt_a = alt_t;
    for (int iter = 0; iter < 2; iter++) {
        double est_t = alt_a;
        starfix_unrefract(1010.0, 10.0, alt_a, &est_t);
        alt_a += (alt_t - est_t);
    }
    double theta_a = M_PI / 2.0 - alt_a;

    double ux = cam.x - cos_theta * g_cam.x;
    double uy = cam.y - cos_theta * g_cam.y;
    double uz = cam.z - cos_theta * g_cam.z;
    double u_norm = hypot(hypot(ux, uy), uz);
    if (u_norm <= 1e-12) return cam;
    ux /= u_norm;
    uy /= u_norm;
    uz /= u_norm;

    starfix_vector3_t refracted = {
        cos(theta_a) * g_cam.x + sin(theta_a) * ux,
        cos(theta_a) * g_cam.y + sin(theta_a) * uy,
        cos(theta_a) * g_cam.z + sin(theta_a) * uz,
    };
    double r_norm = hypot(hypot(refracted.x, refracted.y), refracted.z);
    if (r_norm > 1e-12) {
        refracted.x /= r_norm;
        refracted.y /= r_norm;
        refracted.z /= r_norm;
    }
    return refracted;
}

static uint32_t mock_rng = 0x5A17F123U;

static double mock_rand_uniform(void) {
    mock_rng = mock_rng * 1664525U + 1013904223U;
    return (double)(mock_rng & 0x00FFFFFFU) / (double)0x01000000U;
}

static double mock_rand_jitter(void) { return (mock_rand_uniform() + mock_rand_uniform() - 1.0); }

int camera_init(const char* dev_name, int width, int height) {
    (void)dev_name;
    if (width <= 0 || height <= 0) return -1;
    free(mock_buffer);
    mock_buffer = (uint8_t*)malloc((size_t)width * (size_t)height);
    if (!mock_buffer) return -1;
    mock_width = width;
    mock_height = height;
    mock_rng = 0x5A17F123U;
    return 0;
}

int camera_capture_frame(uint8_t** buffer, double* captured_utc) {
    if (!mock_buffer || !buffer) return -1;
    const prototype_sim_state_t* s = prototype_sim_get_state();
    if (!s || !s->db || !s->db->stars) return -1;
    if (captured_utc) {
        struct timespec now;
        if (clock_gettime(CLOCK_REALTIME, &now) < 0) return -1;
        *captured_utc = (double)now.tv_sec + (double)now.tv_nsec / 1e9;
    }

    /* fill sensor background with dark sky glow + read noise */
    size_t num_pixels = (size_t)mock_width * (size_t)mock_height;
    memset(mock_buffer, 12, num_pixels);

    /* daytime glare or heavy overcast obscures stars from the sensor */
    if (s->is_daytime || s->cloud_cover >= 0.85) {
        *buffer = mock_buffer;
        return 0;
    }

    double R[3][3];
    starfix_vector3_t g_cam;
    build_camera_attitude(s, R, &g_cam);
    double focal = starfix_focal_length(mock_width, 12.0);
    double field_cos = cos(10.0 * M_PI / 180.0);

    for (uint32_t i = 0; i < s->db->num_stars; i++) {
        const starfix_catalog_star_t* star = &s->db->stars[i];
        starfix_vector3_t v = starfix_celestial_ray(star->ra, star->dec);
        starfix_vector3_t cam = starfix_rotate((const double (*)[3])R, v);
        if (cam.z < field_cos) continue;

        /* patchy cloud cover obscures individual stars */
        if (s->cloud_cover > 0.0 && mock_rand_uniform() < s->cloud_cover) continue;

        cam = refract_camera_ray(cam, g_cam);
        if (cam.z <= 0.0) continue;

        /* atmospheric seeing + angular-rate exposure smear centroid jitter */
        double rate_deg_s = hypot(s->roll_rate_rad_s, s->pitch_rate_rad_s) * 180.0 / M_PI;
        double jitter_px = 0.20 + 0.025 * rate_deg_s;
        double u =
            (double)mock_width / 2.0 + focal * (cam.x / cam.z) + mock_rand_jitter() * jitter_px;
        double v_px =
            (double)mock_height / 2.0 + focal * (cam.y / cam.z) + mock_rand_jitter() * jitter_px;
        if (u < 8.0 || u >= (double)mock_width - 8.0 || v_px < 8.0 ||
            v_px >= (double)mock_height - 8.0)
            continue;

        double scint = (1.0 - 0.5 * s->cloud_cover) * (1.0 + mock_rand_jitter() * 0.12);
        double peak = 240.0 * pow(10.0, -0.25 * (star->mag - 1.0)) * scint;
        if (peak < 90.0) peak = 90.0;
        if (peak > 240.0) peak = 240.0;

        int cx = (int)floor(u + 0.5);
        int cy = (int)floor(v_px + 0.5);
        for (int dy = -3; dy <= 3; dy++) {
            for (int dx = -3; dx <= 3; dx++) {
                int px = cx + dx;
                int py = cy + dy;
                if (px < 0 || px >= mock_width || py < 0 || py >= mock_height) continue;
                double dist_sq =
                    ((double)px - u) * ((double)px - u) + ((double)py - v_px) * ((double)py - v_px);
                int noise = (int)(mock_rand_jitter() * 4.0);
                int intensity = (int)floor(peak * exp(-dist_sq / 1.2) + 0.5) + noise;
                if (intensity < 0) intensity = 0;
                int idx = py * mock_width + px;
                int total = (int)mock_buffer[idx] + intensity;
                mock_buffer[idx] = (uint8_t)(total > 255 ? 255 : total);
            }
        }
    }

    /* inject 3 spurious hot-pixel/cosmic-ray blobs per frame to exercise outlier rejection */
    for (int k = 0; k < 3; k++) {
        int hx = 30 + (int)(mock_rand_uniform() * (double)(mock_width - 60));
        int hy = 30 + (int)(mock_rand_uniform() * (double)(mock_height - 60));
        for (int dy = -1; dy <= 1; dy++) {
            for (int dx = -1; dx <= 1; dx++) {
                int idx = (hy + dy) * mock_width + (hx + dx);
                mock_buffer[idx] = (uint8_t)(dx == 0 && dy == 0 ? 180 : 95);
            }
        }
    }

    *buffer = mock_buffer;
    return 0;
}

void camera_close(void) {
    free(mock_buffer);
    mock_buffer = NULL;
    mock_width = 0;
    mock_height = 0;
}
