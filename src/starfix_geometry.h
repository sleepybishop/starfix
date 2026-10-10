#ifndef STARFIX_GEOMETRY_H
#define STARFIX_GEOMETRY_H

#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct {
    double x, y, z;
} starfix_vector3_t;

static inline double starfix_dot(starfix_vector3_t a, starfix_vector3_t b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

static inline double starfix_row_dot(const double row[3], starfix_vector3_t v) {
    return row[0] * v.x + row[1] * v.y + row[2] * v.z;
}

static inline starfix_vector3_t starfix_rotate(const double R[3][3], starfix_vector3_t v) {
    starfix_vector3_t result = {starfix_row_dot(R[0], v), starfix_row_dot(R[1], v),
                                starfix_row_dot(R[2], v)};
    return result;
}

/* RA and dec in radians, in the catalog frame */
static inline starfix_vector3_t starfix_celestial_ray(double ra, double dec) {
    starfix_vector3_t result = {cos(dec) * cos(ra), cos(dec) * sin(ra), sin(dec)};
    return result;
}

/* horizontal FOV in degrees; focal length in pixels */
static inline double starfix_focal_length(double width, double fov_deg) {
    return width / (2.0 * tan(fov_deg * M_PI / 360.0));
}

/* centered pixel offsets to a unit ray: +X right, +Y down, +Z forward */
static inline starfix_vector3_t starfix_camera_ray(double u, double v, double focal) {
    double x = u / focal, y = v / focal;
    double norm = hypot(hypot(x, y), 1.0);
    starfix_vector3_t result = {x / norm, y / norm, 1.0 / norm};
    return result;
}

#endif
