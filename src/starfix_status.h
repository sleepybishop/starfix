#ifndef STARFIX_STATUS_H
#define STARFIX_STATUS_H

#include <stdint.h>

/* starfix_status: Comprehensive status codes for starfix telemetry */

typedef enum {
    /* Success */
    STARFIX_SUCCESS = 0,

    /* General Errors */
    STARFIX_ERR_NULL_POINTER = -1,
    STARFIX_ERR_MEMORY = -2,
    STARFIX_ERR_INVALID_PARAM = -3,
    STARFIX_ERR_DATA_CORRUPTION = -4, /* NaN or Inf detected */

    /* Centroiding Errors */
    STARFIX_ERR_CENTROID_SATURATION = -10,
    STARFIX_ERR_CENTROID_TOO_FEW = -11,
    STARFIX_ERR_CENTROID_BUFFER_FULL = -12,

    /* Identification Errors */
    STARFIX_ERR_IDENTIFY_NO_MATCH = -20,
    STARFIX_ERR_IDENTIFY_AMBIGUOUS = -21,

    /* Attitude Errors */
    STARFIX_ERR_ATTITUDE_DEGENERATE = -30, /* Collinear vectors */
    STARFIX_ERR_RANSAC_INSUFFICIENT_INLIERS = -31,

    /* Solver / Graph Errors */
    STARFIX_ERR_SOLVER_DIVERGENCE = -40,
    STARFIX_ERR_SOLVER_SINGULAR = -41,

    /* EKF Errors */
    STARFIX_ERR_EKF_DIVERGENCE = -50
} starfix_status_t;

/* starfix_telemetry: Diagnostic metrics for system health monitoring */
typedef struct {
    /* Centroiding Metrics */
    uint32_t num_stars_detected;
    uint32_t saturated_pixels;

    /* ID Metrics */
    uint32_t hash_collisions;
    uint32_t identify_matches;

    /* Attitude Metrics */
    uint32_t ransac_iterations;
    uint32_t ransac_inliers;
    double ransac_inlier_ratio;

    /* Solver Metrics */
    uint32_t solver_iterations;
    double final_residual;

    /* EKF / Graph Metrics */
    double ekf_innovation_norm;
} starfix_telemetry_t;

#endif /* STARFIX_STATUS_H */

#include <math.h>
#define STARFIX_ASSERT_VALID_DOUBLE(val)        \
    do {                                        \
        if (isnan(val) || isinf(val)) {         \
            return STARFIX_ERR_DATA_CORRUPTION; \
        }                                       \
    } while (0)
