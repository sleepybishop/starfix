#ifndef STARFIX_STATUS_H
#define STARFIX_STATUS_H

#include <stdint.h>

typedef enum {
    STARFIX_SUCCESS = 0,

    STARFIX_ERR_NULL_POINTER = -1,
    STARFIX_ERR_MEMORY = -2,
    STARFIX_ERR_INVALID_PARAM = -3,
    STARFIX_ERR_DATA_CORRUPTION = -4,

    STARFIX_ERR_CENTROID_BUFFER_FULL = -12,

    STARFIX_ERR_IDENTIFY_NO_MATCH = -20,

    STARFIX_ERR_ATTITUDE_DEGENERATE = -30, /* Collinear vectors */
    STARFIX_ERR_RANSAC_INSUFFICIENT_INLIERS = -31,

    STARFIX_ERR_SOLVER_DIVERGENCE = -40,
    STARFIX_ERR_SOLVER_SINGULAR = -41,

    STARFIX_ERR_EKF_DIVERGENCE = -50
} starfix_status_t;

typedef struct {
    uint32_t solver_iterations;
    double ekf_innovation_norm;
} starfix_telemetry_t;

#endif /* STARFIX_STATUS_H */
