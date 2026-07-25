#ifndef STARFIX_GRAPH_H
#define STARFIX_GRAPH_H

#include "starfix_arena.h"
#include "starfix_status.h"

/* starfix_graph: C API for factor graph trajectory smoothing. performs batch
   non-linear least-squares trajectory optimization over a set of pose keyframes,
   fusing dead reckoning and celestial fixes, and calibrating sensor biases. */

typedef struct {
    double lat; /* latitude in degrees */
    double lon; /* longitude in degrees */
} starfix_pose_t;

typedef struct {
    double speed;   /* measured speed in knots */
    double heading; /* measured heading in degrees */
} starfix_odometry_t;

typedef struct {
    int node_idx; /* keyframe node index (0 to num_nodes) */
    double lat;   /* celestial latitude in degrees */
    double lon;   /* celestial longitude in degrees */
} starfix_fix_t;

/* solves the factor graph optimization using gauss-newton iterations:
   - num_nodes: number of keyframe poses to solve (excluding start pose)
   - steps_per_node: number of odometry readings integrated per keyframe interval
   - dt_sec: time interval of each odometry reading in seconds
   - initial_lat: true starting latitude in degrees
   - initial_lon: true starting longitude in degrees
   - odom: raw odometry readings (size num_nodes * steps_per_node)
   - num_fixes: number of celestial fixes
   - fixes: array of celestial fixes
   - poses: output array of optimized keyframe poses (size num_nodes + 1)
   - est_scale: input/output pointer to speed log scale factor
   - est_bias: input/output pointer to heading bias in radians
   - returns: number of iterations on success, or a negative value on error */
starfix_status_t starfix_solve_graph(int num_nodes, int steps_per_node, double dt_sec,
                                     double initial_lat, double initial_lon,
                                     const starfix_odometry_t* odom, int num_fixes,
                                     const starfix_fix_t* fixes, starfix_pose_t* poses,
                                     double* est_scale, double* est_bias, starfix_arena_t* arena,
                                     starfix_telemetry_t* telem);

#endif /* STARFIX_GRAPH_H */
