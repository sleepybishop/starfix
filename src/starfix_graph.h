#ifndef STARFIX_GRAPH_H
#define STARFIX_GRAPH_H

#include "starfix_arena.h"
#include "starfix_status.h"

#ifdef __cplusplus
extern "C" {
#endif

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

/* smooth num_nodes intervals with gauss-newton; poses has num_nodes + 1 entries.
 * odom has num_nodes * steps_per_node readings spaced dt_sec seconds apart.
 * initial lat/lon are degrees; scale (measured/true speed) and bias (radians) are in/out.
 * returns zero or a negative error */
starfix_status_t starfix_solve_graph(int num_nodes, int steps_per_node, double dt_sec,
                                     double initial_lat, double initial_lon,
                                     const starfix_odometry_t* odom, int num_fixes,
                                     const starfix_fix_t* fixes, starfix_pose_t* poses,
                                     double* est_scale, double* est_bias, starfix_arena_t arena);

#ifdef __cplusplus
}
#endif
#endif /* STARFIX_GRAPH_H */
