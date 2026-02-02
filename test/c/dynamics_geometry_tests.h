#ifndef DYNAMICS_GEOMETRY_TESTS_H_
#define DYNAMICS_GEOMETRY_TESTS_H_

#include <stdbool.h>

bool c_test_line_from_2_points();
bool c_test_line_from_2_planes();
bool c_test_line_eval();
bool c_test_plane_from_3_points();
bool c_test_plane_from_point_and_normal();
bool c_test_plane_normal();
bool c_test_is_parallel_vectors();
bool c_test_is_parallel_lines();
bool c_test_is_parallel_planes();
bool c_test_is_point_on_plane();
bool c_test_point_to_line_distance();
bool c_test_is_point_on_line();
bool c_test_point_to_plane_distance();
bool c_test_vector_plane_projection();
bool c_test_fit_line_to_many_points();
bool c_test_fit_plane_to_many_points();
bool c_test_flip_plane_normal();

#endif
