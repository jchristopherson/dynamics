#include "dynamics_geometry_tests.h"
#include "dynamics_c_test_helper.h"
#include "dynamics.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

bool c_test_line_from_2_points()
{
    bool rst;
    int i;
    double pt1[3], pt2[3], v[3];
    c_line ln;
    const double tol = 1.0e-8;

    // Initialization
    rst = true;
    create_random_vector(3, pt1);
    create_random_vector(3, pt2);
    for (i = 0; i < 3; ++i) v[i] = pt2[i] - pt1[i];
    c_line_from_2_points(pt1, pt2, &ln);

    // Tests
    if (!compare_arrays(3, pt1, ln.r0, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_line_from_2_points -1\n");
    }
    if (!compare_arrays(3, v, ln.v, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_line_from_2_points -2\n");
    }
    return rst;
}

bool c_test_line_from_2_planes()
{
    bool rst;
    double n1[3], n2[3], v[3], pt[3];
    c_plane p1, p2;
    c_line ln;
    const double tol = 1.0e-8;

    // Initialization
    rst = true;
    n1[0] = 1.0;    n1[1] = 1.0;    n1[2] = 0.0;
    n2[0] = 0.0;    n2[1] = 0.0;    n2[2] = 1.0;
    v[0] = 1.0;     v[1] = -1.0;    v[2] = 0.0;
    pt[0] = 0.0;    pt[1] = 0.0;    pt[2] = 0.0;
    c_vector_normalize(3, n1);
    c_vector_normalize(3, n2);
    c_vector_normalize(3, v);
    c_plane_from_point_and_normal(pt, n1, &p1);
    c_plane_from_point_and_normal(pt, n2, &p2);
    c_line_from_2_planes(&p1, &p2, &ln);

    // Tests
    if (!compare_arrays(3, pt, ln.r0, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_line_from_2_planes -1\n");
    }
    if (!compare_arrays(3, v, ln.v, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_line_from_2_planes -2\n");
    }
    return rst;
}

bool c_test_line_eval()
{
    bool rst;
    int i;
    double t, ans[3], s[3];
    c_line ln;
    const double tol = 1.0e-8;

    // Initialization
    rst = true;
    create_random_vector(1, &t);
    create_random_vector(3, ln.v);
    create_random_vector(3, ln.r0);
    for (i = 0; i < 3; ++i) ans[i] = ln.r0[i] + t * ln.v[i];
    c_evaluate_line_position(&ln, t, s);

    // Test
    if (!compare_arrays(3, ans, s, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_line_eval -1\n");
    }
    return rst;
}

bool c_test_plane_from_3_points()
{
    bool rst;
    int i;
    double pt1[3], pt2[3], pt3[3], nrm[3], pt21[3], pt31[3];
    c_plane pln;
    const double tol = 1.0e-8;

    // Initialization
    rst = true;
    create_random_vector(3, pt1);
    create_random_vector(3, pt2);
    create_random_vector(3, pt3);
    for (i = 0; i < 3; ++i)
    {
        pt21[i] = pt2[i] - pt1[i];
        pt31[i] = pt3[i] - pt1[i];
    }
    c_cross_product(pt21, pt31, nrm);
    c_vector_normalize(3, nrm);
    c_plane_from_3_points(pt1, pt2, pt3, &pln);

    // Test
    if (fabs(pln.a - nrm[0]) > tol)
    {
        rst = false;
        printf("TEST FAILED: c_test_plane_from_3_points -1\n");
    }
    if (fabs(pln.b - nrm[1]) > tol)
    {
        rst = false;
        printf("TEST FAILED: c_test_plane_from_3_points -2\n");
    }
    if (fabs(pln.c - nrm[2]) > tol)
    {
        rst = false;
        printf("TEST FAILED: c_test_plane_from_3_points -3\n");
    }
    if (fabs(pln.a * pt1[0] + pln.b * pt1[1] + pln.c * pt1[2] + pln.d) > tol)
    {
        rst = false;
        printf("TEST FAILED: c_test_plane_from_3_points -4\n");
    }

    return rst;
}

bool c_test_plane_from_point_and_normal()
{
    bool rst;
    double pt1[3], nrm[3];
    c_plane pln;
    const double tol = 1.0e-8;

    // Initialization
    rst = true;
    create_random_vector(3, nrm);
    create_random_vector(3, pt1);
    c_vector_normalize(3, nrm);
    c_plane_from_point_and_normal(pt1, nrm, &pln);

    // Test
    if (fabs(pln.a - nrm[0]) > tol)
    {
        rst = false;
        printf("TEST FAILED: test_plane_from_point_and_normal -1\n");
    }
    if (fabs(pln.b - nrm[1]) > tol)
    {
        rst = false;
        printf("TEST FAILED: test_plane_from_point_and_normal -2\n");
    }
    if (fabs(pln.c - nrm[2]) > tol)
    {
        rst = false;
        printf("TEST FAILED: test_plane_from_point_and_normal -3\n");
    }
    if (fabs(pln.a * pt1[0] + pln.b * pt1[1] + pln.c * pt1[2] + pln.d) > tol)
    {
        rst = false;
        printf("TEST FAILED: test_plane_from_point_and_normal -4\n");
    }
    return rst;
}

bool c_test_plane_normal()
{
    bool rst;
    double pt1[3], nrm[3], pnrm[3];
    c_plane pln;
    const double tol = 1.0e-8;

    // Initialization
    rst = true;
    create_random_vector(3, pt1);
    create_random_vector(3, nrm);
    c_plane_from_point_and_normal(pt1, nrm, &pln);
    c_vector_normalize(3, nrm);
    c_plane_normal(&pln, pnrm);

    // Test
    if (!compare_arrays(3, nrm, pnrm, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_plane_normal -1\n");
    }
    return rst;
}

bool c_test_is_parallel_vectors()
{
    bool rst;
    int i;
    double p1[3], p2[3], np1[3], np2[3];
    const double tol = 1.0e-8;

    // Initialization
    rst = true;
    create_random_vector(3, p1);
    create_random_vector(3, np1);
    create_random_vector(3, np2);
    for (i = 0; i < 3; ++i) p2[i] = 3.0 * p1[i];  // ensuring p1 and p2 are parallel
    np2[2] += 1.5;  // ensuring np1 and np2 are not parallel ever

    if (!c_is_parallel_vectors(3, p1, p2, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_is_parallel_vectors -1\n");
    }
    if (c_is_parallel_vectors(3, np1, np2, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_is_parallel_vectors -2\n");
    }
    return rst;
}

bool c_test_is_parallel_lines()
{
    bool rst;
    int i;
    c_line ln1, ln2, ln3, ln4;
    const double tol = 1.0e-8;

    // Initialization
    rst = true;
    create_random_vector(3, ln1.r0);
    create_random_vector(3, ln1.v);
    create_random_vector(3, ln3.r0);
    create_random_vector(3, ln3.v);
    create_random_vector(3, ln4.r0);
    create_random_vector(3, ln4.v);
    for (i = 0; i < 3; ++i)
    {
        ln2.r0[i] = ln1.r0[i];
        ln2.v[i] = 3.0 * ln1.v[i];
    }
    ln4.v[2] += 1.5;

    // Test
    if (!c_is_parallel_lines(&ln1, &ln2, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_is_parallel_lines -1\n");
    }
    if (c_is_parallel_lines(&ln3, &ln4, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_is_parallel_lines -2\n");
    }
    return rst;
}

bool c_test_is_parallel_planes()
{
    bool rst;
    int i;
    double p1[3], p2[3], np1[3], np2[3], pt[3];
    c_plane pln1, pln2, pln3, pln4;
    const double tol = 1.0e-8;

    // Initialization
    rst = true;
    create_random_vector(3, p1);
    create_random_vector(3, np1);
    create_random_vector(3, np2);
    create_random_vector(3, pt);
    for (i = 0; i < 3; ++i) p2[i] = 3.0 * p1[i];
    np2[2] += 3.0;
    c_plane_from_point_and_normal(pt, p1, &pln1);
    c_plane_from_point_and_normal(pt, p2, &pln2);
    c_plane_from_point_and_normal(pt, np1, &pln3);
    c_plane_from_point_and_normal(pt, np2, &pln4);

    // Tests
    if (!c_is_parallel_planes(&pln1, &pln2, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_is_parallel_planes -1\n");
    }
    if (c_is_parallel_planes(&pln3, &pln4, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_is_parallel_planes -2\n");
    }
    return rst;
}

bool c_test_is_point_on_plane()
{
    bool rst;
    int i;
    double pt0[3], pt1[3], pt2[3], nrm[3];
    c_plane pln;
    const double tol = 1.0e-8;

    // Initialization
    rst = true;
    create_random_vector(3, pt0);
    create_random_vector(3, pt1);
    create_random_vector(3, nrm);
    c_plane_from_point_and_normal(pt0, nrm, &pln);

    // Define a point on the plane and a point off the plane
    pt1[2] = -(pln.a * pt1[0] + pln.b * pt1[1] + pln.d) / pln.c;
    for (i = 0; i < 3; ++i) pt2[i] = pt1[i] + 0.5;

    // Test
    if (!c_is_point_on_plane(pt1, &pln, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_is_point_on_plane -1\n");
    }
    if (c_is_point_on_plane(pt2, &pln, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_is_point_on_plane -2\n");
    }
    return rst;
}

bool c_test_point_to_line_distance()
{
    bool rst;
    int i;
    double pt1[3], pt2[3], pt[3], t, dist, d;
    c_line ln;
    const double t_ans = 0.5;
    const double tol = 1.0e-8;

    // Initialization
    rst = true;
    create_random_vector(1, &dist);
    pt1[0] = 0.0;   pt1[1] = 0.0;   pt1[2] = 0.0;
    pt2[0] = 1.0;   pt2[1] = 1.0;   pt2[2] = 0.0;
    for (i = 0; i < 3; ++i) pt[i] = t_ans * pt2[i];
    pt[2] = dist;
    c_line_from_2_points(pt1, pt2, &ln);

    // Find the nearest point on the line
    t = c_nearest_point_on_line(pt, &ln);
    if (fabs(t - t_ans) > tol)
    {
        rst = false;
        printf("TEST FAILED: c_test_point_to_line_distance -1\n");
    }

    // Compute the distance
    d = c_point_to_line_distance(pt, &ln);
    if (fabs(dist - d) > tol)
    {
        rst = false;
        printf("TEST FAILED: c_test_point_to_line_distance -2\n");
    }
    return rst;
}

bool c_test_is_point_on_line()
{
    bool rst;
    int i;
    double pt1[3], pt2[3], ptOn[3], ptOff[3], t;
    c_line ln;
    const double tol = 1.0e-8;

    // Initialization
    rst = true;
    create_random_vector(1, &t);
    create_random_vector(3, pt1);
    create_random_vector(3, pt2);
    c_line_from_2_points(pt1, pt2, &ln);
    c_evaluate_line_position(&ln, t, ptOn);
    for (i = 0; i < 3; ++i) ptOff[i] = ptOn[i] + 0.5;

    // Test
    if (!c_is_point_on_line(ptOn, &ln, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_is_point_on_line -1\n");
    }
    if (c_is_point_on_line(ptOff, &ln, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_is_point_on_line -2\n");
    }
    return rst;
}

bool c_test_point_to_plane_distance()
{
    bool rst;
    int i;
    double pt[3], pt0[3], nrm[3], dist, d;
    c_plane pln;
    const double tol = 1.0e-8;

    // Initialization
    rst = true;
    create_random_vector(3, pt0);
    create_random_vector(3, nrm);
    c_vector_normalize(3, nrm);
    create_random_vector(1, &dist);
    for (i = 0; i < 3; ++i) pt[i] = pt0[i] + dist * nrm[i];
    c_plane_from_point_and_normal(pt0, nrm, &pln);

    // Test
    d = c_point_to_plane_distance(pt, &pln);
    if (fabs(d - dist) > tol)
    {
        rst = false;
        printf("TEST FAILED: c_test_point_to_plane_distance -1\n");
    }
    return rst;
}

bool c_test_vector_plane_projection()
{
    bool rst;
    int i;
    double v[3], pt[3], nrm[3], prj[3], ans[3], vn, nn;
    c_plane pln;
    const double tol = 1.0e-8;

    // Initialization
    rst = true;
    pt[0] = 0.0;    pt[1] = 0.0;    pt[2] = 0.0;
    nrm[0] = 1.0;   nrm[1] = 1.0;   nrm[2] = 1.0;
    c_plane_from_point_and_normal(pt, nrm, &pln);
    v[0] = 1.0;     v[1] = 0.0;     v[2] = 0.0;
    vn = c_dot_product(3, v, nrm);
    nn = c_dot_product(3, nrm, nrm);
    for (i = 0; i < 3; ++i) ans[i] = v[i] - nrm[i] * vn / nn;
    c_vector_plane_projection(v, &pln, prj);

    if (!compare_arrays(3, ans, prj, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_vector_plane_projection -1\n");
    }
    return rst;
}

bool c_test_fit_line_to_many_points()
{
    bool rst;
    double pts[9], v[3];
    c_line ln;
    const double tol = 1.0e-8;

    // Initialization
    rst = true;
    pts[0] = 0.0;       pts[3] = 0.0;       pts[6] = 0.0;
    pts[1] = 0.5;       pts[4] = 0.5;       pts[7] = 0.5;
    pts[2] = 1.0;       pts[5] = 1.0;       pts[8] = 1.0;
    v[0] = 1.0;         v[1] = 1.0;         v[2] = 1.0;
    c_vector_normalize(3, v);
    c_line_from_points(3, pts, 3, &ln);

    // Test
    c_vector_normalize(3, ln.v);
    if (!compare_arrays(3, v, ln.v, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_fit_line_to_many_points -1\n");
    }
    return rst;
}

bool c_test_fit_plane_to_many_points()
{
    bool rst;
    double pts[12], nrm[3], pn[3];
    c_plane pln;
    const double tol = 1.0e-8;

    // Initialization
    rst = true;
    pts[0] = 0.0;       pts[4] = 0.0;       pts[8] = 0.0;
    pts[1] = 1.0;       pts[5] = 1.0;       pts[9] = 0.0;
    pts[2] = 1.0;       pts[6] = 1.0;       pts[10] = 1.0;
    pts[3] = 0.0;       pts[7] = 0.0;       pts[11] = 1.0;
    nrm[0] = 1.0;       nrm[1] = -1.0;      nrm[2] = 0.0;
    c_vector_normalize(3, nrm);
    c_plane_from_points(4, pts, 4, &pln);

    // Test
    c_plane_normal(&pln, pn);
    if (!compare_arrays(3, pn, nrm, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_fit_plane_to_many_points -1\n");
    }
    return rst;
}

bool c_test_flip_plane_normal()
{
    bool rst;
    int i;
    double pt[3], nrm[3], n[3], ans[3];
    c_plane pln;
    const double tol = 1.0e-8;

    // Initialization
    rst = true;
    create_random_vector(3, pt);
    create_random_vector(3, nrm);
    c_vector_normalize(3, nrm);
    c_plane_from_point_and_normal(pt, nrm, &pln);
    c_flip_plane_normal(&pln);
    c_plane_normal(&pln, n);
    for (i = 0; i < 3; ++i) ans[i] = -nrm[i];

    // Test
    if (!compare_arrays(3, n, ans, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_flip_plane_normal -1\n");
    }
    return rst;
}
