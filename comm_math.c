#include <math.h>
#include <stdbool.h>
#include <float.h>
#include <assert.h>
#include "comm_math.h"

/*----------------------------------------------------------------------------*/
/** @file comm_math.c
    common math utilities
    @author rafael grompone von gioi <grompone@gmail.com>
    @author Chiu Yue Chun [aka BrianOn99] <chiu6700@gmail.com>
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/** Computes Euclidean distance between point (x1,y1) and point (x2,y2).
 */
double dist(double x1, double y1, double x2, double y2)
{
	return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}

/*----------------------------------------------------------------------------*/
/** Absolute value angle difference.
 */
double angle_diff(double a, double b)
{
	a -= b;
	while (a <= -M_PI)
		a += M_2__PI;
	while (a > M_PI)
		a -= M_2__PI;
	if (a < 0.0)
		a = -a;
	return a;
}

/*----------------------------------------------------------------------------*/
/** Signed angle difference.
 */
double angle_diff_signed(double a, double b)
{
	a -= b;
	while (a <= -M_PI)
		a += M_2__PI;
	while (a > M_PI)
		a -= M_2__PI;
	return a;
}

/*----------------------------------------------------------------------------*/
/** Interpolate y value corresponding to 'x' value given, in
    the line 'x1,y1' to 'x2,y2'; if 'x1=x2' return the smaller
    of 'y1' and 'y2'.

    The following restrictions are required:
    - x1 <= x2
    - x1 <= x
    - x  <= x2
 */
double inter_low(double x, double x1, double y1, double x2, double y2)
{
	/* check parameters */
	assert(x1 <= x2 && x1 <= x && x <= x2);

	/* interpolation */
	if (double_equal(x1, x2) && y1 < y2)
		return y1;
	if (double_equal(x1, x2) && y1 > y2)
		return y2;
	return y1 + (x - x1) * (y2 - y1) / (x2 - x1);
}

/*----------------------------------------------------------------------------*/
/** Interpolate y value corresponding to 'x' value given, in
    the line 'x1,y1' to 'x2,y2'; if 'x1=x2' return the larger
    of 'y1' and 'y2'.

    The following restrictions are required:
    - x1 <= x2
    - x1 <= x
    - x  <= x2
 */
double inter_hi(double x, double x1, double y1, double x2, double y2)
{
	/* check parameters */
	assert(x1 <= x2 && x1 <= x && x <= x2);

	/* interpolation */
	if (double_equal(x1, x2) && y1 < y2)
		return y2;
	if (double_equal(x1, x2) && y1 > y2)
		return y1;
	return y1 + (x - x1) * (y2 - y1) / (x2 - x1);
}

/*----------------------------------------------------------------------------*/
/** Doubles relative error factor
 */
#define RELATIVE_ERROR_FACTOR 100.0

/*----------------------------------------------------------------------------*/
/** Compare doubles by relative error.

    The resulting rounding error after floating point computations
    depend on the specific operations done. The same number computed by
    different algorithms could present different rounding errors. For a
    useful comparison, an estimation of the relative rounding error
    should be considered and compared to a factor times EPS. The factor
    should be related to the cumulated rounding error in the chain of
    computation. Here, as a simplification, a fixed factor is used.
 */
int double_equal(double a, double b)
{
	double abs_diff, aa, bb, abs_max;

	/* trivial case */
	if (a == b)
		return true;

	abs_diff = fabs(a - b);
	aa = fabs(a);
	bb = fabs(b);
	abs_max = aa > bb ? aa : bb;

	/* DBL_MIN is the smallest normalized number, thus, the smallest
	   number whose relative error is bounded by DBL_EPSILON. For
	   smaller numbers, the same quantization steps as for DBL_MIN
	   are used. Then, for smaller numbers, a meaningful "relative"
	   error should be computed by dividing the difference by DBL_MIN. */
	if (abs_max < DBL_MIN)
		abs_max = DBL_MIN;

	/* equal if relative error <= factor x eps */
	return (abs_diff / abs_max) <= (RELATIVE_ERROR_FACTOR * DBL_EPSILON);
}

void make_mag_vector_by_line(struct mag_vector *vec, struct line* l)
{
	make_mag_vector(vec, l->x2 - l->x1, l->y2 - l->y1);
}

void make_unit_vector(struct vector *unit_vector, double x_len, double y_len)
{
	double magnitude = sqrt(y_len * y_len + x_len * x_len);
	unit_vector->x = x_len / magnitude;
	unit_vector->y = y_len / magnitude;
}

void make_mag_vector(struct mag_vector *vec, double x_len, double y_len)
{
	double magnitude = sqrt(y_len * y_len + x_len * x_len);
	vec->unit_vec.x = x_len / magnitude;
	vec->unit_vec.y = y_len / magnitude;
	vec->magnitude = magnitude;
}

void reverse_line(struct line *l)
{
	double tmp;
	tmp = l->x1;
	l->x1 = l->x2;
	l->x2 = tmp;
	tmp = l->y1;
	l->y1 = l->y2;
	l->y2 = tmp;
}

/*
 * Inner product (aka dot product) of 2D vectors.
 */
inline double dot_product(struct vector *v1, struct vector *v2)
{
	return v1->x * v2->x + v1->y * v2->y;
}

/*
 * 2x2 determinant
 */
double det2by2(double M11, double M12, double M21, double M22)
{
	return M11 * M22 - M12 * M21;
}

/*
 * From https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection
 * y3 y4 are the y cordinate of the ends of line2
 * |...| is determinant
 *
 *        x-coordinate           y-coordinaate
 *    ||x1  y1|  |x1   1||    ||x1  y1| |y1   1||
 *    ||x2  y2|  |x2   1||    ||x2  y2| |y2   1||
 *    |                  |    |                 |
 *    ||x3  y3|  |x3   1||    ||x3  y3| |y3   1||
 *    ||x4  y4|  |x4   1||    ||x4  y4| |y4   1||
 *    ---------------------  --------------------
 *    ||x1   1|  |y1   1||    ||x1  y1| |y1   1||
 *    ||x2   1|  |y2   1||    ||x2  y2| |y2   1||
 *    |                  |    |                 |
 *    ||x3   1|  |y3   1||    ||x3  y3| |y3   1||
 *    ||x4   1|  |y4   1||    ||x4  y4| |y4   1||
 */
void line_intersection(struct point_d *res, struct line *l1, struct line *l2)
{
	double a, b, c, d, e, f;
	double dividend;
	a = det2by2(l1->x1, 1, l1->x2, 1);
	b = det2by2(l1->y1, 1, l1->y2, 1);
	c = det2by2(l2->x1, 1, l2->x2, 1);
	d = det2by2(l2->y1, 1, l2->y2, 1);
	e = det2by2(l1->x1, l1->y1, l1->x2, l1->y2);
	f = det2by2(l2->x1, l2->y1, l2->x2, l2->y2);
	dividend = det2by2(a, b, c, d);

	res->x = det2by2(e, a, f, c) / dividend;
	res->y = det2by2(e, b, f, d) / dividend;
}
