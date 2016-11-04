/*----------------------------------------------------------------------------*/
/** @file lsd_square.c
    @author Chiu Yue Chun [aka BrianOn99] <chiu6700@gmail.com>
 */
/*----------------------------------------------------------------------------*/

#ifdef DEBUG
# include <stdio.h>
# define DEBUG_PRINT(...) do{ fprintf( stderr, __VA_ARGS__ ); } while( false )
#else
# define DEBUG_PRINT(...) do{ } while ( false )
#endif

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <stdbool.h>
#include "comm_math.h"
#include "lsd.h"
#include "lsd_square.h"

#define COS_1_OVER_25_PI 0.9921
#define COS_1_OVER_8_PI 0.9239
#define COS_3_OVER_8_PI 0.3827

enum line_class {TOP = 0, RIGHT, BOTTOM, LEFT, OTHERS};

struct vector horizontal_vector = { .x = 1, .y = 0 };

struct classified_lines {
	struct line *line_group[4];
	unsigned count[4];
	double side_len;
};

void init_classified_lines(struct classified_lines *cl_lines, int max_count)
{
	unsigned i;
	memset(cl_lines, 0, sizeof(struct classified_lines));
	for (i = 0; i < 4; i++)
		cl_lines->line_group[i] = malloc(max_count * sizeof(struct line));
}

void free_classified_lines(struct classified_lines *cl_lines)
{
	unsigned i;
	for (i = 0; i < 4; i++)
		free(cl_lines->line_group[i]);
}

void add_classified_line(enum line_class class, struct line *l, struct classified_lines *cl_lines)
{
	cl_lines->line_group[class][cl_lines->count[class]] = *l;
	(cl_lines->count[class])++;
}

/*
 * 0 if they are not (almost) parallel, else their dot product
 * In other words, we treat the 2 lines as perpendicular if they have "some angle"
 */
double vectors_parallel(struct vector *v1, struct vector *v2) {
	double dot_prod;
	dot_prod = dot_product(v1, v2);
	return (fabs(dot_prod) < COS_1_OVER_25_PI) ? 0.0 : dot_prod;
}

/*
 * join adjacent lines which seems to belongs to 1 line
 * in_lines should be 7tuple return by lsd() family.  It will be modified.
 */
struct line* join_lines(unsigned in_count, struct line  out_lines[], double in_lines[], unsigned *out_count, double thres_dist)
{
	unsigned i, j;
	unsigned num_items;  // how many items in in_lines (not removed).  It is <= in_count
	int num_joined;  // number of processed and joined lines
	struct line current_line;
	struct line target_line;
	struct vector *unit_vectors;

	unit_vectors = malloc(sizeof(struct vector) * in_count);

	num_joined = 0;
	num_items = in_count;

	if (out_lines == NULL || unit_vectors == NULL)
		return NULL;

	for (i = 0; i < in_count; i++) {
		double *tup;
		tup = in_lines + i * 7;
		make_unit_vector(unit_vectors + i,  tup[2] - tup[0], tup[3] - tup[1]);
	}

	for (i = 0; i < num_items; i++) {
		double *tup;
		tup = in_lines + i * 7;
		current_line.x1 = tup[0];
		current_line.y1 = tup[1];
		current_line.x2 = tup[2];
		current_line.y2 = tup[3];

redo:
		for (j = i + 1; j < num_items; j++) {
			double test_res;
			double d1, d2;
			double *inner_tup;
			struct mag_vector joining_vector;
			struct line candidate;

			test_res = vectors_parallel(&unit_vectors[i], &unit_vectors[j]);
			if (!test_res)
				continue;

			DEBUG_PRINT("!! %d %d passed parallel test\n", i, j);

			inner_tup = in_lines + j * 7;

			if (test_res < 0) {
				/* align line direction if not pointing same way */
				target_line.x1 = inner_tup[2];
				target_line.y1 = inner_tup[3];
				target_line.x2 = inner_tup[0];
				target_line.y2 = inner_tup[1];
			} else {
				target_line.x1 = inner_tup[0];
				target_line.y1 = inner_tup[1];
				target_line.x2 = inner_tup[2];
				target_line.y2 = inner_tup[3];
			}

			d1 = dist(current_line.x2, current_line.y2, target_line.x1, target_line.y1);
			d2 = dist(current_line.x1, current_line.y1, target_line.x2, target_line.y2);

			if (d1 > d2) {
				candidate = (struct line) {current_line.x2, current_line.y2,
								target_line.x1, target_line.y1};
				make_mag_vector(&joining_vector, target_line.x2 - current_line.x1,
								target_line.y2 - current_line.y1);
			} else {
				candidate = (struct line) {current_line.x1, current_line.y1,
								target_line.x2, target_line.y2};
				make_mag_vector(&joining_vector, target_line.x1 - current_line.x2,
								target_line.y1 - current_line.y2);
			}

			DEBUG_PRINT("!!  candidate %f %f %f %f\n",
					candidate.x1, candidate.y1, candidate.x2, candidate.y2);

                        if (joining_vector.magnitude > thres_dist ||
				!vectors_parallel(&joining_vector.unit_vec, &unit_vectors[j]))
				continue;

			DEBUG_PRINT("!!! passed interpolate test\n");

			/* Here, the 2 lines are joinable */
			current_line = candidate;

			/* remove the original line */
			if (j != num_items -1) {
				unit_vectors[j] = unit_vectors[num_items-1];
				memcpy(inner_tup, in_lines + (num_items-1)*7, sizeof(double) * 7);
			}
			num_items--;
			goto redo;
		}
		num_joined++;
		out_lines[i] = current_line;
	}

	*out_count = num_joined;

	free(unit_vectors);
	return out_lines;
}

enum line_class classify_line(struct line *l, int img_x, int img_y)
{
	struct mag_vector vec;
	make_mag_vector_by_line(&vec, l);
	if (vec.magnitude < img_x/6)
		return OTHERS;

	double dot_prod = dot_product(&(vec.unit_vec), &horizontal_vector);
	if (fabs(dot_prod) > COS_1_OVER_8_PI) {
		/* horizontal */
		double y_center = (l->y1 + l->y2) / 2;
		if (y_center > img_y * 3/5) {
			if  (l->x1 > l->x2)
				reverse_line(l);
			return TOP;
		} else if (y_center < img_y * 2/5) {
			if  (l->x2 > l->x1)
				reverse_line(l);
			return BOTTOM;
		}
	} else if (fabs(dot_prod) < COS_3_OVER_8_PI) {
		/* vertical */
		double x_center = (l->x1 + l->x2) / 2;
		if (x_center > img_x * 3/5) {
			if  (l->y2 > l->y1)
				reverse_line(l);
			return RIGHT;
		} else if (x_center < img_x * 2/5) {
			if  (l->y1 > l->y2)
				reverse_line(l);
			return LEFT;
		}
	}

	return OTHERS;
}

/*
 * some ad-hoc rule to score how similar the lines are to a rectangle.
 * lower is better
 * gap between lines should be small (positive score)
 * area should be large (negative score) (calculated approximately)
 */
double score_square(struct line lines[4], int side_len)
{
	int i;
	double approx_width, approx_height;
	double score = 0;
	for (i = 0; i < 4; i++) {
		struct line *l1 = &lines[i];
		struct line *l2 = &lines[(i+1) % 4];
		score += dist(l1->x2, l1->y2, l2->x1, l2->y1) * (4.0/side_len);
	}
	approx_height = ( (lines[0].y1 + lines[0].y2) - (lines[2].y1 + lines[2].y2) ) / 2;
	approx_width = ( (lines[1].x1 + lines[1].x2) - (lines[3].x1 + lines[3].x2) ) / 2;
	score -= sqrt(approx_height * approx_width);
	return score;
}

void res_filter_best_square(struct scored_square *res, struct line test_square[4], struct classified_lines *cl_lines, int rec_level)
{
	unsigned i;
	for (i=0; i < cl_lines->count[rec_level]; i++) {
		DEBUG_PRINT("%*s recursion i=%d\n", rec_level, "", i);

		test_square[rec_level] = cl_lines->line_group[rec_level][i];
		if (rec_level != 3) {
			res_filter_best_square(res, test_square, cl_lines, rec_level + 1);
		} else {
			double score = score_square(test_square, cl_lines->side_len);
			DEBUG_PRINT("%*s score %7.5f\n", rec_level, "", score);

			if (score < res->score) {
				memcpy(res->lines, test_square, sizeof(struct line) * 4);
				res->score = score;
			}
		}
	}
}

/* find square with brute force */
void filter_best_square(struct scored_square *res, struct classified_lines *cl_lines)
{
	struct line test_square[4];
	res_filter_best_square(res, test_square, cl_lines, 0);
}

void find_square(struct scored_square *res_square, struct line *lines, unsigned n, int X, int Y)
{
	unsigned i;
	struct classified_lines cl_lines;
	res_square->score = INT_MAX;

	init_classified_lines(&cl_lines, n);
	cl_lines.side_len = X;

	for (i = 0; i < n; i++) {
		enum line_class cls = classify_line(&lines[i], X, Y);
		DEBUG_PRINT("found line type %d\n", cls);
		if (cls != OTHERS)
			add_classified_line(cls, &lines[i], &cl_lines);
	}

	filter_best_square(res_square, &cl_lines);
	free_classified_lines(&cl_lines);
}

int find_square_corner(struct point_d res_pt[4], struct line *lines, unsigned n, int X, int Y)
{
	struct scored_square best_square;
	struct line *sq_lines;

	find_square(&best_square, lines, n, X, Y);
	if (best_square.score != INT_MAX) {
		sq_lines = best_square.lines;
		line_intersection(res_pt, sq_lines, sq_lines+1);
		line_intersection(res_pt + 1, sq_lines + 1, sq_lines + 2);
		line_intersection(res_pt + 2, sq_lines + 2, sq_lines + 3);
		line_intersection(res_pt + 3, sq_lines + 3, sq_lines);
		return true;
	} else {
		return false;
	}
}

/*
 * A wrapper to find_square, for bitmap input
 */
int find_square_corner_bitmap(struct point_d res_pt[4], unsigned char *image, int X, int Y)
{
	int n;
	unsigned m;
	struct line *joined_lines;

	double *segs = lsd(&n, image, X, Y);
	joined_lines = malloc(sizeof(struct line) * n);
	join_lines(n, joined_lines, segs, &m, X/8);
	int is_success = find_square_corner(res_pt, joined_lines, m, X, Y);
	free(joined_lines);
	free(segs);
	return is_success;
}
