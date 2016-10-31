/*----------------------------------------------------------------------------*/
/** @file lsd_square.c
    @author Chiu Yue Chun [aka BrianOn99] <chiu6700@gmail.com>
 */
/*----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "comm_math.h"

#define COS_1_OVER_25_PI 0.9921

/*
 * 0 if they are not (almost) parallel, else their dot product
 * In other words, we treat the 2 lines as perpendicular if they have "some angle"
 */
double vectors_parallel(struct vector *v1, struct vector *v2) {
	double dot_prod;
	dot_prod = dot_product(v1, v2);
	return (fabs(dot_prod) < COS_1_OVER_50_PI) ? 0.0 : dot_prod;
}

/*
 * join adjacent lines which seems to belongs to 1 line
 * in_lines should be 7tuple return by lsd() family.  It will be modified.
 */
struct line* join_lines(int in_count, double *in_lines, int *out_count, double thres_dist)
{
	int i, j;
	int num_joined;
	int num_items;
	struct line current_line;
	struct line target_line;
	struct line *out_lines;
	struct vector *unit_vectors;

	/* malloc the worst case size */
	out_lines = malloc(sizeof(struct line) * in_count);
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
								target_line.y1 - current_line.x2);
			}

                        if (joining_vector.magnitude > thres_dist ||
				!vectors_parallel(&joining_vector.unit_vec, &unit_vectors[j]))
				continue;

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
