#include "comm_math.h"

#ifndef __LSD_SQUARE_H
#define __LSD_SQUARE_H

struct scored_square {
	struct line lines[4];
	double score;
};

struct line* join_lines(unsigned in_count, struct line  out_lines[], double in_lines[], unsigned *out_count, double thres_dist);
void find_square(struct scored_square *res_square, struct line *lines, unsigned n, int X, int Y);
int find_square_corner(struct point_d res_pt[4], struct line *lines, unsigned n, int X, int Y);
int find_square_corner_bitmap(struct point_d res_pt[4], struct image_char *image, double scale);

#endif
