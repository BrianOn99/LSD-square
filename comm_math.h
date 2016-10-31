/** ln(10) */
#ifndef M_LN10
#define M_LN10 2.30258509299404568402
#endif				/* !M_LN10 */

/** PI */
#ifndef M_PI
#define M_PI   3.14159265358979323846
#endif				/* !M_PI */

/** 3/2 pi */
#define M_3_2_PI 4.71238898038

/** 2 pi */
#define M_2__PI  6.28318530718

struct line {
	double x1, y1, x2, y2;
};

struct vector {
	double x, y;
};

struct mag_vector {
	struct vector unit_vec;
	double magnitude;
};

double dist(double x1, double y1, double x2, double y2);
double angle_diff(double a, double b);
double angle_diff_signed(double a, double b);
double inter_low(double x, double x1, double y1, double x2, double y2);
double inter_hi(double x, double x1, double y1, double x2, double y2);
int double_equal(double a, double b);
void make_unit_vector(struct vector *unit_vector, double y_len, double x_len);
double dot_product(struct vector *v1, struct vector *v2);
void make_mag_vector(struct mag_vector *vec, double x_len, double y_len);
