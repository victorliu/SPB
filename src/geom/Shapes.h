#ifndef _SHAPES_H_
#define _SHAPES_H_

// Shape representations
// =====================
//  Each shape is defined relative to its own local coordinate frame.
// This is prefered for numerical robustness and ease of translation.

// For the purposes of mutual intersection tests, there are only the following types:
//   convex, ellipsoid, frustum
typedef enum{
	SHAPE3_SPHERE,
	SHAPE3_BLOCK, // axes not necessarily perpendicular
	SHAPE3_FRUSTUM, // encompasses cylinder and cone
	SHAPE3_ELLIPSOID, // axes not necessarily perpendicular
	SHAPE3_TETRAHEDRON,
	SHAPE3_CONVEX // convex hull of a set of points
} shape3_type;

typedef struct{
	double r;
} shape3_sphere;

typedef struct{
	double Q[9]; // orthogonal matrix
	double len;
	double r_base;
	double r_tip_base;
	// The last column of Q is the direction of the axis of the frustum.
	// The axis starts at the local origin and ends at the tip cap.
	// Let the matrix T be defined as
	//   T = Q.diag(r_base,r_base,len)
	// Then
	//   inv(T) = inv(diag(len,r_base,r_base)) . Q^T
	// Let r = inv(T).{x,y,z}
	// The frustum is defined by the set
	//   max(
	//     2*norm_inf(r.z - 0.5),
	//     norm_2(r.xy)^2 / (1 + (r_tip_base-1)*r.z)
	//   ) <= 1
} shape3_frustum;

typedef struct{
	double A[9]; // semi axes in columns
	double B[9];
	// The block is defined by the set
	//   norm(B . {x,y,z})_inf = 1
	//   B = inv(A)
} shape3_block;

typedef struct{
	double A[9]; // semi axes in columns
	double B[9];
	// The ellipsoid is defined by the set
	//   (B . {x,y,z})^2 = 1
	//   B = inv(A)
} shape3_ellipsoid;

typedef struct{
	double p[9]; // other points
} shape3_tetrahedron;

typedef struct{
	unsigned int np;
	double *p; // plane points (triplets)
	double *n; // plane normals (triplets)
} shape3_convex;

typedef struct{
	shape3_type type;
	double o[3]; // origin point
	union{
		shape3_sphere sphere;
		shape3_frustum frustum;
		shape3_block block;
		shape3_ellipsoid ellipsoid;
		shape3_tetrahedron tetrahedron;
		shape3_convex convex;
	};
} shape3;



// For the purposes of mutual intersection tests, there are only 2 types:
//   ellipse, arcspline
typedef enum{
	SHAPE2_CIRCLE,
	SHAPE2_ELLIPSE,
	SHAPE2_QUAD,
	SHAPE2_POLYGON,
	SHAPE2_ARCSPLINE
} shape2_type;

typedef struct{
	double r;
} shape2_circle;

// An ellipse may not be an actual ellipse if the axes are not orthogonal
// In this case, the ellipse is the locus of points p such that
// [e[0]/|e[0]|^2 . (p-o)]^2 + [e[1]/|e[1]|^2 . (p-o)]^2 = 1
typedef struct{
	double A[4]; // semi axes of quad in ellipse
	double B[4];
	// The block is defined by the set
	//   norm(B . {x,y,z})_inf = 1
	//   B = inv(A)
} shape2_ellipse;

typedef struct{
	double A[4]; // semi axes of quad in columns
	double B[4];
	// The block is defined by the set
	//   norm(B . {x,y,z})_inf = 1
	// B = inv(A)
} shape2_quad;

typedef struct{
	unsigned int n;
	double *p; // points relative to local origin (pairs)
} shape2_polygon;

typedef struct{
	unsigned int n;
	double *p; // points relative to local origin (pairs)
	double *h; // heights h[i] between p[i] and p[i+1]
} shape2_arcspline;

typedef struct{
	shape2_type type;
	double o[2]; // origin point
	union{
		shape2_circle circle;
		shape2_ellipse ellipse;
		shape2_quad quad;
		shape2_polygon polygon;
		shape2_arcspline arcspline;
	};
} shape2;

typedef struct{
	double c[3]; // center
	double h[3]; // extents
} aabb3;

typedef struct{
	double c[2];
	double h[2];
} aabb2;

// These functions do not check input arguments for NULL.
// All operations return false if not supported

// Expands the box b1 to include the box b2.
void aabb3_union(aabb3 *b1, const aabb3 *b2);
void aabb2_union(aabb2 *b1, const aabb2 *b2);

// Expands the box to include the point.
void aabb3_union_pt(aabb3 *b1, const double p[3]);
void aabb2_union_pt(aabb2 *b1, const double p[2]);

// Determines if p lies within the box. Returns 0 if no, 1 if yes.
int aabb3_contains(const aabb3 *b, const double p[3]);
int aabb2_contains(const aabb2 *b, const double p[2]);

// Determines if two boxes intersect. Returns 0 if no, 1 if yes.
int aabb3_intersects(const aabb3 *a, const aabb3 *b);
int aabb2_intersects(const aabb2 *a, const aabb2 *b);

// Determines if p lies within the shape. Returns 0 if no, 1 if yes.
int shape3_contains(const shape3 *s, const double p[3]);
int shape2_contains(const shape2 *s, const double p[2]); // arcspline not implemented

// Compute the bounding box of the shape.
int shape3_get_aabb(const shape3 *s, aabb3 *b);
int shape2_get_aabb(const shape2 *s, aabb2 *b); // arcspline not implemented

// Determines if two shapes intersect. Returns 0 if no, 1 if yes.
int shape3_intersects(const shape3 *s, const shape3 *t); // not implemented
int shape2_intersects(const shape2 *s, const shape2 *t); // not implemented

// Computes the overlapping volume/area between a shape and
// the given box.
int shape3_aabb_overlap(const shape3 *s, const aabb3 *b); // not implemented
int shape2_aabb_overlap(const shape2 *s, const aabb2 *b); // not implemented

// Computes the overlapping volume/area between a shape and
// the given simplex. The orientation of the simplex determines
// the sign of the returned volume.
double shape3_simplex_overlap(const shape3 *s, const double t[3*4]); // not implemented
double shape2_simplex_overlap(const shape2 *s, const double t[2*3]); // not implemented

// Returns an approximate outward normal vector to the shape at the point
// given by p. p should be "near" the surface of the shape, although
// any p should produce some n.
int shape3_normal(const shape3 *s, const double p[3], double n[3]); // not implemented
int shape2_normal(const shape2 *s, const double p[2], double n[2]); // not implemented

// Determines if a shape intersects a given line segment defined by the point
// and vector. Returns the number of intersections (up to 2) in t. The values
// in t are the offsets along v, and are always in the range [0,1].
int shape3_segment_intersect(const shape3 *s,
	const double a[3], const double v[3], double t[2]); // not implemented
int shape2_segment_intersect(const shape2 *s,
	const double a[2], const double v[2], double t[2]); // not implemented

// Returns the fourier transform (real and imag part in ft) of the shape
// at the k-point 2*pi*f. There is no normalization factor to the Fourier
// integral.
int shape3_fourier_transform(const shape3 *s, const double f[3], double ft[2]);
int shape2_fourier_transform(const shape2 *s, const double f[2], double ft[2]);

#include <stdio.h>

// Output a 3D shape description to a POVRay block. The content string is output
// after the shape description to allow setting material properties.
int shape3_output_POVRay(const shape3 *s, FILE *fp, const char *content);
int aabb3_output_POVRay(const aabb3 *b, FILE *fp, const char *content);

#define SHAPE_OUTPUT_POSTSCRIPT_FILL        1 // apply a fill to the figure
#define SHAPE_OUTPUT_POSTSCRIPT_NOSTROKE    2 // do not stroke the figure
#define SHAPE_OUTPUT_POSTSCRIPT_NOTRANSLATE 4 // draw using only local coords

// Output a 2D shape description to PostScript commands. The options
// are defined above.
int shape2_output_postscript(const shape2 *s, FILE *fp, int opts);
int aabb2_output_postscript(const aabb2 *b, FILE *fp, int opts);


static const char Shape2_typename[] = "SharedNumerics.shape2";
static const char Shape3_typename[] = "SharedNumerics.shape3";


#endif // _SHAPES_H_
