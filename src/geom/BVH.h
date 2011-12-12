#ifndef _BVH_H_
#define _BVH_H_

// Bounding volume hierarchy structures
//
// The underlying implementation of the BVH is an R-tree. We assume
// that the structures are _static_; that is, they are created and
// never modified. This assumption allows us to generate more
// efficient trees by bulk-loading the trees up front.
//  The bulk loading method we use is Sort-Tile-Recursive (STR) due
// to its simplicity and reasonable effectiveness.

typedef struct BVH2_* BVH2;
typedef struct BVH3_* BVH3;

// b should be a box given by min and max coordinates
// In 2D, b[4] is {x-min, x-max, y-min, y-max}
// In 3D, b[6] is {x-min, x-max, y-min, y-max, z-min, z-max}
BVH2 BVH2_new(int n, int (*shape_iterator)(double b[4], int *tag, void *data), void *data);
BVH3 BVH3_new(int n, int (*shape_iterator)(double b[6], int *tag, void *data), void *data);

void BVH2_destroy(BVH2 bvh);
void BVH3_destroy(BVH3 bvh);

// p is the query point
// In 2D, p[2] is {x,y}. In 3D, p[3] is {x,y,z}
// The query function is passed boxes which contain the point p, along with the tag.
// The function should return nonzero to continue the query.
// A zero return value terminates the query.
int BVH2_query_pt(BVH2 bvh, const double p[2], int (*query_func)(int tag, const double b[4], void *data), void *data);
int BVH3_query_pt(BVH3 bvh, const double p[3], int (*query_func)(int tag, const double b[6], void *data), void *data);

int BVH2_query_box(BVH2 bvh, const double b[4], int (*query_func)(int tag, const double b[4], void *data), void *data);
int BVH3_query_box(BVH3 bvh, const double b[6], int (*query_func)(int tag, const double b[6], void *data), void *data);

int BVH2_traverse(BVH2 bvh, int (*func)(int tag, const double b[4], int internal, void *data), void *data);
int BVH3_traverse(BVH3 bvh, int (*func)(int tag, const double b[6], int internal, void *data), void *data);

#endif // _BVH_H_
