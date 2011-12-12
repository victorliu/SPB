#ifndef _SHAPE_SET_H_
#define _SHAPE_SET_H_

#include "Shapes.h"

// The source file depends on BVH.c/h

typedef struct ShapeSet2_* ShapeSet2;
typedef struct ShapeSet3_* ShapeSet3;

ShapeSet2 ShapeSet2_new(const double lattice[4]);
ShapeSet3 ShapeSet3_new(const double lattice[6]);

void ShapeSet2_destroy(ShapeSet2 ss);
void ShapeSet3_destroy(ShapeSet3 ss);

int ShapeSet2_add(ShapeSet2 ss, const shape2 *s, int tag);
int ShapeSet3_add(ShapeSet3 ss, const shape3 *s, int tag);

// Returns item of largest tag
int ShapeSet2_query_pt(ShapeSet2 ss, const double p[2], const shape2 **s, int *tag);
int ShapeSet3_query_pt(ShapeSet3 ss, const double p[3], const shape3 **s, int *tag);

static const char ShapeSet2_typename[] = "SharedNumerics.ShapeSet2";
static const char ShapeSet3_typename[] = "SharedNumerics.ShapeSet3";

#endif // _SHAPE_SET_H_
