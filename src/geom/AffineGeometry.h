#ifndef _AFFINE_GEOMETRY_H_
#define _AFFINE_GEOMETRY_H_

typedef struct{ double v[2]; } vec2;
typedef struct{ double v[3]; } vec3;
typedef struct{ double v[4]; } vec4;
typedef struct{ double r[2]; } pt2;
typedef struct{ double r[3]; } pt3;

int vec2_normalize(vec2 *v);
int vec3_normalize(vec3 *v);
double vec2_dot(const vec2 *u, const vec2 *v);
double vec3_dot(const vec3 *u, const vec3 *v);
double vec2_cross(const vec2 *u, const vec2 *v);
void vec3_cross(const vec3 *u, const vec3 *v, vec3 *result);

static const char vec2_typename[] = "AffineGeometry.vec2";
static const char vec3_typename[] = "AffineGeometry.vec3";
static const char pt2_typename[] = "AffineGeometry.pt2";
static const char pt3_typename[] = "AffineGeometry.pt3";

// void AffineGeometry_init(lua_State *L);

#endif // _AFFINE_GEOMETRY_H_
