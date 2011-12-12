#include <stdlib.h>
#include <stdio.h>
#include "Shapes.h"
#include <math.h>
#include "lp_solve.h"

static double inline hypot22(double a, double b){
	return a*a+b*b;
}
static double inline hypot2v(const double v[2]){
	return hypot(v[0],v[1]);
}
static double inline hypot2v2(const double v[2]){
	return v[0]*v[0] + v[1]*v[1];
}
static double inline hypot3(double a, double b, double c){
	double w;
	double aa = fabs(a);
	double ab = fabs(b);
	double ac = fabs(c);
	if(aa > ab){
		if(aa > ac){
			w = aa;
		}else{
			w = ac;
		}
	}else{
		if(ab > ac){
			w = ab;
		}else{
			w = ac;
		}
	}
	if(0 == w){ return aa+ab+ac; }
	aa /= w; ab /= w; ac /= w;
	return w * sqrt(aa*aa + ab*ab + ac*ac);
}
static double inline hypot32(double a, double b, double c){
	return a*a+b*b+c*c;
}
static double inline hypot3v(const double v[3]){
	return hypot3(v[0],v[1],v[2]);
}
static double inline hypot3v2(const double v[3]){
	return v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
}

static double polygon_area(int n, const double *v){
	int p, q;
	double area = 0;
	for(p=n-1, q=0; q < n; p = q++){
		area += v[2*p+0]*v[2*q+1] - v[2*q+0]*v[2*p+1];
	}
	return 0.5*area;
}
static double mat3_det(const double m[9]){
	return
		m[0]*m[4]*m[8] - m[6]*m[4]*m[2] +
		m[3]*m[7]*m[2] + m[6]*m[1]*m[5] -
		m[3]*m[1]*m[8] - m[0]*m[7]*m[5];
}
static double mat3_invert(const double m[9], double mi[9]){
	double ret = mat3_det(m);
	double detinv = 1./ret;

	mi[0] = detinv * (m[4] * m[8] - m[7] * m[5]);
	mi[4] = detinv * (m[0] * m[8] - m[6] * m[2]);
	mi[8] = detinv * (m[4] * m[0] - m[1] * m[3]);
	mi[2] = detinv * (m[1] * m[5] - m[4] * m[2]);
	mi[1] = -detinv * (m[1] * m[8] - m[7] * m[2]);
	mi[5] = -detinv * (m[0] * m[5] - m[3] * m[2]);
	mi[6] = detinv * (m[3] * m[7] - m[4] * m[6]);
	mi[3] = -detinv * (m[3] * m[8] - m[5] * m[6]);
	mi[7] = -detinv * (m[0] * m[7] - m[1] * m[6]);
	return ret;
}
static void mat3_multv(const double m[9], const double v[3], double r[3]){
	r[0] = m[0]*v[0] + m[3]*v[1] + m[6]*v[2];
	r[1] = m[1]*v[0] + m[4]*v[1] + m[7]*v[2];
	r[2] = m[2]*v[0] + m[5]*v[1] + m[8]*v[2];
}
static void mat3_multTv(const double m[9], const double v[3], double r[3]){
	r[0] = m[0]*v[0] + m[1]*v[1] + m[2]*v[2];
	r[1] = m[3]*v[0] + m[4]*v[1] + m[5]*v[2];
	r[2] = m[6]*v[0] + m[7]*v[1] + m[8]*v[2];
}
static double mat2_invert(const double m[4], double mi[4]){
	double ret = m[0]*m[3]-m[1]*m[2];
	double detinv = 1./ret;

	mi[0] = detinv * (m[3]);
	mi[2] = detinv * (-m[2]);
	mi[1] = detinv * (-m[1]);
	mi[3] = detinv * (m[0]);
	return ret;
}
















void aabb3_union(aabb3 *b1, const aabb3 *b2){
	int i;
	for(i = 0; i < 3; ++i){
		double mn1 = b1->c[i] - b1->h[i];
		double mn2 = b2->c[i] - b2->h[i];
		double mx1 = b1->c[i] + b1->h[i];
		double mx2 = b2->c[i] + b2->h[i];
		if(mn2 < mn1){ mn1 = mn2; }
		if(mx2 > mx1){ mx1 = mx2; }
		b1->c[i] = 0.5*(mx1+mn1);
		b1->h[i] = 0.5*(mx1-mn1);
	}
}
void aabb2_union(aabb2 *b1, const aabb2 *b2){
	int i;
	for(i = 0; i < 3; ++i){
		double mn1 = b1->c[i] - b1->h[i];
		double mn2 = b2->c[i] - b2->h[i];
		double mx1 = b1->c[i] + b1->h[i];
		double mx2 = b2->c[i] + b2->h[i];
		if(mn2 < mn1){ mn1 = mn2; }
		if(mx2 > mx1){ mx1 = mx2; }
		b1->c[i] = 0.5*(mx1+mn1);
		b1->h[i] = 0.5*(mx1-mn1);
	}
}
void aabb3_union_pt(aabb3 *b1, const double p[3]){
	int i;
	for(i = 0; i < 3; ++i){
		double mn1 = b1->c[i] - b1->h[i];
		double mx1 = b1->c[i] + b1->h[i];
		if(p[i] < mn1){ mn1 = p[i]; }
		if(p[i] > mx1){ mx1 = p[i]; }
		b1->c[i] = 0.5*(mx1+mn1);
		b1->h[i] = 0.5*(mx1-mn1);
	}
}
void aabb2_union_pt(aabb2 *b1, const double p[2]){
	int i;
	for(i = 0; i < 2; ++i){
		double mn1 = b1->c[i] - b1->h[i];
		double mx1 = b1->c[i] + b1->h[i];
		if(p[i] < mn1){ mn1 = p[i]; }
		if(p[i] > mx1){ mx1 = p[i]; }
		b1->c[i] = 0.5*(mx1+mn1);
		b1->h[i] = 0.5*(mx1-mn1);
	}
}



int aabb3_contains(const aabb3 *b, const double p[3]){
	return
		(fabs(p[0]-b->c[0]) < b->h[0]) &&
		(fabs(p[1]-b->c[1]) < b->h[1]) &&
		(fabs(p[2]-b->c[2]) < b->h[2]);
}
int aabb2_contains(const aabb2 *b, const double p[2]){
	return
		(fabs(p[0]-b->c[0]) < b->h[0]) &&
		(fabs(p[1]-b->c[1]) < b->h[1]);
}

int aabb3_intersects(const aabb3 *a, const aabb3 *b){
	// separating axis test
	return
		(fabs(b->c[0] - a->c[0]) <= (a->h[0] + b->h[0])) &&
		(fabs(b->c[1] - a->c[1]) <= (a->h[1] + b->h[1])) &&
		(fabs(b->c[2] - a->c[2]) <= (a->h[2] + b->h[2]));
}

int aabb2_intersects(const aabb2 *a, const aabb2 *b){
	// separating axis test
	return
		(fabs(b->c[0] - a->c[0]) <= (a->h[0] + b->h[0])) &&
		(fabs(b->c[1] - a->c[1]) <= (a->h[1] + b->h[1]));
}

int shape3_contains(const shape3 *s, const double p[3]){
	const double po[3] = {p[0]-s->o[0], p[1]-s->o[1], p[2]-s->o[2]};
	switch(s->type){
	case SHAPE3_SPHERE:
		return hypot3(po[0], po[1], po[2]) <= s->sphere.r;
	case SHAPE3_ELLIPSOID:
		{
			double v[3];
			int i ,j;
			for(i = 0; i < 3; ++i){
				v[i] = 0;
				for(j = 0; j < 3; ++j){
					v[i] += s->ellipsoid.B[i+j*3]*po[j];
				}
			}
			return hypot3v(v) <= 1.;
		}
	case SHAPE3_BLOCK:
		{
			int i ,j;
			for(i = 0; i < 3; ++i){
				double v = 0;
				for(j = 0; j < 3; ++j){
					v += s->block.B[i+j*3]*po[j];
				}
				if(fabs(v) > 1){ return 0; }
			}
			return 1;
		}
	case SHAPE3_FRUSTUM:
		{
			double proj = (po[0]*s->frustum.Q[6] + po[1]*s->frustum.Q[7] + po[2]*s->frustum.Q[8]);
			double rz = proj/s->frustum.len;
			if(rz < 0. || rz > 1.){ return 0; }
			return hypot3(
				po[0] - proj*s->frustum.Q[6],
				po[1] - proj*s->frustum.Q[7],
				po[2] - proj*s->frustum.Q[8]
			) <= (1. + (s->frustum.r_tip_base-1.) * rz);
		}
	case SHAPE3_TETRAHEDRON:
		{
			const double org[3] = {0,0,0};
			double d0 = orient3d(org, s->tetrahedron.p[0], s->tetrahedron.p[1], s->tetrahedron.p[2]);
			double d1 = orient3d(po, s->tetrahedron.p[0], s->tetrahedron.p[1], s->tetrahedron.p[2]);
			double d2 = orient3d(org, po, s->tetrahedron.p[1], s->tetrahedron.p[2]);
			double d3 = orient3d(org, s->tetrahedron.p[0], po, s->tetrahedron.p[2]);
			double d4 = orient3d(org, s->tetrahedron.p[0], s->tetrahedron.p[1], po);
			return (d0 > 0 && d1 > 0 && d2 > 0 && d3 > 0 && d4 > 0)
				|| (d0 < 0 && d1 < 0 && d2 < 0 && d3 < 0 && d4 < 0);
		}
	case SHAPE3_CONVEX:
		{
			int i;
			for(i = 0; i < s->convex.np; ++i){
				double pp[3] = {
					po[0] - s->convex.p[3*i+0],
					po[1] - s->convex.p[3*i+1],
					po[2] - s->convex.p[3*i+2]};
				if(pp[0]*s->convex.n[3*i+0] +
				   pp[1]*s->convex.n[3*i+1] +
				   pp[2]*s->convex.n[3*i+2] > 0){ return 0; }
			}
			return 1;
		}
	default:
		return 0;
	}
}
int shape2_contains(const shape2 *s, const double p[2]){
	// offset vector from local origin
	const double po[2] = {p[0]-s->o[0], p[1]-s->o[1]};
	switch(s->type){
	case SHAPE2_CIRCLE:
		return hypot2v(po) <= s->circle.r;
	case SHAPE2_ELLIPSE:
		{
			double v[2];
			v[0] = s->ellipse.B[0]*po[0] + s->ellipse.B[2]*po[1];
			v[1] = s->ellipse.B[1]*po[0] + s->ellipse.B[3]*po[1];
			return hypot2v(v) <= 1.;
		}
	case SHAPE2_QUAD:
		{
			double v[2];
			v[0] = s->quad.B[0]*po[0] + s->quad.B[2]*po[1];
			v[1] = s->quad.B[1]*po[0] + s->quad.B[3]*po[1];
			return fabs(v[0]) <= 1. && fabs(v[1]) <= 1.;
		}
	case SHAPE2_POLYGON:
		{
			int i, j;
			int c = 0;
			for(i = 0, j = s->polygon.n-1; i < s->polygon.n; j = i++){
				double vix = s->polygon.p[2*i+0];
				double viy = s->polygon.p[2*i+1];
				double vjx = s->polygon.p[2*j+0];
				double vjy = s->polygon.p[2*j+1];
				if ( ((viy>po[1]) != (vjy>po[1]))
				&& (po[0] < (vjx-vix) * (po[1]-viy) / (vjy-viy) + vix) ){ c = !c; }
			}
			return c;
		}
	case SHAPE2_ARCSPLINE:
		{
			return 0;
		}
	default:
		return 0;
	}
}

int shape3_get_aabb(const shape3 *s, aabb3 *b){
	switch(s->type){
	case SHAPE3_SPHERE:
		b->c[0] = s->o[0]; b->c[1] = s->o[1]; b->c[2] = s->o[2];
		b->h[2] = b->h[1] = b->h[0] = s->sphere.r;
		return 0;
	case SHAPE3_BLOCK:
		{
			static const double sgn[8] = {
				1.,1.,
				1.,-1.,
				-1.,1.,
				-1.,-1.
			};
			int i;
			b->c[0] = s->o[0]; b->c[1] = s->o[1]; b->c[2] = s->o[2];
			b->h[2] = b->h[1] = b->h[0] = 0;
			
			for(i = 0; i < 4; ++i){
				const double u[3] = {
					fabs(s->block.A[3*0+0] + sgn[2*i+0]*s->block.A[3*1+0] + sgn[2*i+1]*s->block.A[3*2+0]),
					fabs(s->block.A[3*0+1] + sgn[2*i+0]*s->block.A[3*1+1] + sgn[2*i+1]*s->block.A[3*2+1]),
					fabs(s->block.A[3*0+2] + sgn[2*i+0]*s->block.A[3*1+2] + sgn[2*i+1]*s->block.A[3*2+2])};
				if(u[0] > b->h[0]){ b->h[0] = u[0]; }
				if(u[1] > b->h[1]){ b->h[1] = u[1]; }
				if(u[2] > b->h[2]){ b->h[2] = u[2]; }
			}
			return 0;
		}
	case SHAPE3_FRUSTUM:
		{
			// compute bounding boxes of end caps and union
			// The end caps are circles with radius R and normal n:
			//   (r-c).n = 0, (r-c).(r-c) = R^2, (r-c) = {x,y,z}
			// Solving for x,y,z and setting discriminant equal to zero gives
			//   x = R/|n| sqrt(ny^2+nz^2)
			//   y = R/|n| sqrt(nx^2+nz^2)
			//   z = R/|n| sqrt(nx^2+ny^2)
			// The centers of the boxes are just the centers of the circles.
			{
				double sbase = s->frustum.r_base/s->frustum.len;
				double stip= s->frustum.r_tip_base*sbase;
				aabb3 btip;
				b->c[0] = s->o[0]; b->c[1] = s->o[1]; b->c[2] = s->o[2];
				btip.c[0] = s->o[0]+s->frustum.Q[6];
				btip.c[1] = s->o[1]+s->frustum.Q[7];
				btip.c[2] = s->o[2]+s->frustum.Q[8];
				b->h[0] = sbase * hypot2(s->frustum.Q[6], s->frustum.Q[7]);
				b->h[1] = sbase * hypot2(s->frustum.Q[6], s->frustum.Q[8]);
				b->h[2] = sbase * hypot2(s->frustum.Q[7], s->frustum.Q[8]);
				btip.h[0] = stip * hypot2(s->frustum.Q[6], s->frustum.Q[7]);
				btip.h[1] = stip * hypot2(s->frustum.Q[6], s->frustum.Q[8]);
				btip.h[2] = stip * hypot2(s->frustum.Q[7], s->frustum.Q[8]);
				aabb3_union(b, &btip);
			}
			return 0;
		}
	case SHAPE3_ELLIPSOID:
		{
			// The ellipse is defined by the equation
			//   {x,y,z}.(P.{x,y,z}) == 1
			// The bounds are: max = sqrt(diag(inv(P)))
			// Since P = B^T B, and B = inv({e}),
			// inv(P) = inv(B)*inv(B^T) = {e}*{e}^T
			b->c[0] = s->o[0]; b->c[1] = s->o[1]; b->c[2] = s->o[2];
			b->h[0] = sqrt(s->ellipsoid.A[0]*s->ellipsoid.A[0] + s->ellipsoid.A[3]*s->ellipsoid.A[3] + s->ellipsoid.A[6]*s->ellipsoid.A[6]);
			b->h[1] = sqrt(s->ellipsoid.A[1]*s->ellipsoid.A[1] + s->ellipsoid.A[4]*s->ellipsoid.A[4] + s->ellipsoid.A[7]*s->ellipsoid.A[7]);
			b->h[2] = sqrt(s->ellipsoid.A[2]*s->ellipsoid.A[2] + s->ellipsoid.A[5]*s->ellipsoid.A[5] + s->ellipsoid.A[8]*s->ellipsoid.A[8]);
		}
		return 0;
	case SHAPE3_TETRAHEDRON:
		{
			int i, j;
			b->c[0] = s->o[0]; b->c[1] = s->o[1]; b->c[2] = s->o[2];
			b->h[2] = b->h[1] = b->h[0] = 0;
			for(j = 0; j < 3; ++j){
				double p[3];
				for(i = 0; i < 3; ++i){
					p[i] = s->o[i] + s->tetrahedron.p[3*j+i];
				}
				aabb3_union_pt(b, p);
			}
		}
		return 0;
	case SHAPE3_CONVEX:
		// given a set of hyperplanes, find the supporting hyperplanes with axis aligned normals
		// The convex region is defined as
		//   [  n0     n1    ... ]^T [ r ] <= 0
		//   [ -n0.p0 -n1.p1 ... ]   [ 1 ]
		// or we can write C^T x <= 0
		// We seek to find hyperplanes P = [ ei; -ei.x ], where ei is an axis aligned unit vector, such that
		//   min -ei.x
		//    x
		//   s.t. C^T x <= 0
		// Let N denote the matrix of normals, one per column, and P the column vector of normal.planepoint for each hyperplane
		// The problem can be written as
		//   min -ei.r
		//    r
		//   s.t. N^T r <= P
		// We introduce slack variables d, which correspond to the distance away from each hyperplane,
		//   min -ei.(r+ - r-)
		//   s.t. N^T (r+ - r-) + d = P
		//        r+, r-, d >= 0
		//        r = r+ - r-
		// To place this in standard form, let
		//   x^T = [ d r+ r- ]
		//   c^T = [ 0 -ei^T ei^T ]
		//   A = [ I N^T -N^T ]
		//   b = P
		// so that
		//   min c^T x
		//   s.t. A x = b, x >= 0
		// We must solve this problem 6 times, one for each side of the bounding box.
		// Also, the dimension of A is #planes x #planes+6
		{
			const int m = s->convex.np;
			const int n = m+6;
			int i, j, ret = 0;
			double *x = (double*)malloc(sizeof(double)*n);
			enum lp_solution_status status;
			int nsteps;
			lp_t lp;
			lp_init(&lp, m,n);
			// A and b never change
			for(j = 0; j < m; ++j){
				for(i = 0; i < m; ++i){
					lp.A[i+j*m] = (i==j);
				}
			}
			for(j = 0; j < 3; ++j){
				for(i = 0; i < m; ++i){
					lp.A[i+(m+j)*m] = s->convex.n[3*i+j];
				}
			}
			for(j = 0; j < 3; ++j){
				for(i = 0; i < m; ++i){
					lp.A[i+(m+3+j)*m] = -lp.A[i+(m+j)*m];
				}
			}
			for(i = 0; i < m; ++i){
				lp.b[i] = s->convex.n[3*i+0]*s->convex.p[3*i+0] +
				          s->convex.n[3*i+1]*s->convex.p[3*i+1] +
				          s->convex.n[3*i+2]*s->convex.p[3*i+2];
			}
			for(i = 0; i < m; ++i){
				lp.c[i] = 0;
			}
			for(i = 0; i < 3; ++i){
				double mn, mx;
				lp.c[m+5] = lp.c[m+4] = lp.c[m+3] = lp.c[m+2] = lp.c[m+1] = lp.c[m+0] = 0;
				
				lp.c[m+i] = 1;
				lp.c[m+3+i] = -1;
				if(0 != lp_solve(&lp, x, NULL, NULL, &status, &nsteps)){ ret = 2; goto lp_error; }
				mn = x[m+i] - x[m+3+i];
				
				lp.c[m+i] = -1;
				lp.c[m+3+i] = 1;
				if(0 != lp_solve(&lp, x, NULL, NULL, &status, &nsteps)){ ret = 2; goto lp_error; }
				mx = x[m+i] - x[m+3+i];
				b->c[i] = 0.5*(mx+mn);
				b->h[i] = 0.5*(mx-mn);
			}
lp_error:
			free(x);
			lp_destroy(&lp);
			return ret;
		}
	default:
		return 1;
	}
}
int shape2_get_aabb(const shape2 *s, aabb2 *b){
	switch(s->type){
	case SHAPE2_CIRCLE:
		b->c[0] = s->o[0]; b->c[1] = s->o[1];
		b->h[1] = b->h[0] = s->circle.r;
		return 0;
	case SHAPE2_ELLIPSE:
		{
			// The ellipse is defined by the equation
			//   {x,y}.(A.{x,y}) == 1
			// Expanding and setting the descriminant equal to zero gives the following simple solution.
			b->c[0] = s->o[0]; b->c[1] = s->o[1];
			b->h[0] = sqrt(s->ellipse.A[0]*s->ellipse.A[0] + s->ellipse.A[2]*s->ellipse.A[2]);
			b->h[1] = sqrt(s->ellipse.A[1]*s->ellipse.A[1] + s->ellipse.A[3]*s->ellipse.A[3]);
			return 0;
		}
	case SHAPE2_QUAD:
		{
			b->c[0] = s->o[0]; b->c[1] = s->o[1];
			const double u[2] = {
				fabs(s->quad.A[2*0+0] + s->quad.A[2*1+0]),
				fabs(s->quad.A[2*0+1] + s->quad.A[2*1+1])};
			const double v[2] = {
				fabs(s->quad.A[2*0+0] - s->quad.A[2*1+0]),
				fabs(s->quad.A[2*0+1] - s->quad.A[2*1+1])};
			b->h[0] = u[0] > v[0] ? u[0] : v[0];
			b->h[1] = u[1] > v[1] ? u[1] : v[1];
			return 0;
		}
	case SHAPE2_POLYGON:
		{
			int i;
			double mn[2] = {s->polygon.p[0],s->polygon.p[1]};
			double mx[2] = {s->polygon.p[0],s->polygon.p[1]};
			for(i = 1; i < s->polygon.n; ++i){
				if(s->polygon.p[2*i+0] < mn[0]){
					mn[0] = s->polygon.p[2*i+0];
				}
				if(s->polygon.p[2*i+1] < mn[1]){
					mn[1] = s->polygon.p[2*i+1];
				}
				if(s->polygon.p[2*i+0] > mx[0]){
					mx[0] = s->polygon.p[2*i+0];
				}
				if(s->polygon.p[2*i+1] > mx[1]){
					mx[1] = s->polygon.p[2*i+1];
				}
			}
			b->c[0] = 0.5*mn[0] + 0.5*mx[0];
			b->c[1] = 0.5*mn[1] + 0.5*mx[1];
			b->h[0] = mx[0] - mn[0];
			b->h[1] = mx[1] - mn[1];
			return 0;
		}
	case SHAPE2_ARCSPLINE:
		return 1;
	default:
		return 1;
	}
}

int shape3_intersects(const shape3 *s, const shape3 *t){
	return 0;
}
int shape2_intersects(const shape2 *s, const shape2 *t){
	const double r[2] = {t->o[0]-s->o[0], t->o[1]-s->o[1]};
	switch(s->type){
	case SHAPE2_CIRCLE:
		switch(t->type){
		case SHAPE2_CIRCLE:
			{
				double d = hypot2(r[0], r[1]);
				return d <= s->circle.r + t->circle.r;
			}
		default:
			return 0;
		}
	default:
		return 0;
	}
	return 0;
}

double shape3_simplex_overlap(const shape3 *s, const double t[4*3]);
double shape2_simplex_overlap(const shape2 *s, const double t[3*2]){
	return 0;
}

int shape3_aabb_overlap(const shape3 *s, const aabb3 *b);
int shape2_aabb_overlap(const shape2 *s, const aabb2 *b);

int shape3_normal(const shape3 *s, const double p[3], double n[3]){
	const double po[3] = {p[0]-s->o[0], p[1]-s->o[1], p[2]-s->o[2]};
	switch(s->type){
	case SHAPE3_SPHERE:
		{
			n[0] = po[0];
			n[1] = po[1];
			n[2] = po[2];
		}
		break;
	case SHAPE3_ELLIPSOID:
		{
			// The ellipsoid is defined by
			//   (B.po)^2 = 1
			// The gradient is
			//   B^T.B.po
			double Bpo[3];
			int i,j;
			for(i = 0; i < 3; ++i){
				Bpo[3] = 0;
				for(j = 0; j < 3; ++j){
					Bpo[i] += s->ellipsoid.B[2*j+i]*po[j];
				}
			}
			for(i = 0; i < 3; ++i){
				n[i] = 0;
				for(j = 0; j < 3; ++j){
					n[i] += s->ellipsoid.B[2*i+j]*po[j];
				}
			}
		}
		break;
	case SHAPE3_BLOCK:
		{
			// The cube is defined by
			//   norm(B.po, infty) = 1
			// where the columns of B are s->block.e. Let Bpo = B.po.
			//   norm(Bpo,infty) = max(abs(Bpo[0]),abs(Bpo[1]),abs(Bpo[2]))
			// We compute the gradient:
			//   Let I = argmax(max(abs(Bpo[i])))
			// Then
			//   grad = sgn(Bpo[I]) {I-th row of B}
			double Bpo[3];
			int i,j;
			for(i = 0; i < 3; ++i){
				Bpo[3] = 0;
				for(j = 0; j < 3; ++j){
					Bpo[i] += s->block.B[i+j*3]*po[j];
				}
			}
			const double ab[3] = {fabs(Bpo[0]),fabs(Bpo[1]),fabs(Bpo[2])};
			int I = 0;
			if(ab[0] > ab[1]){
				if(ab[2] > ab[0]){
					I = 2;
				}else{
					// I = 0;
				}
			}else{
				if(ab[2] > ab[1]){
					I = 2;
				}else{
					I = 1;
				}
			}
			double sgn = (Bpo[I] > 0 ? 1. : -1);
			n[0] = sgn * s->block.B[I+0];
			n[1] = sgn * s->block.B[I+3];
			n[2] = sgn * s->block.B[I+6];
		}
		break;
	case SHAPE3_FRUSTUM:
		// Let Q = {v1, v2, axis} be an orthogonal set.
		// Let the un-transformed point r = inv(Q).po;
		// If Q = U.diag{S} where U is orthogonal, then inv(Q) = inv(diag(S)).U^T
		// The frustum is defined by
		//   max(
		//     2*norm_inf(r.z - 0.5),
		//     norm_2(r.xy) / (1 + (r_tip_base-1)*r.z)
		//   ) <= 1
		// Let denom = (1 + (r_tip_base-1)*r.z)
		// To compute the gradient, determine which of the two arguments
		// to max is active, then,
		//   grad = Q.{
		//     sgn(r.z - 0.5) * <0,0,1>
		//     2 r.xy / (|r.xy| denom) + (1-r_tip_base) norm_2(r.xy) / denom
		//   }
		{
			double r[3]; mat3_multTv(s->frustum.Q, po, r);
			r[0] /= s->frustum.r_base;
			r[1] /= s->frustum.r_base;
			r[2] /= s->frustum.len;
			const double denom = 1. + (s->frustum.r_tip_base-1.)*r[2];
			const double nrxy = hypot2v(r);
			if(2*fabs(r[2]-0.5) > nrxy/denom){
				double sgn = (r[2] > 0.5 ? 1. : -1.);
				n[0] = sgn*s->frustum.Q[6];
				n[1] = sgn*s->frustum.Q[7];
				n[2] = sgn*s->frustum.Q[8];
			}else{
				r[0] *= 2./nrxy;
				r[1] *= 2./nrxy;
				r[2] = 1.-s->frustum.r_tip_base * nrxy;
				mat3_multv(s->frustum.Q, r, n);
			}
		}
		break;
	case SHAPE3_TETRAHEDRON:
		{
			// See convex below.
			double N[3*4];
			int i;
		}
		break;
	case SHAPE3_CONVEX:
		// Arrange all the points and normals as
		//   P = [ n0    n1    ... ] in R^{3 x N}
		//   d = [ n0.p0 n1.p1 ... ]^T column vector
		// The polyhedron is defined by
		//   P^T.po <= d
		// or
		//   max(P^T.po - d) <= 0
		// The gradient is then
		//   grad = n[argmax(P^T.po-d)];
		break;
	default:
		break;
	}
	double ilen = 1./hypot3v(n);
	n[0] *= ilen;
	n[1] *= ilen;
	n[2] *= ilen;
	return 0;
}

int shape2_normal(const shape2 *s, const double p[2], double n[2]);

int shape3_segment_intersect(const shape3 *s, const double p[3], const double v[3], double t[2]){
	const double po[3] = {p[0] - s->o[0], p[1] - s->o[1], p[2] - s->o[2]};
	int nret = 0;
	switch(s->type){
	case SHAPE3_SPHERE:
		{
			// The ellipsoid surface is defined by
			//   norm_w(r)^2 = 1
			// Parameterize the segment by 
			//   r = po + t v
			//   norm_w(po + t v)^2 = 1
			//   po.po + 2 t v.po + t^2 v.v = 1
			const double a = vec3_dot(v,v);
			const double b = 2*vec3_dot(po,v) / a;
			const double c = (vec3_dot(po,po)-1.) / a;
			// Solve t^2 + 2 b t + c == 0
			double disc = b*b-c;
			if(disc < 0){ nret = 0; }
			else if(0 == disc){ t[0] = -b; nret = 1; }
			else{
				disc = sqrt(disc);
				t[0] = -b - disc;
				t[1] = -b + disc;
				nret = 2;
			}
		}
		break;
	case SHAPE3_ELLIPSOID:
		{
			// The ellipsoid surface is defined by
			//   norm_w(B.r)^2 = 1
			// Parameterize the segment by 
			//   r = po + t v
			//   norm_w(B.(po + t v))^2 = 1
			// Let Bpo = B.po, Bv = B.v
			//   Bpo.Bpo + 2 t Bv.Bpo + t^2 Bv.Bv = 1
			double Bpo[3]; mat3_multv(s->ellipsoid.B, po, Bpo);
			double Bv[3];  mat3_multv(s->ellipsoid.B, v , Bv );
			const double a = vec3_dot(Bv,Bv);
			const double b = 2*vec3_dot(Bpo,Bv) / a;
			const double c = (vec3_dot(Bpo,Bpo)-1.) / a;
			// Solve t^2 + 2 b t + c == 0
			double disc = b*b-c;
			if(disc < 0){ nret = 0; }
			else if(0 == disc){ t[0] = -b; nret = 1; }
			else{
				disc = sqrt(disc);
				t[0] = -b - disc;
				t[1] = -b + disc;
				nret = 2;
			}
		}
		break;
	case SHAPE3_BLOCK:
		{
			// The block is defined by
			//   norm_inf(B.po + t B.v) = 1
			double Bpo[3]; mat3_multv(s->ellipsoid.B, po, Bpo);
			double Bv[3];  mat3_multv(s->ellipsoid.B, v , Bv );
			t[0] = -1.; t[1] = 2.;
			int i;
			for(i = 0; i < 3; ++i){
				double s;
				if(0 == Bv[i]){ continue; }
				s = ( 1 - Bpo[i])/Bv[i];
				if(s > 0){
					// todo finish me
				}
				if(t[0] < s && s < t[1]){ t[1] = s; }
				s = (-1 - Bpo[i])/Bv[i];
				if(s > t[0]){ t[0] = s; }
			}
			nret = 2; // have the cleanup code fix things
		}
		break;
	case SHAPE3_FRUSTUM:
		{
			// Let Q = {v1, v2, axis} be an orthogonal set.
			// Let the un-transformed point r = inv(Q).po;
			// If Q = U.diag{S} where U is orthogonal, then inv(Q) = inv(diag(S)).U^T
			// The frustum is defined by
			//   max(
			//     2*norm_inf(r.z - 0.5),
			//     norm_2(r.rho)^2 / (r_base*(1-r.z) + r_tip*r.z)
			//   ) = 1
			double r[3]; mat3_multTv(s->frustum.Q, po, r);
			r[0] /= s->frustum.r_base;
			r[1] /= s->frustum.r_base;
			r[2] /= s->frustum.len;
		}
		break;
	}
	// The following conditions should handle NaNs as well.
	if(nret == 1){
		if(0 <= t[0] && t[0] <= 1.){
			return 1;
		}else{
			return 0;
		}
	}else if(nret == 2){
		if(0 <= t[0] && t[0] <= 1.){
			if(0 <= t[1] && t[1] <= 1.){
				if(t[0] == t[1]){
					return 1;
				}else{
					return 2;
				}
			}else{
				return 1;
			}
		}else{
			if(0 <= t[1] && t[1] <= 1.){
				t[0] = t[1];
				return 1;
			}else{
				return 0;
			}
		}
	}
	return 0;
}
int shape2_segment_intersect(const shape2 *s, const double p[2], const double v[2], double t[2]){
}




int shape3_fourier_transform(const shape3 *s, const double f[3], double ft[2]){
	// Fourier transform of a sphere is
	//   len(f)^{-3/2} J_{3/2}(2 pi len(f))
	//   where J_{3/2}(r) = sqrt(2/(pi r)) [sin(r)/r - cos(r)]
}
int shape2_fourier_transform(const shape2 *s, const double f[2], double ft[2]){
	// Fourier transform of a disc is
	//   J_1(2 pi len(f))/len(f)
	// For k = 2 pi f != 0, FT of a polygon is
	//   S(k) = i/|k|^2 * Sum_{i=0,n-1} z.((v_{i+1}-v_{i}) x k) j0(k.(v_{i+1}-v_{i})/2) e^{ik.(v_{i+1}+v_{i})/2}
	// where j0(x) = sin(x)/x
}



#define FFMT "%f"

int shape3_output_POVRay(const shape3 *s, FILE *fp, const char *content){
	switch(s->type){
	case SHAPE3_SPHERE:
		fprintf(fp, "sphere{ <" FFMT "," FFMT "," FFMT ">," FFMT "\n",
			s->o[0], s->o[1], s->o[2], s->sphere.r);
		break;
	case SHAPE3_BLOCK:
		fprintf(fp, "box{ <-1,-1,-1>,<1,1,1>\n"
			"matrix <" FFMT "," FFMT "," FFMT ",\n"
			FFMT "," FFMT "," FFMT ",\n"
			FFMT "," FFMT "," FFMT ",\n"
			FFMT "," FFMT "," FFMT ">\n",
			s->block.A[3*0+0], s->block.A[3*1+0], s->block.A[3*2+0],
			s->block.A[3*0+1], s->block.A[3*1+1], s->block.A[3*2+1],
			s->block.A[3*0+2], s->block.A[3*1+2], s->block.A[3*2+2],
			s->o[0], s->o[1], s->o[2]
			);
		break;
	case SHAPE3_FRUSTUM:
		fprintf(fp, "cone{ <" FFMT "," FFMT "," FFMT ">," FFMT ",\n"
			"<" FFMT "," FFMT "," FFMT ">," FFMT "\n",
			s->o[0], s->o[1], s->o[2], s->frustum.r_base,
			s->o[0] + s->frustum.len * s->frustum.Q[3*2+0],
			s->o[1] + s->frustum.len * s->frustum.Q[3*2+1],
			s->o[2] + s->frustum.len * s->frustum.Q[3*2+2], s->frustum.r_tip_base*s->frustum.r_base
			);
		break;
	case SHAPE3_ELLIPSOID:
		{
			fprintf(fp, "sphere{ <0,0,0>," FFMT "\n"
				"matrix <"
				FFMT "," FFMT "," FFMT ",\n"
				FFMT "," FFMT "," FFMT ",\n"
				FFMT "," FFMT "," FFMT ",\n"
				FFMT "," FFMT "," FFMT ">\n",
				s->ellipsoid.A[3*0+0],s->ellipsoid.A[3*1+0],s->ellipsoid.A[3*2+0],
				s->ellipsoid.A[3*0+1],s->ellipsoid.A[3*1+1],s->ellipsoid.A[3*2+1],
				s->ellipsoid.A[3*0+2],s->ellipsoid.A[3*1+2],s->ellipsoid.A[3*2+2],
				s->o[0], s->o[1], s->o[2]
				);
		}
		break;
	case SHAPE3_TETRAHEDRON:
		fprintf(fp, "triangle{ <" FFMT "," FFMT "," FFMT ">,"
			"<" FFMT "," FFMT "," FFMT ">,"
			"<" FFMT "," FFMT "," FFMT ">\n",
			s->o[0], s->o[1], s->o[2],
			s->o[0] + s->tetrahedron.p[3*0+0],
			s->o[1] + s->tetrahedron.p[3*0+1],
			s->o[2] + s->tetrahedron.p[3*0+2],
			s->o[0] + s->tetrahedron.p[3*1+0],
			s->o[1] + s->tetrahedron.p[3*1+1],
			s->o[2] + s->tetrahedron.p[3*1+2]);
		fprintf(fp, "triangle{ <" FFMT "," FFMT "," FFMT ">,"
			"<" FFMT "," FFMT "," FFMT ">,"
			"<" FFMT "," FFMT "," FFMT ">\n",
			s->o[0], s->o[1], s->o[2],
			s->o[0] + s->tetrahedron.p[3*1+0],
			s->o[1] + s->tetrahedron.p[3*1+1],
			s->o[2] + s->tetrahedron.p[3*1+2],
			s->o[0] + s->tetrahedron.p[3*2+0],
			s->o[1] + s->tetrahedron.p[3*2+1],
			s->o[2] + s->tetrahedron.p[3*2+2]);
		fprintf(fp, "triangle{ <" FFMT "," FFMT "," FFMT ">,"
			"<" FFMT "," FFMT "," FFMT ">,"
			"<" FFMT "," FFMT "," FFMT ">\n",
			s->o[0], s->o[1], s->o[2],
			s->o[0] + s->tetrahedron.p[3*2+0],
			s->o[1] + s->tetrahedron.p[3*2+1],
			s->o[2] + s->tetrahedron.p[3*2+2],
			s->o[0] + s->tetrahedron.p[3*0+0],
			s->o[1] + s->tetrahedron.p[3*0+1],
			s->o[2] + s->tetrahedron.p[3*0+2]);
		fprintf(fp, "triangle{ <" FFMT "," FFMT "," FFMT ">,"
			"<" FFMT "," FFMT "," FFMT ">,"
			"<" FFMT "," FFMT "," FFMT ">\n",
			s->o[0] + s->tetrahedron.p[3*0+0],
			s->o[1] + s->tetrahedron.p[3*0+1],
			s->o[2] + s->tetrahedron.p[3*0+2],
			s->o[0] + s->tetrahedron.p[3*2+0],
			s->o[1] + s->tetrahedron.p[3*2+1],
			s->o[2] + s->tetrahedron.p[3*2+2],
			s->o[0] + s->tetrahedron.p[3*1+0],
			s->o[1] + s->tetrahedron.p[3*1+1],
			s->o[2] + s->tetrahedron.p[3*1+2]);
		break;
	case SHAPE3_CONVEX:
		{
			int i;
			fprintf(fp, "intersection{\n");
			for(i = 0; i < s->convex.np; ++i){
				double inv = 1./hypot3v(&s->convex.n[3*i]);
				fprintf(fp, "plane{<" FFMT "," FFMT "," FFMT ">," FFMT " }\n",
					s->convex.n[3-i+0], s->convex.n[3*i+1], s->convex.n[3*i+2],
					inv * (s->convex.n[3*i+0]*s->convex.p[3*i+0] +
					s->convex.n[3*i+1]*s->convex.p[3*i+1] +
					s->convex.n[3*i+2]*s->convex.p[3*i+2]));
			}
		}
		break;
	default:
		break;
	}
	fprintf(fp, "%s\n}\n", content);
	return 0;
}
int aabb3_output_POVRay(const aabb3 *b, FILE *fp, const char *content){
	fprintf(fp, "box{ <" FFMT "," FFMT "," FFMT ">,<" FFMT "," FFMT "," FFMT ">\n",
		b->c[0] - b->h[0], b->c[1] - b->h[1], b->c[2] - b->h[2],
		b->c[0] + b->h[0], b->c[1] + b->h[1], b->c[2] + b->h[2]
		);
	fprintf(fp, "%s\n}\n", content);
	return 0;
}

int shape2_output_postscript(const shape2 *s, FILE *fp, int opts){
	if(opts & SHAPE_OUTPUT_POSTSCRIPT_NOTRANSLATE){
	}else{
		fprintf(fp, "gsave " FFMT " " FFMT " translate\n", s->o[0], s->o[1]);
	}
	switch(s->type){
	case SHAPE2_CIRCLE:
		fprintf(fp, "0 0 " FFMT " 0 360 arc\n", s->circle.r);
		break;
	case SHAPE2_ELLIPSE:
		{
			double s00 = sqrt(s->ellipse.B[2*0+0]);
			double s01 = s->ellipse.B[2*0+1]/s00;
			// todo fixme
			fprintf(fp, "[" FFMT " " FFMT " " FFMT " " FFMT " 0 0] concat\n",
				s->ellipse.A[2*0+0], s->ellipse.A[2*1+0],
				s->ellipse.A[2*0+1], s->ellipse.A[2*1+1]
				);
			fprintf(fp, "0 0 1 0 360 arc\n");
		}
		break;
	case SHAPE2_QUAD:
		fprintf(fp, FFMT " " FFMT " moveto "
			FFMT " " FFMT " lineto "
			FFMT " " FFMT " lineto "
			FFMT " " FFMT " lineto closepath\n",
			-s->quad.A[2*0+0]-s->quad.A[2*1+0],
			-s->quad.A[2*0+1]-s->quad.A[2*1+1],
			+s->quad.A[2*0+0]-s->quad.A[2*1+0],
			+s->quad.A[2*0+1]-s->quad.A[2*1+1],
			+s->quad.A[2*0+0]+s->quad.A[2*1+0],
			+s->quad.A[2*0+1]+s->quad.A[2*1+1],
			-s->quad.A[2*0+0]+s->quad.A[2*1+0],
			-s->quad.A[2*0+1]+s->quad.A[2*1+1]);
		break;
	case SHAPE2_POLYGON:
		{
			int i;
			fprintf(fp, FFMT " " FFMT " moveto ", s->polygon.p[2*0+0], s->polygon.p[2*0+1]);
			for(i = 1; i < s->polygon.n; ++i){
				fprintf(fp, FFMT " " FFMT " lineto ", s->polygon.p[2*i+0], s->polygon.p[2*i+1]);
			}
			fprintf(fp, "closepath\n");
		}
		break;
	case SHAPE2_ARCSPLINE:
		break;
	default:
		break;
	}
	if((SHAPE_OUTPUT_POSTSCRIPT_FILL & opts) && !(SHAPE_OUTPUT_POSTSCRIPT_NOSTROKE & opts)){
		fprintf(fp, "gsave fill grestore\n");
	}else if(SHAPE_OUTPUT_POSTSCRIPT_FILL & opts){
		fprintf(fp, "fill\n");
	}else if(!(SHAPE_OUTPUT_POSTSCRIPT_NOSTROKE & opts)){
		fprintf(fp, "stroke\n");
	}
	
	if(opts & SHAPE_OUTPUT_POSTSCRIPT_NOTRANSLATE){
	}else{
		fprintf(fp, "grestore\n");
	}

	return 0;
}

int aabb2_output_postscript(const aabb2 *s, FILE *fp, int opts){
	if(opts & SHAPE_OUTPUT_POSTSCRIPT_NOTRANSLATE){
	}else{
		fprintf(fp, "gsave " FFMT " " FFMT " translate\n", s->c[0], s->c[1]);
	}
	fprintf(fp, FFMT " " FFMT " moveto "
		FFMT " " FFMT " lineto "
		FFMT " " FFMT " lineto "
		FFMT " " FFMT " lineto closepath\n",
		-s->h[0], -s->h[1],
		+s->h[0], -s->h[1],
		+s->h[0], +s->h[1],
		-s->h[0], +s->h[1]);
		
	if((SHAPE_OUTPUT_POSTSCRIPT_FILL & opts) && !(SHAPE_OUTPUT_POSTSCRIPT_NOSTROKE & opts)){
		fprintf(fp, "gsave fill grestore\n");
	}else if(SHAPE_OUTPUT_POSTSCRIPT_FILL & opts){
		fprintf(fp, "fill\n");
	}else if(!(SHAPE_OUTPUT_POSTSCRIPT_NOSTROKE & opts)){
		fprintf(fp, "stroke\n");
	}
	
	if(opts & SHAPE_OUTPUT_POSTSCRIPT_NOTRANSLATE){
	}else{
		fprintf(fp, "grestore\n");
	}

	return 0;
}





#ifdef HAVE_LUA


#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
#include <LuaCommon.h>


static shape2* make_shape2(lua_State *L, shape2_type type){
	shape2 *s = (shape2*)lua_newuserdata(L, sizeof(shape2));
	luaL_getmetatable(L, Shape2_typename);
	lua_setmetatable(L, -2);
	s->type = type;
	return s;
}
static shape3* make_shape3(lua_State *L, shape3_type type){
	shape3 *s = (shape3*)lua_newuserdata(L, sizeof(shape3));
	luaL_getmetatable(L, Shape3_typename);
	lua_setmetatable(L, -2);
	s->type = type;
	return s;
}

static int Shapes2D_circle(lua_State *L){
	shape2 *s = make_shape2(L, SHAPE2_CIRCLE);
	LP_named_arg argtab[] = {
		{"Center", LP_ARG_VEC3, 1, &(s->o[0])},
		{"Radius", LP_ARG_REAL_POS, 1, &(s->circle.r)},
		LP_ARG_END
	};
	if(0 != LP_get_args(L, 1, argtab)){
		LP_error(L, "Error in Shapes2D.circle");
	}
	return 1;
}

static int Shapes2D_ellipse(lua_State *L){
	shape2 *s = make_shape2(L, SHAPE2_ELLIPSE);
	int i;
	LP_named_arg argtab[] = {
		{"Center", LP_ARG_VEC3, 1, &(s->o[0])},
		{"Axis1", LP_ARG_VEC3, 1, &(s->ellipse.A[2*0+0])},
		{"Axis2", LP_ARG_VEC3, 1, &(s->ellipse.A[2*1+0])},
		LP_ARG_END
	};
	if(0 != LP_get_args(L, 1, argtab)){
		LP_error(L, "Error in Shapes2D.ellipse");
	}
	mat2_invert(s->ellipse.A, s->ellipse.B);
	return 1;
}

static int Shapes2D_quad(lua_State *L){
	shape2 *s = make_shape2(L, SHAPE2_QUAD);
	LP_named_arg argtab[] = {
		{"Center", LP_ARG_VEC3, 1, &(s->o[0])},
		{"Axis1", LP_ARG_VEC3, 1, &(s->quad.A[2*0+0])},
		{"Axis2", LP_ARG_VEC3, 1, &(s->quad.A[2*1+0])},
		LP_ARG_END
	};
	if(0 != LP_get_args(L, 1, argtab)){
		LP_error(L, "Error in Shapes2D.quad");
	}
	mat2_invert(s->quad.A, s->quad.B);
	return 1;
}

static int Shape2__gc(lua_State *L){
	shape2 *s = (shape2*)luaL_checkudata(L, 1, Shape2_typename);
	if(SHAPE2_POLYGON == s->type){
	}else if(SHAPE2_ARCSPLINE == s->type){
	}
	return 0;
}


static int Shapes3D_sphere(lua_State *L){
	shape3 *s = make_shape3(L, SHAPE3_SPHERE);
	LP_named_arg argtab[] = {
		{"Center", LP_ARG_VEC3, 1, &(s->o[0])},
		{"Radius", LP_ARG_REAL_POS, 1, &(s->sphere.r)},
		LP_ARG_END
	};
	if(0 != LP_get_args(L, 1, argtab)){
		LP_error(L, "Error in Shapes3D.sphere");
	}
	return 1;
}

static int Shapes3D_ellipsoid(lua_State *L){
	shape3 *s = make_shape3(L, SHAPE3_ELLIPSOID);
	int i;
	LP_named_arg argtab[] = {
		{"Center", LP_ARG_VEC3, 1, &(s->o[0])},
		{"Axis1", LP_ARG_VEC3, 1, &(s->ellipsoid.A[3*0+0])},
		{"Axis2", LP_ARG_VEC3, 1, &(s->ellipsoid.A[3*1+0])},
		{"Axis3", LP_ARG_VEC3, 1, &(s->ellipsoid.A[3*2+0])},
		LP_ARG_END
	};
	if(0 != LP_get_args(L, 1, argtab)){
		LP_error(L, "Error in Shapes3D.ellipse");
	}
	mat3_invert(s->ellipsoid.A, s->ellipsoid.B);
	return 1;
}

static int Shapes3D_block(lua_State *L){
	shape3 *s = make_shape3(L, SHAPE3_BLOCK);
	LP_named_arg argtab[] = {
		{"origin", LP_ARG_VEC3, 1, &(s->o[0])},
		{"axis1", LP_ARG_VEC3, 1, &(s->block.A[3*0+0])},
		{"axis2", LP_ARG_VEC3, 1, &(s->block.A[3*1+0])},
		{"axis3", LP_ARG_VEC3, 1, &(s->block.A[3*2+0])},
		LP_ARG_END
	};
	if(0 != LP_get_args(L, 1, argtab)){
		LP_error(L, "Error in Shapes3D.block");
	}
	mat3_invert(s->block.A, s->block.B);
	return 1;
}

static int Shapes3D_frustum(lua_State *L){
	shape3 *s = make_shape3(L, SHAPE3_FRUSTUM);
	double axis[3], r_tip;
	LP_named_arg argtab[] = {
		{"Base", LP_ARG_VEC3, 1, &(s->o[0])},
		{"Axis", LP_ARG_VEC3, 1, &(axis[0])},
		{"BaseRadius", LP_ARG_REAL_POS, 1, &(s->frustum.r_base)},
		{"TipRadius", LP_ARG_REAL_POS, 1, &(r_tip)},
		LP_ARG_END
	};
	if(0 != LP_get_args(L, 1, argtab)){
		LP_error(L, "Error in Shapes3D.frustum");
	}
	// todo: finish this
	return 1;
}

static int Shapes3D_tetrahedron(lua_State *L){
	shape3 *s = make_shape3(L, SHAPE3_FRUSTUM);
	LP_named_arg argtab[] = {
		{"v1", LP_ARG_VEC3, 1, &(s->o[0])},
		{"v2", LP_ARG_VEC3, 1, &(s->tetrahedron.p[3*0+0])},
		{"v3", LP_ARG_VEC3, 1, &(s->tetrahedron.p[3*1+0])},
		{"v4", LP_ARG_VEC3, 1, &(s->tetrahedron.p[3*2+0])},
		LP_ARG_END
	};
	if(0 != LP_get_args(L, 1, argtab)){
		LP_error(L, "Error in Shapes3D.tetrahedron");
	}
	return 1;
}

static int Shape3__gc(lua_State *L){
	shape3 *s = (shape3*)luaL_checkudata(L, 1, Shape3_typename);
	if(SHAPE3_CONVEX == s->type){
	}
	return 0;
}

static int Shape2_Output(lua_State *L){
	shape2 *s = (shape2*)luaL_checkudata(L, 1, Shape2_typename);
	luaL_argcheck(L, s != NULL, 1, "Shape2.shape:Output: object expected.");
	
	shape2_output_postscript(s, stdout, 0);
	
	return 0;
}

static int Shape3_Output(lua_State *L){
	shape3 *s = (shape3*)luaL_checkudata(L, 1, Shape3_typename);
	luaL_argcheck(L, s != NULL, 1, "Shape2.shape:Output: object expected.");
	
	shape3_output_POVRay(s, stdout, "");
	
	return 0;
}

void Shapes_init(lua_State *L){
	static const struct luaL_reg Shapes2D_lib[] = {
		{"circle"   , Shapes2D_circle},
		{"ellipse"  , Shapes2D_ellipse},
		{"quad"     , Shapes2D_quad},
		//{"polygon"  , Shapes2D_polygon},
		//{"arcspline", Shapes2D_arcspline},
		{NULL,NULL}
	};
	static const struct luaL_Reg Shape2_obj[] = {
		{"Output", Shape2_Output},
		{NULL, NULL}
	};
	static const struct luaL_reg Shapes3D_lib[] = {
		{"sphere"     , Shapes3D_sphere},
		{"ellipsoid"  , Shapes3D_ellipsoid},
		{"block"      , Shapes3D_block},
		{"frustum"    , Shapes3D_frustum},
		{"tetrahedron", Shapes3D_tetrahedron},
		//{"convex"     , Shapes3D_convex},
		{NULL,NULL}
	};
	static const struct luaL_Reg Shape3_obj[] = {
		{"Output", Shape3_Output},
		{NULL, NULL}
	};
	luaL_register(L, "Shapes2D", Shapes2D_lib);
	luaL_register(L, "Shapes3D", Shapes3D_lib);
	
	
	
	luaL_newmetatable(L, Shape2_typename);
	lua_pushvalue(L, -1);
	lua_setfield(L, -2, "__index");
	luaL_register(L, NULL, Shape2_obj);
	lua_pushstring(L, "__gc");
	lua_pushcfunction(L, Shape2__gc);
	lua_settable(L, -3);
	
	luaL_newmetatable(L, Shape3_typename);
	lua_pushvalue(L, -1);
	lua_setfield(L, -2, "__index");
	luaL_register(L, NULL, Shape3_obj);
	lua_pushstring(L, "__gc");
	lua_pushcfunction(L, Shape3__gc);
	lua_settable(L, -3);

}


#endif // HAVE_LUA
