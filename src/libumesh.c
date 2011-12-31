#include "libumesh.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>

int LibUMesh2_Create(const double u[2], const double v[2], UMesh2 *mesh){
	if(NULL == u){ return -1; }
	if(NULL == v){ return -2; }
	if(NULL == mesh){ return -3; }
	
	const double uv = u[0]*v[0] + u[1]*v[1];
	const double uxv = u[0]*v[1] + u[1]*v[0];
	double len[3] = {
		hypot(u[0], u[1]),
		hypot(v[0], v[1]),
		0
	};
	
	// Error: degenerate lattice
	if(0 == len[0]){ return 1; }
	if(0 == len[1]){ return 2; }
	if(0 == uxv ){ return 3; }
	
	mesh->Lr[0] = u[0];
	mesh->Lr[1] = u[1];
	mesh->Lr[2] = v[0];
	mesh->Lr[3] = v[1];
	double *Lr = mesh->Lr;
	
	// u edge
	mesh->inc01[4*0+2*0+0] = 0; // from u
	mesh->inc01[4*0+2*0+1] = 0; // from v
	mesh->inc01[4*0+2*1+0] = 1; // to u
	mesh->inc01[4*0+2*1+1] = 0; // to v
	// v edge
	mesh->inc01[4*1+2*0+0] = 0; // from u
	mesh->inc01[4*1+2*0+1] = 0; // from v
	mesh->inc01[4*1+2*1+0] = 0; // to u
	mesh->inc01[4*1+2*1+1] = 1; // to v
	
	if(fabs(uv) < DBL_EPSILON * uxv){
		mesh->type = 0;
		mesh->n_edges = 2;
		mesh->n_faces = 1;
		
		mesh->inc12[5*0+0] =  1; // u0
		mesh->inc12[5*0+1] = -1; // v0
		mesh->inc12[5*0+2] =  0; // w
		mesh->inc12[5*0+3] = -1; // u1
		mesh->inc12[5*0+4] =  1; // v1
		
		mesh->star0 = uxv;
		mesh->star1[0] = len[1]/len[0];
		mesh->star1[1] = len[0]/len[1];
		mesh->star2[0] = 1./uxv;
	}else{
		mesh->type = 0;
		mesh->n_edges = 3;
		mesh->n_faces = 2;
		
		mesh->star0 = 0.5*uxv;
		mesh->star2[0] = 2./uxv;
		
		// Edge Hodge stars: example for u: dual_ulen/ulen
		//   dual_ulen = circum(0,u,u+v) - circum(0,-v,u)
		//  where circum(0,a,b) = ((|a|^2 b - |b|^2 a) x (axb))/(2 |axb|^2)
		//  This simplifies to
		//   dual_ulen = ||u|^2(v+w)-(w^2+v^2)u| / (2|uxv|)
		//  Generally, this is
		//   dual_*len = (|*|^2(u+v+w) - (v^2+u^2+w^2)*) / (2 |uxv|)
		if(uv < 0){ // wide angle between u and v
			mesh->type = 1;
			Lr[4] = Lr[0] + Lr[2];
			Lr[5] = Lr[1] + Lr[3];
			
			mesh->inc01[4*2+2*0+0] = 0;
			mesh->inc01[4*2+2*0+1] = 0;
			mesh->inc01[4*2+2*1+0] = 1;
			mesh->inc01[4*2+2*1+1] = 1;
		}else{
			mesh->type = 2;
			Lr[4] = Lr[2] - Lr[0];
			Lr[5] = Lr[3] - Lr[1];
			
			mesh->inc01[4*2+2*0+0] = 1;
			mesh->inc01[4*2+2*0+1] = 0;
			mesh->inc01[4*2+2*1+0] = 0;
			mesh->inc01[4*2+2*1+1] = 1;
		}
		len[2] = hypot(Lr[4], Lr[5]);
		const double nlen2[3] = {
			len[0]*len[0] / (2. * uxv),
			len[1]*len[1] / (2. * uxv),
			len[2]*len[2] / (2. * uxv)
		};
		int i;
		for(i = 0; i < 3; ++i){
			const int i1 = (i+1)%3;
			const int i2 = (i+1)%3;
			double q[2] = {
				nlen2[i]*(Lr[2*i1+0]+Lr[2*i2+0]) - (nlen2[i1]+nlen2[i2])*Lr[2*i+0],
				nlen2[i]*(Lr[2*i1+1]+Lr[2*i2+1]) - (nlen2[i1]+nlen2[i2])*Lr[2*i+1]
			};
			mesh->star1[i] = hypot(q[0], q[1]) / len[i];
		}
	}
	
	return 0;
}

// circumcenter of triangle = c = (origin, a, b)
static void circumcenter2(const double a[2], const double b[2], double c[2]){
	const double axb2 = 2.*fabs(a[0]*b[1] - a[1]*b[0]);
	const double a2 = a[0]*a[0] + a[1]*a[1];
	const double b2 = b[0]*b[0] + b[1]*b[1];
	c[1] =-(a2*b[0] - b2*a[0]) / axb2;
	c[0] = (a2*b[1] - b2*a[1]) / axb2;
}

int LubUMesh2_Neighborhood0(const UMesh2 *mesh, double p[12]){
	if(NULL == mesh){ return -1; }
	if(NULL == p){ return -2; }
	
	const double *Lr = mesh->Lr;
	const double *u = &Lr[0];
	const double *v = &Lr[2];
	double q[2];
	switch(mesh->type){
	case 0:
		p[2*0+0] = -0.5*Lr[0]-0.5*Lr[2];
		p[2*0+1] = -0.5*Lr[1]-0.5*Lr[3];
		p[2*1+0] =  0.5*Lr[0]-0.5*Lr[2];
		p[2*1+1] =  0.5*Lr[1]-0.5*Lr[3];
		return 2;
	case 1:
		q[0] = -v[0];
		q[1] = -v[1];
		circumcenter2(q, u, &p[2*0+0]);
		q[0] = u[0]+v[0];
		q[1] = u[1]+v[1];
		circumcenter2(u, q, &p[2*1+0]);
		circumcenter2(q, v, &p[2*2+0]);
		return 3;
	case 2:
		q[0] = u[0]-v[0];
		q[1] = u[1]-v[1];
		circumcenter2(q, u, &p[2*0+0]);
		circumcenter2(u, v, &p[2*1+0]);
		q[0] = -q[0];
		q[1] = -q[1];
		circumcenter2(v, q, &p[2*2+0]);
		return 3;
	default:
		return -1;
	}
}

int LibUMesh2_Neighborhood1(const UMesh2 *mesh, int which, double p[4]){
	if(NULL == mesh){ return -1; }
	if(which < 0 || which >= mesh->n_edges){ return -2; }
	if(NULL == p){ return -3; }
	
	const double *Lr = mesh->Lr;
	const double *u = &Lr[0];
	const double *v = &Lr[2];
	double q[2];
	if(0 == mesh->type){
		if(0 == which){
			p[2*0+0] = 0.5*Lr[0]+0.5*Lr[2];
			p[2*0+1] = 0.5*Lr[1]+0.5*Lr[3];
			p[2*1+0] = 0.5*Lr[0]-0.5*Lr[2];
			p[2*1+1] = 0.5*Lr[1]-0.5*Lr[3];
		}else{
			p[2*0+0] = -0.5*Lr[0]+0.5*Lr[2];
			p[2*0+1] = -0.5*Lr[1]+0.5*Lr[3];
			p[2*1+0] =  0.5*Lr[0]+0.5*Lr[2];
			p[2*1+1] =  0.5*Lr[1]+0.5*Lr[3];
		}
	}else if(1 == mesh->type){
		if(0 == which){
			q[0] = u[0]+v[0];
			q[1] = u[1]+v[1];
			circumcenter2(u, q, &p[2*0+0]);
			q[0] = -v[0];
			q[1] = -v[1];
			circumcenter2(q, u, &p[2*1+0]);
		}else if(1 == which){
			q[0] = u[0]+v[0];
			q[1] = u[1]+v[1];
			circumcenter2(q, v, &p[2*0+0]);
			circumcenter2(u, q, &p[2*1+0]);
		}else{
			q[0] = -u[0];
			q[1] = -u[1];
			circumcenter2(v, q, &p[2*0+0]);
			q[0] = u[0]+v[0];
			q[1] = u[1]+v[1];
			circumcenter2(q, v, &p[2*1+0]);
		}
	}else{
		if(0 == which){
			circumcenter2(u, v, &p[2*0+0]);
			q[0] = u[0]-v[0];
			q[1] = u[1]-v[1];
			circumcenter2(q, u, &p[2*1+0]);
		}else if(1 == which){
			q[0] = v[0]-u[0];
			q[1] = v[1]-u[1];
			circumcenter2(v, q, &p[2*0+0]);
			circumcenter2(u, v, &p[2*1+0]);
		}else{
			circumcenter2(u, v, &p[2*0+0]);
			q[0] = v[0]-u[0];
			q[1] = v[1]-u[1];
			circumcenter2(v, q, &p[2*1+0]);
			p[2] += u[0];
			p[3] += u[1];
		}
	}
	return 0;
}

