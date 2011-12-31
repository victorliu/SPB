#include <libumesh.h>
#include <stdio.h>

int main(int argc, char *argv[]){
	UMesh2 M;
	double p[12];
	double u[2] = {1,0};
	double v[2] = {0.7,1.9};
	
	int i, n;
	LibUMesh2_Create(u, v, &M);
	n = LubUMesh2_Neighborhood0(&M, p);
	for(i = 0; i < n; ++i){
		p[2*(n+i)+0] = -p[2*i+0];
		p[2*(n+i)+1] = -p[2*i+1];
	}
	printf("72 72 scale 4 4 translate 0.01 setlinewidth\n");
	
	printf("%.14g %.14g 0.05 0 360 arc fill\n", 0,0);
	for(i = 0; i < M.n_edges; ++i){
		printf("%.14g %.14g 0.05 0 360 arc fill\n", M.Lr[2*i+0], M.Lr[2*i+1]);
	}
	for(i = 0; i < M.n_edges; ++i){
		printf("%.14g %.14g 0.05 0 360 arc fill\n", -M.Lr[2*i+0], -M.Lr[2*i+1]);
	}
	
	for(i = 0; i < 2*n; ++i){
		printf("%.14g %.14g", p[2*i+0], p[2*i+1]);
		if(0 == i){
			printf(" moveto\n");
		}else{
			printf(" lineto\n");
		}
	}
	printf("closepath stroke\nshowpage\n");
	return 0;
}
