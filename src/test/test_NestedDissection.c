#include "NestedDissection.h"
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char* argv[]){
	const int dim = 2;
	const int n[2] = {8, 8};
	const int a[2] = {8, 1};
	const int p[2] = {1, 1};
	int ntot = 1;
	int i[2],d;
	
	for(d = 0; d < dim; ++d){
		ntot *= n[d];
	}
	int *perm = (int*)malloc(sizeof(int) * ntot);
		
	NestedDissectionCartesian1v(dim, n, a, p, perm);
	for(i[0] = 0; i[0] < n[0]; ++i[0]){
		for(i[1] = 0; i[1] < n[1]; ++i[1]){
			int k = 0; for(d = 0; d < dim; ++d){ k += a[d]*i[d]; }
			printf("%d\t%d\t%d\n", i[0], i[1], perm[k]);
		}
	}
	
	for(i[1] = n[1]-1; i[1] >= 0; --i[1]){
		for(i[0] = 0; i[0] < n[0]; ++i[0]){
			int k = 0; for(d = 0; d < dim; ++d){ k += a[d]*i[d]; }
			if(0 != i[0]){ putchar('\t'); }
			printf("%d", perm[k]);
			
			int nk;
			NestedDissectionCartesian1(dim, n, a, p, i, &nk);
			//printf("(%d)", nk);
			if(nk != perm[k]){
				putchar('*');
			}
		} printf("\n");
	}
	
	free(perm);
	return 0;
}
