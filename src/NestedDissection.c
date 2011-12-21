#include <stdlib.h>
#ifdef DEBUG_NESDIS
# include <stdio.h>
# include <assert.h>
#endif
#define MAX_DIM 3

// returns index i of largest a[i]
static int argmax(int dim, const int a[]){
	int ret = 0;
	int i = 0;
	for(i = 1; i < dim; ++i){
		if(a[i] > a[ret]){ ret = i; }
	}
	return ret;
}

#ifdef DEBUG_NESDIS
static void print_level(int level){
	int i;
	for(i = 0; i < level; ++i){
		putchar(' ');
	}
}
#endif

static void NestedDissectionCartesian1v_r(
	int dim,
	const int n[],
	const int o[],
	const int a[],
	const int p[],
	int *perm,
	int *next,
	int level
){
	int i, m;
	int nsub[MAX_DIM];
	int osub[MAX_DIM];
	int psub[MAX_DIM];
	int maxdim = argmax(dim, n);
	
	// we have to shave off an additional row off the end of
	// the current interval if it is periodic
	int pshave;
	
#ifdef DEBUG_NESDIS
	print_level(level);
	printf("n = {");
	for(i = 0; i < dim; ++i){
		printf(" %d", n[i]);
	} printf(" }\n");
	
	print_level(level);
	printf("o = {");
	for(i = 0; i < dim; ++i){
		printf(" %d", o[i]);
	} printf(" }\n");
#endif
	
	if(1 == n[maxdim]){
		int k = 0;
		for(i = 0; i < dim; ++i){ k += a[i]*o[i]; }
		
#ifdef DEBUG_NESDIS
		print_level(level);
		printf("Assigning {");
		for(i = 0; i < dim; ++i){ printf(" %d",o[i]); }
		printf(" }(k=%d) value %d\n", k, *next);
#endif
		
		perm[k] = (*next);
		(*next)++;
		return;
	}
	for(i = 0; i < dim; ++i){ nsub[i] = n[i]; }
	for(i = 0; i < dim; ++i){ osub[i] = o[i]; }
	for(i = 0; i < dim; ++i){ psub[i] = p[i]; }
	
	pshave = p[maxdim] ? 1 : 0;
	psub[maxdim] = 0;
	
	m = o[maxdim] + n[maxdim]/2;
#ifdef DEBUG_NESDIS
	print_level(level); printf("m = %d\n\n", m);
#endif
	if(m - o[maxdim] > 0){
#ifdef DEBUG_NESDIS
		print_level(level); printf("branch 1:\n", m);
#endif
		nsub[maxdim] = m - o[maxdim];
		osub[maxdim] = o[maxdim];
		NestedDissectionCartesian1v_r(dim, nsub, osub, a, psub, perm, next, level+1);
	}
	if(o[maxdim] + n[maxdim] - (m+1) - pshave > 0){
#ifdef DEBUG_NESDIS
		print_level(level); printf("branch 2:\n", m);
#endif
		nsub[maxdim] = o[maxdim] + n[maxdim] - (m+1) - pshave;
		osub[maxdim] = m+1;
		NestedDissectionCartesian1v_r(dim, nsub, osub, a, psub, perm, next, level+1);
	}
#ifdef DEBUG_NESDIS
	print_level(level); printf("middle:\n", m);
#endif
	nsub[maxdim] = 1;
	osub[maxdim] = m;
	NestedDissectionCartesian1v_r(dim, nsub, osub, a, psub, perm, next, level+1);
	if(pshave){
#ifdef DEBUG_NESDIS
		print_level(level); printf("end:\n", m);
#endif
		nsub[maxdim] = 1;
		osub[maxdim] = o[maxdim] + n[maxdim] - 1;
		NestedDissectionCartesian1v_r(dim, nsub, osub, a, psub, perm, next, level+1);
	}
}

// Assume that the index of element (i,j,k) is
//   idx(i,j,k) = a[0]*i + a[1]*j + a[2]*k
// where 0 <= i <= n[0]-1, and similarly for j and k.
// p[i] should be nonzero if the mesh is periodic in
// the i-th direction.
// On output, perm[q] contains the new index for q=idx(i,j,k)
int NestedDissectionCartesian1v(
	int dim,
	const int n[],
	const int a[],
	const int p[],
	int *perm
){
	int i;
	int is_periodic[MAX_DIM];
	int next = 0;
	int org[MAX_DIM];
	
	if(dim < 1 || dim > MAX_DIM){ return -1; }
	if(NULL == n){ return -2; }
	if(NULL == a){ return -3; }
	if(NULL == p){ return -4; }
	if(NULL == perm){ return -5; }
	
	for(i = 0; i < MAX_DIM; ++i){ is_periodic[i] = p[i]; }
	for(i = 0; i < MAX_DIM; ++i){ org[i] = 0; }
	NestedDissectionCartesian1v_r(dim, n, org, a, is_periodic, perm, &next, 0);
	return 0;
}

static void NestedDissectionCartesian1_r(
	int dim,
	const int n[],
	const int o[],
	const int a[],
	const int p[],
	const int oldind[],
	int *newind,
	int level
){
	int i, m;
	int nsub[MAX_DIM];
	int osub[MAX_DIM];
	int psub[MAX_DIM];
	int maxdim = argmax(dim, n);
	int pshave;
	int prod = 1;
	
	for(i = 0; i < dim; ++i){
		if(i != maxdim){ prod *= n[i]; }
	}
#ifdef DEBUG_NESDIS
	print_level(level);
	printf("n = {");
	for(i = 0; i < dim; ++i){
		printf(" %d", n[i]);
	} printf(" }\n");

	print_level(level);
	printf("o = {");
	for(i = 0; i < dim; ++i){
		printf(" %d", o[i]);
	} printf(" }\n");
#endif
	if(1 == n[maxdim]){
		int k = 0;
		for(i = 0; i < dim; ++i){ k += a[i]*o[i]; }
#ifdef DEBUG_NESDIS
		print_level(level);
		printf("Assigning {");
		for(i = 0; i < dim; ++i){ printf(" %d",o[i]); }
		printf(" }(k=%d) value %d\n", k, *newind);
#endif
		return;
	}
	for(i = 0; i < dim; ++i){ nsub[i] = n[i]; }
	for(i = 0; i < dim; ++i){ osub[i] = o[i]; }
	for(i = 0; i < dim; ++i){ psub[i] = p[i]; }
	
	pshave = p[maxdim] ? 1 : 0;
	psub[maxdim] = 0;
	
	m = o[maxdim] + n[maxdim]/2;
#ifdef DEBUG_NESDIS
	print_level(level); printf("m = %d\n\n", m);
#endif
	if(oldind[maxdim] < m){
#ifdef DEBUG_NESDIS
		assert(m - o[maxdim] > 0);
		print_level(level); printf("branch 1:\n", m);
#endif
		nsub[maxdim] = m - o[maxdim];
		osub[maxdim] = o[maxdim];
		NestedDissectionCartesian1_r(dim, nsub, osub, a, psub, oldind, newind, level+1);
	}else if(oldind[maxdim] > m && oldind[maxdim] < o[maxdim]+n[maxdim]-pshave){
#ifdef DEBUG_NESDIS
		assert(o[maxdim] + n[maxdim] - (m+1) - pshave > 0);
#endif
		*newind += prod*(m - o[maxdim]);
#ifdef DEBUG_NESDIS
		print_level(level); printf("branch 2:\n", m);
#endif
		nsub[maxdim] = o[maxdim] + n[maxdim] - (m+1) - pshave;
		osub[maxdim] = m+1;
		NestedDissectionCartesian1_r(dim, nsub, osub, a, psub, oldind, newind, level+1);
	}else if(oldind[maxdim] == m){
		*newind += prod*(n[maxdim] - 1 - pshave);
#ifdef DEBUG_NESDIS
		print_level(level); printf("middle:\n", m);
#endif
		nsub[maxdim] = 1;
		osub[maxdim] = m;
		NestedDissectionCartesian1_r(dim, nsub, osub, a, psub, oldind, newind, level+1);
	}else{ // if(pshave && oldind[maxdim] == o[maxdim]+n[maxdim]-1){
#ifdef DEBUG_NESDIS
		assert(pshave);
#endif
		*newind += prod*(n[maxdim] - 1);
#ifdef DEBUG_NESDIS
		print_level(level); printf("end:\n", m);
#endif
		nsub[maxdim] = 1;
		osub[maxdim] = o[maxdim] + n[maxdim] - 1;
		NestedDissectionCartesian1_r(dim, nsub, osub, a, psub, oldind, newind, level+1);
	}
}

int NestedDissectionCartesian1(
	int dim,
	const int n[],
	const int a[],
	const int p[],
	const int oldind[],
	int *newind
){
	int i;
	int is_periodic[MAX_DIM];
	int next = 0;
	int org[MAX_DIM];
	
	if(dim < 1 || dim > MAX_DIM){ return -1; }
	if(NULL == n){ return -2; }
	if(NULL == a){ return -3; }
	if(NULL == p){ return -4; }
	if(NULL == oldind){ return -5; }
	if(NULL == newind){ return -6; }
	
	for(i = 0; i < MAX_DIM; ++i){ is_periodic[i] = p[i]; }
	for(i = 0; i < MAX_DIM; ++i){ org[i] = 0; }
	*newind = 0;
	NestedDissectionCartesian1_r(dim, n, org, a, is_periodic, oldind, newind, 0);
	return 0;
}
