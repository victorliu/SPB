// Nested dissection ordering for a radius 1 stencil
// (5 point stencil in 2D, 7 point stencil in 3D) on
// a regular cartesian grid.

// Assume that the index of element (i,j,k) is
//   idx(i,j,k) = a[0]*i + a[1]*j + a[2]*k
// where 0 <= i <= n[0]-1, and similarly for j and k.
// p[i] should be nonzero if the mesh is periodic in
// the i-th direction.
// On output, perm[q] contains the new index for q=idx(i,j,k)
// Returns 0 on success, -i if the i-th argument is invalid.
int NestedDissectionCartesian1v(
	int dim,
	const int n[],
	const int a[],
	const int p[],
	int *perm
);

// Obtain the permuted index newind for just a single tuple
// oldind=(i,j,k). Performs essentially a binary search, so
// this is quite efficient if storing the entire permutation
// array is not practical.
// Returns 0 on success, -i if the i-th argument is invalid.
int NestedDissectionCartesian1(
	int dim,
	const int n[],
	const int a[],
	const int p[],
	const int oldind[],
	int *newind
);
