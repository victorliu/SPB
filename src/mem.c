#include "stdio.h"
#include "stdlib.h"
void *SPBmalloc(size_t n, const char *reason){
	void *ptr = malloc(n);
//	fprintf(stderr, "Allocating %u bytes at %p: %s\n", n, ptr, reason); fflush(stderr);
	return ptr;
}
void SPBfree(void *ptr, const char *reason){
//	fprintf(stderr, "Freeing %p: %s\n", ptr, reason); fflush(stderr);
	free(ptr);
}
