#define SPB_Alloc(TYPE, COUNT, REASON) (TYPE*)SPBmalloc(sizeof(TYPE) * (COUNT), (REASON))
#define SPB_Free(PTR, REASON) SPBfree((PTR), (REASON))
#define STR(x) XSTR(x)
#define XSTR(x) #x

void *SPBmalloc(size_t n, const char *reason);
void SPBfree(void *ptr, const char *reason);
