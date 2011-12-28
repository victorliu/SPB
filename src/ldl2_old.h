#define CPRINT(A) do{ \
		printf("%.14g+i%.14g", (A)[0], (A)[1]); \
	}while(0)
#define C2PRINT(A) do{ \
		printf("{"); CPRINT(&((A)[2*0])); printf(" "); CPRINT(&((A)[2*2])); printf("}\n"); \
		printf("{"); CPRINT(&((A)[2*1])); printf(" "); CPRINT(&((A)[2*3])); printf("}\n"); \
	}while(0)
#define DABS(A) (((A)<0.)?(-(A)):(A))
#define CMUL(A,B,RES) do{ \
		(RES)[0] = (A)[0]*(B)[0]-(A)[1]*(B)[1]; \
		(RES)[1] = (A)[0]*(B)[1]+(A)[1]*(B)[0]; \
	}while(0)
#define CADD(A,B) do{ \
		(A)[0] += (B)[0]; \
		(A)[1] += (B)[1]; \
	}while(0)
#define CSUB(A,B) do{ \
		(A)[0] -= (B)[0]; \
		(A)[1] -= (B)[1]; \
	}while(0)
#define CADDPROD(C,A,B) do{ \
		(C)[0] += ((A)[0]*(B)[0]-(A)[1]*(B)[1]); \
		(C)[1] += ((A)[0]*(B)[1]+(A)[1]*(B)[0]); \
	}while(0)
#define CSUBPROD(C,A,B) do{ \
		(C)[0] -= ((A)[0]*(B)[0]-(A)[1]*(B)[1]); \
		(C)[1] -= ((A)[0]*(B)[1]+(A)[1]*(B)[0]); \
	}while(0)
#define CSUBCPROD(C,A,B) do{ \
		(C)[0] -= ((A)[0]*(B)[0]+(A)[1]*(B)[1]); \
		(C)[1] -= ((A)[0]*(B)[1]-(A)[1]*(B)[0]); \
	}while(0)
#define CSET(A,B) do{ \
		(A)[0] = (B)[0]; \
		(A)[1] = (B)[1]; \
	}while(0)
#define CNEG(A) do{ \
		(A)[0] = -(A)[0]; \
		(A)[1] = -(A)[1]; \
	}while(0)
#define CSWP(A,B) do{ \
		double t; \
		t = (A)[0]; \
		(A)[0] = (B)[0]; \
		(B)[0] = t; \
		t = (A)[1]; \
		(A)[1] = (B)[1]; \
		(B)[1] = t; \
	}while(0)
/*
1/(R+iI) = R-iI/(RR+II)
1/(R+I/R) - iI/R/(R+I/R)
R/I/(I+R/I) - i/(I+R/I)
*/
#define CINV(A) do{ \
		double r,d; \
		if(DABS(A[0]) > DABS(A[1])){ \
			r = (A)[1]/(A)[0]; \
			d = r + (A)[0]; \
			(A)[0] = 1./d; \
			(A)[1] = -r/d; \
		}else{ \
			r = (A)[0]/(A)[1]; \
			d = r + (A)[1]; \
			(A)[0] = r/d; \
			(A)[1] = -1./d; \
		} \
	}while(0)
#define C2ZERO(A) do{ \
		(A)[0] = 0.; (A)[1] = 0.; \
		(A)[2] = 0.; (A)[3] = 0.; \
		(A)[4] = 0.; (A)[5] = 0.; \
		(A)[6] = 0.; (A)[7] = 0.; \
	}while(0)
#define C2ADD(A,B) do{ \
		(A)[0] += (B)[0]; (A)[1] += (B)[1]; \
		(A)[2] += (B)[2]; (A)[3] += (B)[3]; \
		(A)[4] += (B)[4]; (A)[5] += (B)[5]; \
		(A)[6] += (B)[6]; (A)[7] += (B)[7]; \
	}while(0)
#define C2SET(A,B) do{ \
		(A)[0] = (B)[0]; (A)[1] = (B)[1]; \
		(A)[2] = (B)[2]; (A)[3] = (B)[3]; \
		(A)[4] = (B)[4]; (A)[5] = (B)[5]; \
		(A)[6] = (B)[6]; (A)[7] = (B)[7]; \
	}while(0)
#define C2SETH(A,B) do{ \
		(A)[0] = (B)[0]; (A)[1] = -(B)[1]; \
		(A)[2] = (B)[4]; (A)[3] = -(B)[5]; \
		(A)[4] = (B)[2]; (A)[5] = -(B)[3]; \
		(A)[6] = (B)[6]; (A)[7] = -(B)[7]; \
	}while(0)
#define C2INV(A) do{ \
		double idet[2], t[2]; \
		CMUL(&((A)[2*0]), &((A)[2*3]), idet); \
		CMUL(&((A)[2*1]), &((A)[2*2]), t); \
		CSUB(idet, t); CINV(idet); \
		CSWP(&((A)[2*0]), &((A)[2*3])); \
		CNEG(&((A)[2*1])); \
		CNEG(&((A)[2*2])); \
		CMUL(idet, &((A)[2*0]), t); CSET(&((A)[2*0]), t); \
		CMUL(idet, &((A)[2*1]), t); CSET(&((A)[2*1]), t); \
		CMUL(idet, &((A)[2*2]), t); CSET(&((A)[2*2]), t); \
		CMUL(idet, &((A)[2*3]), t); CSET(&((A)[2*3]), t); \
	}while(0)
// A -= B*C
#define C2SUBPROD(A,B,C) do{ \
		CSUBPROD(&((A)[2*0]),&((B)[2*0]),&((C)[2*0])); \
		CSUBPROD(&((A)[2*0]),&((B)[2*2]),&((C)[2*1])); \
		CSUBPROD(&((A)[2*1]),&((B)[2*1]),&((C)[2*0])); \
		CSUBPROD(&((A)[2*1]),&((B)[2*3]),&((C)[2*1])); \
		CSUBPROD(&((A)[2*2]),&((B)[2*0]),&((C)[2*2])); \
		CSUBPROD(&((A)[2*2]),&((B)[2*2]),&((C)[2*3])); \
		CSUBPROD(&((A)[2*3]),&((B)[2*1]),&((C)[2*2])); \
		CSUBPROD(&((A)[2*3]),&((B)[2*3]),&((C)[2*3])); \
	}while(0)
#define C2MUL(A,B,C) do{ \
		CMUL    (&((A)[2*0]), &((B)[2*0]), &((C)[2*0])); \
		CADDPROD(&((C)[2*0]), &((A)[2*2]), &((B)[2*1])); \
		CMUL    (&((A)[2*1]), &((B)[2*0]), &((C)[2*1])); \
		CADDPROD(&((C)[2*1]), &((A)[2*3]), &((B)[2*1])); \
		CMUL    (&((A)[2*0]), &((B)[2*2]), &((C)[2*2])); \
		CADDPROD(&((C)[2*2]), &((A)[2*2]), &((B)[2*3])); \
		CMUL    (&((A)[2*1]), &((B)[2*2]), &((C)[2*3])); \
		CADDPROD(&((C)[2*3]), &((A)[2*3]), &((B)[2*3])); \
	}while(0)

/* ========================================================================== */
/* === ldl.h:  include file for the LDL package ============================= */
/* ========================================================================== */

/* LDL Copyright (c) Timothy A Davis,
 * University of Florida.  All Rights Reserved.  See README for the License.
 */


/* ========================================================================== */
/* === int version ========================================================== */
/* ========================================================================== */

void LDL_symbolic(
	int n,		/* A and L are n-by-n, where n >= 0 */
	
	int max_col_nnz,
	void (*Acol)(int col, int *nnz, int *row_ind, void *data),

	int Lp[],	/* output of size n+1, not defined on input */
	int Parent[],	/* output of size n, not defined on input */
	// on output, Lnz contains number of non-zeros in each column
	int Lnz[],	/* output of size n, not defined on input */
	
	int Flag[],	/* workspace of size n, not defn. on input or output */
	
	int rowind[], // workspace of size max_col_nnz
	void *data
);

int LDL_numeric(	/* returns n if successful, k if D (k,k) is zero */
	int n,		/* A and L are n-by-n, where n >= 0 */
	
	int max_col_nnz,
	void (*Acol)(int col, int *nnz, int *row_ind, double *row_val, void *data),

	const int Lp[],	/* input of size n+1, not modified */
	int Parent[],	/* input of size n, not modified */
	
	int Lnz[],	/* output of size n, not defn. on input */
	int Li[],	/* output of size lnz=Lp[n], not defined on input */
	double Lx[],	/* output of size lnz=Lp[n], not defined on input */
	double D[],	/* output of size n, not defined on input */
	
	double Y[],	/* workspace of size 8*n, not defn. on input or output */
	int Pattern[],/* workspace of size n, not defn. on input or output */
	int Flag[],	/* workspace of size n, not defn. on input or output */
	int rowind[], // workspace size max_col_nnz
	double rowval[], // workspace size 8*max_col_nnz
	
	void *data
);

void LDL_lsolve(int n, double X[], int Lp[], int Li[], double Lx[]);
void LDL_dsolve(int n, double X[], double D[]);
void LDL_ltsolve(int n, double X[], int Lp[], int Li[], double Lx[]);

