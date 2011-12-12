#include "lp_solve.h"
#include "lapack_decl.h"
#include <math.h>
#include <float.h>
#include <string.h>
#include <stdio.h>

// Here we allocate one extra column in A so that the solver does not
// need to reallocate an entire matrix just to perform an augmentation.
int lp_init(lp_t *lp, int n, int m){
	if(NULL == lp){ return -1; }
	if(0 == n){ return -2; }
	if(0 == m){ return -3; }
	if(m >= n){ return -2; }
	
	lp->n = n;
	lp->m = m;
	int alloc_size = (n+1)*m + m + n; // note the n+1, see above
	lp->A = (double*)malloc(alloc_size*sizeof(double));
	if(NULL == lp->A){ return 1; }
	lp->b = lp->A + (n+1)*m; // size m
	lp->c = lp->b + m;       // size n
	
	memset(lp->A, 0, alloc_size*sizeof(double));
	return 0;
}

void lp_destroy(lp_t *lp){
	if(NULL != lp){
		free(lp->A);
	}
}

double lp_eval(const lp_t *lp, double *x){
	return FCALL(ddot)(lp->n, lp->c, 1, x, 1);
}

int lp_in_domain(const lp_t *lp, double *x){
	register int i;
	for(i = 0; i < lp->n; ++i){
		if(x[i] <= 0){ return 0; }
	}
	return 1;
}


////
// Notes:
//   Gradient of objective: c - 1./x
//   Hessian of objective:  1./x.^2   (diagonal)
int lp_centering_newton(const lp_t *lp,
                        const double *x0,
                        const double *barrier_t,
                        double *x_opt, double *nu_opt,
                        int *n_steps)
{
	// Parameters
	const double stopping_newton_decrement = 1e-6;
	const double beta = 0.2; // in (0, 1)
	const double alpha = 0.1; // in (0, 0.5)
	const int max_iterations = 1000; // to prevent out of control loops
	
	if(NULL == lp       ){ return -1; }
	if(NULL == x0       ){ return -2; }
	if(NULL == barrier_t){ return -3; }
	if(NULL == x_opt    ){ return -4; }
	if(NULL == nu_opt   ){ return -5; }
	
	const int m = lp->m, n = lp->n;
	const double kappa = 1.0/(*barrier_t);
	int info, ret = 0;
	register int i, j;
	double f;
	int iter;
	
	double lam2; // lambda^2 (Newton decrement)
	double *work; // the overall work array that contains all of the following
	double *g;    // gradient
	double *Hi;   // Hessian inverse (diagonal)
	double *AH;   // AH = A*H^{-1}
	double *S;    // S = -A*H^{-1}*A' = -A*AH'
	double *dx;   // Newton step
	
	// Set up our temporaries
	work = (double*)malloc((
		 m*(m+n) // for AH,S
		+3*n     // for g,Hi,dx
		) * sizeof(double));
	if(NULL == work){
		ret = 1;
		goto error_work;
	}
	int *iwork; // for pivots
	iwork = (int*)malloc(n*sizeof(int));
	if(NULL == iwork){
		ret = 1;
		goto error_iwork;
	}
	AH = work;    // size m x n
	S = AH + n*m; // size m x m
	g = S + m*m;  // size n
	Hi = g + n;   // size n
	dx = Hi + n;  // size m
	
	FCALL(dcopy)(n, x0, 1, x_opt, 1); // x_opt = x0;
	
	f = lp_eval(lp, x_opt); // function value at previous point
	for(i = 0; i < n; ++i){
		f -= kappa*log(x_opt[i]);
	} // f = c'*x0 - kappa*sum(log(x0))
	for(iter = 0; iter < max_iterations; ++iter){
		//// compute gradient and Hessian
		for(i = 0; i < n; ++i){
			g[i] = lp->c[i] - kappa/x_opt[i];
			Hi[i] = (x_opt[i]*x_opt[i])/kappa; // diagonal of Hessian inverse
		}
		
		// Compute Newton step by block elimination of KKT system:
		// [ H  A' ] [ dx ] == [ -g ]   (g is gradient, H is Hessian)
		// [ A  0  ] [ nu ]    [  0 ]
		// H*dx + A'*nu == -g  -->  dx == -H\(g + A'*nu)
		// A*dx == 0
		// Together giving A*H^{-1}*g == -A'*H^{-1}*A*nu == S*nu
		// Therefore nu can be found first from the above symmetric system
		// Then dx can be back solved.
		for(j = 0; j < n; ++j){
			FCALL(dcopy)(m, lp->A+m*j, 1, AH+m*j, 1);
			FCALL(dscal)(m, Hi[j], AH+m*j, 1);
		} // AH = A*H^{-1}
		// FCALL(dgemm)('N','T', m, m, n, -1.0, lp->A, m,
		//              AH, m, 0.0, S, m); // S = -A*AH'
		FCALL(dsyr2k)('L', 'N', m, n, -0.5, lp->A, m,
		              AH, m, 0.0, S, m); // S = -A*AH'
		FCALL(dgemv)('N', m, n, 1.0, AH, m, g, 1,
		             0.0, nu_opt, 1); // AHg = AH*g
		// Solves S*nu = AHg
		FCALL(dsysv)('L', m, 1, S, m, iwork, nu_opt, m,
		             AH /*use as scratch*/, n,
					 &info);
		if(0 != info){ ret = 3; goto error_solve; }
		// Solves H*v = -A'*nu - g
		FCALL(dgemv)('T', m, n, -1.0, lp->A, m, nu_opt, 1,
		             0.0, dx, 1); // dx = -A'*nu
		for(i = 0; i < n; ++i){
			dx[i] = Hi[i]*(dx[i] - g[i]);
		}

		lam2 = -FCALL(ddot)(n, dx, 1, g, 1); // lam2 = dx'*g
		
		// At this point, g is preserved, Hi no longer needed;
		double *xp = Hi; // xp = x + t*dx (reusing storage)
		
		//// stopping criterion
		if(lam2 <= 2*stopping_newton_decrement){ break; }
		
		//// Line search
		double t = 1.0;
		int found = 0;
		while(0 == found){
			if(0 == t){ ret = 3; break; }
			FCALL(dcopy)(n, x_opt, 1, xp, 1);
			FCALL(daxpy)(n, t, dx, 1, xp, 1); // xp = x + t*dx;
			if(0 == lp_in_domain(lp, xp)){ t *= beta; continue; }
			double fp = lp_eval(lp, xp);
			for(i = 0; i < n; ++i){
				fp -= kappa*log(xp[i]);
			} // fp = c'*xp - kappa*sum(log(xp))
			if(fp >= f - alpha*t*lam2){ t *= beta; continue; }
			else{
				found = 1;			
				//// Update
				FCALL(dcopy)(n, xp, 1, x_opt, 1); // x = xp;
				f = fp; // update the objective value
				break;
			}
		}
	}
	if(NULL != n_steps){ *n_steps = iter; }

	if(iter == max_iterations){
		ret = 2;
	}
//	printf("x = "); print_vector(n, x_opt); printf("\n");
//	printf("nu = "); print_vector(m, nu_opt); printf("\n");
//	printf("f = %f\n", f);

error_solve:
	free(iwork);
error_iwork:
	free(work);
error_work:
	return ret;
}


int lp_solve_with_feasible_starting_point(const lp_t *lp,
                                          const double *x0,
                                          double *x_opt,
                                          double *lambda_opt, double *nu_opt,
                                          int *n_steps)
{
	// Parameters
	const double mu = 16;
	const double tolerance = 1e-10;
	const double initial_objective_weight = 1.0;
	const int max_iterations = 1000; // to prevent out of control loops
	
	if(NULL == lp){ return -1; }
	if(NULL == x0){ return -2; }
	if(NULL == x_opt){ return -3; }
	if(NULL != n_steps){ *n_steps = 0; }
	
	int info, ret = 0;
	const int m = lp->m, n = lp->n;
	int i;
	double t = initial_objective_weight;
	
	double *work;
	double *xp, *nu = NULL;
	int lwork = n;
	if(NULL == nu_opt){ lwork += m; } // If nu_opt was NULL, use allocated nu
	else{ nu = nu_opt; }              // otherwise use the provided buffer
	work = (double*)malloc(lwork * sizeof(double));
	if(NULL == work){ return 1; }
	xp = work; // size n
	if(NULL == nu_opt){ nu = xp + n; }
	
	FCALL(dcopy)(n, x0, 1, xp, 1);
	
	int iter;
	for(iter = 0; iter < max_iterations; ++iter){
		//// Centering step
		int ns = 0;
		info = lp_centering_newton(lp, xp, &t, x_opt, nu, &ns);
		if(0 != info){ ret = info; break; }
		if(NULL != n_steps){ *n_steps += ns; }
		
		//// Update
		FCALL(dcopy)(n, x_opt, 1, xp, 1);
		
		//// Stopping criterion
		if((double)n/t < tolerance){ break; }
		t *= mu;
	}
	if(0 == ret && iter == max_iterations){
		ret = 2;
	}
	
	// Compute the dual variables if requested
	t = 1.0/t;
	if(NULL != lambda_opt){
		for(i = 0; i < n; ++i){
			lambda_opt[i] = -t/x_opt[i];
		}
	}
	if(NULL != nu_opt){
		for(i = 0; i < m; ++i){
			nu_opt[i] = t*nu[i];
		}
	}
	
	free(work);
	return ret;
}

int lp_solve(const lp_t *lp,
             double *x,
             double *lambda, double *nu,
             enum lp_solution_status* status,
             int *n_steps)
{
	if(NULL == lp){ return -1; }
	if(NULL == x){ return -2; }
	if(NULL == status){ return -5; }
	if(NULL != n_steps){ *n_steps = 0; }
	
	int info, ns, ret = 0;
	const int m = lp->m, n = lp->n;
	const int mn = m*n;
	register int i, j;
	
	double *work;
	double *b_new, *c_new, *Acopy, *zt, *zt_opt;
	int lwork = m + (n+1) + mn + 2*(n+1);
	work = (double*)malloc(lwork * sizeof(double));
	if(NULL == work){ return 1; }
	b_new = work;          // size m
	c_new = b_new + m;     // size n+1
	Acopy = c_new + (n+1); // size mn
	zt = Acopy + mn;       // size n+1
	zt_opt = zt + (n+1);   // size n+1
	
	// Phase I:
	//   Solve:
	//     minimize   t
	//     subject to A*x == b
	//                x >= (1-t)*ones(n,1), t >= 0
	//   If t < 1, x is strictly feasible
	//   Change variables:
	//     z = x + (t-1)*ones
	//     A*x == b  becomes  A*z == b+(t-1)*A*ones
	//   Modified problem:
	//     minimize   t
	//     subject to [A, -A*ones] * [z; t] == b - A*ones(n,1)
	//                [z; t] >= 0
	//   Any x0 for which A*x0 == b is feasible.
	//   Then if there is an x_i < 0, choose t0 = 2-min_i(x_i)
	//   Otherwise choose t0 = 1
	// We will rely on lp_init to allocate an augmented A matrix
	// so that we don't have to do it here.
	lp_t lp1;
	lp1.n = n+1;
	lp1.m = m;
	lp1.A = lp->A;
	lp1.b = b_new;
	lp1.c = c_new;
	// Set the new b vector
	FCALL(dcopy)(n, lp->b, 1, lp1.b, 1);
	for(i = 0; i <= n; ++i){ // set all of c to be 1
		lp1.c[i] = 1.0;
	}
	const char fcN = 'N';
	FCALL(dgemv)(fcN, m, n, -1.0, lp->A, m, lp1.c, 1,
	             1.0, lp1.b, 1); // b <- -A*ones + b
	memset(lp1.c, 0, n*sizeof(double)); // set the first n elements to zero
	for(i = 0; i < m; ++i){ // set the last column of the augmented A
		lp1.A[i+mn] = 0;
		for(j = 0; j < n; ++j){
			lp1.A[i+mn] -= lp1.A[i+m*j];
		}
	}
	// Find an x0 (place it in x for now)
	double *dgels_work;
	int dgels_lwork = -1;
	FCALL(dgels)(fcN, m, n, 1, Acopy, m, zt, n,
	             Acopy, dgels_lwork, &info); // workspace query
	dgels_lwork = (int)Acopy[0];
	dgels_work = (double*)malloc(dgels_lwork*sizeof(double));
	if(NULL == dgels_work){ ret = 1; goto error_dgels_work; }

	FCALL(dcopy)(mn, lp->A, 1, Acopy, 1);
	FCALL(dcopy)(m, lp->b, 1, zt, 1);
	// Solves A*xt == b, least norm solution
	FCALL(dgels)('N', m, n, 1, Acopy, m, zt, n,
	             dgels_work, dgels_lwork,
				 &info);
	if(0 != info){ ret = 3; goto error_rank; }
	
	zt[n] = 1.0;
	for(i = 0; i < n; ++i){ // determine t
		if(zt[i] < 0){
			double temp = 2.0-zt[i];
			if(temp > zt[n]){ zt[n] = temp; }
		}
	}
	for(i = 0; i < n; ++i){ // transform x into z
		zt[i] += (zt[n]-1.0);
	}
	
	// Solve the feasibility test LP
	ns = 0;
	*status = LP_INFEASIBLE;
	info = lp_solve_with_feasible_starting_point(&lp1, zt, zt_opt,
	                                             NULL, NULL, &ns);
	if(0 == info){
		if(NULL != n_steps){ *n_steps += ns; }
		if(zt_opt[n] < 1.0){ // strictly feasible
			// transform optimal z into x, x is the strictly feasible x0
			for(i = 0; i < n; ++i){
				zt_opt[i] -= (zt_opt[n]-1.0);
			}
			// Solve the actual problem
			ns = 0;
			info = lp_solve_with_feasible_starting_point(lp, zt_opt, x,
			                                             lambda, nu, &ns);
			if(0 != info){
				ret = info;
			}else{
				if(NULL != n_steps){ *n_steps += ns; }
				*status = LP_SOLVED;
			}
		}
	}else{
		ret = info;
	}
	
error_rank:
	free(dgels_work);
error_dgels_work:
	free(work);
	return ret;
}

int lp_check_KKT(const lp_t *lp, double *x, double *nu){
	const int m = lp->m, n = lp->n;
	double nrm, *temp;
	register int i;
	const double tol = 8*(double)n*DBL_EPSILON;
	
	temp = (double*)malloc(n*sizeof(double));
	if(NULL == temp){ return 1; }
	
	// The KKT condition is:
	//   Ax = b
	//   A'nu + c - 1./x = 0
	// First check the condition A*x == b:
	FCALL(dcopy)(m, lp->b, 1, temp, 1); // temp = b
	// temp = b - A*x
	FCALL(dgemv)('N', m, n, 1.0, lp->A, m, x, 1, -1.0, temp, 1);
	// temp should be small
	nrm = FCALL(dnrm2)(m, temp, 1);
	if(nrm > tol){
		return -1;
	}
	
	// Now check the second condition
	// The Hessian is 1./(x.*x), so H*x = 1./x
	for(i = 0; i < n; ++i){
		temp[i] = lp->c[i] - 1.0/x[i];
	}
	// temp = A'*nu + c - H*x
	FCALL(dgemv)('T', m, n, 1.0, lp->A, m, nu, 1, 1.0, temp, 1);
	// temp should be small
	nrm = FCALL(dnrm2)(n, temp, 1);
	if(nrm > tol){
		return -1;
	}
	
	free(temp);
	return 0;
}
