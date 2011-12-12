#ifndef _LP_SOLVER_H_ // written by Victor Liu
#define _LP_SOLVER_H_

#include <stdlib.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

// Represents a standard form linear program
//     minimize   (c' * x)
//     subject to A*x == b, x >= 0 (elementwise)
// x is the variable, of dimension n. A has m rows, with m < n.
// It is assumed that A is full rank, and that the sublevel sets
// {x | A*x == b, x >= 0, c'*x <= gamma} are all bounded.

typedef struct linear_program_{
	int n, m;
	double *A, *b, *c;
	// Note that A is assumed to be in column-major (Fortran) order.
	// i.e. A_{i,j} is located at A[m*j+i]
	// Also, A is allocated with an extra column to fascilitate
	// implementation of lp_solve()
} lp_t;

// Returns negative number for errors in argument (-value is which argument)
// Returns positive number for allocation error
// Returns 0 on success
int  lp_init(lp_t *lp, int n, int m);
void lp_destroy(lp_t *lp);

double lp_eval(const lp_t *lp, double *x); // evaluates the objective function
// Returns nonzero if in domain of inequality constraints
int lp_in_domain(const lp_t *lp, double *x);

enum lp_solution_status{
	LP_SOLVED,
	LP_INFEASIBLE, // Actually, just not strictly feasible
	LP_UNBOUNDED   // Currently never used
};

// Purpose:
//   Solves the LP
//     minimize   c'*x
//     subject to A*x == b, x >= 0 (element-wise)
//   with variable x of length n, and m equality constraints with m < n
// Inputs:
//   lp     : Linear program with full rank A (nonsingular KKT matrix)
// Outputs:
//   x      : Primal optimal point (assumes preallocated)
//   lambda : Optional dual variable to inequality constraints
//   nu     : Optional dual variable to equality constraints
//   status : The solution status; only valid if this function returned 0
//            If the problem is not strictly feasible, LP_INFEASIBLE is
//            returned. Does not detect unbounded problems, so LP_UNBOUNDED
//            is never returned.
// Return value:
//   0      : Success
//   < 0    : Error in a parameter. The negative gives which one.
//   1      : Memory allocation error
//   2      : Internal maximum iteration limit reached, convergence failed.
//   3      : The A matrix is not full rank.
int lp_solve(const lp_t *lp,
             double *x,
             double *lambda, double *nu,
             enum lp_solution_status* status, int *n_steps);

///////////////////////////////////////////////////////////////////////////////
/////////////////////////////// Private methods ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// Purpose:
//   Uses Newton's method to solve the centering problem
//     minimize   c'*x - sum_{i=0}^n log x_i
//     subject to A*x == b
//   with variable x of length n, given a strictly feasible starting point x0.
// Inputs:
//   lp        : Linear program with full rank A (nonsingular KKT matrix)
//   x0        : Initial strictly feasible starting point
//   barrier_t : The (inverse) weight of the barrier:
//                 f0(x) = barrier_t * c'*x - sum_{i=0}^n log x_i
//               is the function that is minimized.
// Outputs:
//   x_opt     : Primal optimal point (assumes preallocated)
//   nu_opt    : A dual optimal point (assumes preallocated)
//   n_steps   : The number of Newton steps if this parameter was not NULL.
// Return value:
//   0         : Success
//   < 0       : Error in a parameter. The negative gives which one.
//   1         : Memory allocation error
//   2         : Internal maximum iteration limit reached, convergence failed.
//   3         : An error occured in symmetric positive definite matrix solve.
int lp_centering_newton(const lp_t *lp,
                        const double *x0,
                        const double *barrier_t,
                        double *x_opt, double *nu_opt,
                        int *n_steps);

// Purpose:
//   Checks the solution x and the dual variable nu against the KKT condition.
// Return value:
//   0       : Valid solution
//  -1       : Invalid solution
//   1       : Memory allocation error
int lp_check_KKT(const lp_t *lp, double *x, double *nu);

// Purpose:
//   Solves
//     minimize   c'*x
//     subject to A*x == b, x >= 0 (element-wise)
//   with variable x of length n, given a strictly feasible starting point x0.
// Inputs:
//   lp         : Linear program with full rank A (nonsingular KKT matrix)
//   x0         : Initial strictly feasible starting point
// Outputs:
//   x_opt      : Primal optimal point (assumes preallocated)
//   lambda_opt : Optional dual variable to inequality constraints
//   nu_opt     : Optional dual variable to equality constraints
// Return value:
//   0   : Success
//   < 0 : Error in a parameter. The negative gives which one.
//   1   : Memory allocation error
//   2   : Internal maximum iteration limit reached, convergence failed.
int lp_solve_with_feasible_starting_point(const lp_t *lp,
                                          const double *x0,
                                          double *x_opt,
                                          double *lambda_opt, double *nu_opt,
                                          int *n_steps);

#ifdef __cplusplus
}
#endif

#endif // _LP_SOLVER_H_
