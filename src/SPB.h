/* Copyright (C) 2011, Stanford University
 * This file is part of SPB
 * Written by Victor Liu (vkl@stanford.edu)
 *
 * SPB is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * SPB is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/* This header provides SPB's C interface definition. This is version 0
 * of the public interface, and subject to change.
 */

#ifndef _SPB_H_INCLUDED_
#define _SPB_H_INCLUDED_


#if __STDC_VERSION__ >= 199901L /* Use C99 complex type */
#include <complex.h>
typedef _Complex double* SPB_complex_ptr;
#else
typedef double* SPB_complex_ptr;
#endif


typedef struct tag_SPB_LorentzPole{
	double omega_0, Gamma, omega_p;
} SPB_LorentzPole;

typedef struct tag_SPB_ConstitutiveTensor{
	enum{
		SPB_ConstitutiveTensor_SCALAR,
		SPB_ConstitutiveTensor_DIAGONAL,
		SPB_ConstitutiveTensor_TENSOR
	} type;
	SPB_complex_ptr value;
} SPB_ConstitutiveTensor;


typedef struct tag_SPB_BandSolver SPB_BandSolver;

SPB_BandSolver* SPB_BandSolver_New(int dim, int pol, double *Lr);
void            SPB_BandSolver_Destroy(SPB_BandSolver *S);

int SPB_BandSolver_AddMaterial(SPB_BandSolver *S, const char *name, const SPB_ConstitutiveTensor *eps);
int SPB_BandSolver_SetMaterial(SPB_BandSolver *S, const char *name, const SPB_ConstitutiveTensor *eps);
int SPB_BandSolver_Material_AddLorentzPole(SPB_BandSolver *S, const char *name, const SPB_LorentzPole *pole);
int SPB_BandSolver_RemoveMaterial(SPB_BandSolver *S, const char *name);

int SPB_BandSolver_AddRectangle(SPB_BandSolver *S,
	const char *material,
	double center[2],
	double halfwidth[2],
	double angle);

int SPB_BandSolver_SetNumWanted(SPB_BandSolver *S, int n);
int SPB_BandSolver_SetTolerance(SPB_BandSolver *S, double tol);
int SPB_BandSolver_SetResolution(SPB_BandSolver *S, int *res);
int SPB_BandSolver_SetTargetFrequency(SPB_BandSolver *S, double freq);
int SPB_BandSolver_SetTargetFrequencyRange(SPB_BandSolver *S, double freq0, double freq1);

int SPB_BandSolver_SolveK(SPB_BandSolver *S, double *k);
int SPB_BandSolver_GetFrequencies(SPB_BandSolver *S, int *n, SPB_complex_ptr z);
int SPB_BandSolver_GetNumFrequencies(SPB_BandSolver *S);
int SPB_BandSolver_GetBand(SPB_BandSolver *S, int n, SPB_complex_ptr z);

#endif /* _SPB_H_INCLUDED_ */
