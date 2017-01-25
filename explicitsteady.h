/** \file explicitsteady.h
 * \brief Struct and functions for explicit steady-state solver.
 * \author Aditya Kashi
 * \date January 2017
 */

#ifndef __EXPLICITSTEADY_H
#define __EXPLICITSTEADY_H

#ifndef __1DEULER_H
#include "1deuler.h"
#endif

#ifndef __RECONSTRUCTION_H
#include "reconstruction.h"
#endif

/// Explicit RK solver for steady-state 1D Euler
typedef struct 
{
	Float tol;
	int maxiter;
} Euler1dSteadyExplicit;

void setup_data_steady(const size_t num_cells, const int bcleft, const int bcright, const Float bcvalleft[NVARS], const Float bcvalright[NVARS], const Float domain_length, 
		const char *const _flux, const Float _cfl, Float _tol, int max_iter, Grid *const grid, Euler1d *const sim, Euler1dSteadyExplicit *const tsim);

void run_steady(const Grid *const grid, Euler1d *const sim, Euler1dSteadyExplicit *const tsim);

void postprocess_steady(const Grid *const grid, const Euler1d *const sim, const char *const outfilename);

#endif
