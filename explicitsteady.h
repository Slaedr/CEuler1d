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

/// Explicit RK solver for steady-state 1D Euler
typedef struct 
{
	Float tol;
	int maxiter;
	Float mws;									///< for computing time steps
} Euler1dSteadyExplicit;

void run_steady(const Grid *const grid, Euler1d *const sim, Euler1dSteadyExplicit *const tsim);

void postprocess_steady(const Grid *const grid, const Euler1d *const sim, const char *const outfilename);

#endif
