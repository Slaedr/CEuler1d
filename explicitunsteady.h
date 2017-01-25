/** \file explicitunsteady.h
 * \brief Struct and functions for explicit unsteady solver.
 * \author Aditya Kashi
 * \date January 2017
 */

#ifndef __EXPLICITUNSTEADY_H
#define __EXPLICITUNSTEADY_H

#ifndef __1DEULER_H
#include "1deuler.h"
#endif

#ifndef __RECONSTRUCTION_H
#include "reconstruction.h"
#endif

/// Parameters for explicit RK solver for time-dependent 1D Euler equations
typedef struct
{
	Float ftime;								///< Physical time for which to simulate
	int temporalOrder;							///< desired temporal order of accuracy
	Float** RKCoeffs;							///< Low-storage multi-stage TVD RK coefficients
} Euler1dUnsteadyExplicit;

/// Allocates memory and initializes structure objects with simulation data
void setup_data_unsteady(const size_t num_cells, const int bcleft, const int bcright, const Float bcvalleft[NVARS], const Float bcvalright[NVARS], const Float domain_length, const char *const _flux, 
		const Float _cfl, Float f_time, int temporal_order, const char* RKfile, Grid *const grid, Euler1d *const sim, Euler1dUnsteadyExplicit *const tsim);

void run_unsteady(const Grid *const grid, Euler1d *const sim, Euler1dUnsteadyExplicit *const tsim);

void postprocess_unsteady(const Grid *const grid, const Euler1d *const sim, const char *const outfilename);

void finalize_unsteady(Euler1dUnsteadyExplicit *const tsim);

#endif
