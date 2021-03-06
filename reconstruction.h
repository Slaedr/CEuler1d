/** \file reconstruction.h
 * \brief Schemes for computing the derivatives by reconstructio, and for computing face values from cell-centred values.
 * \author Aditya Kashi
 * \date Jan 2017
 */

#ifndef __RECONSTRUCTION_H
#define __RECONSTRUCTION_H

#ifndef __DEFINITIONS_H
#include "definitions.h"
#endif

#ifndef __LIMITERS_H
#include "limiters.h"
#endif

#ifndef __1DEULER_H
#include "1deuler.h"
#endif

// Computes states at the faces
//void compute_MUSCLReconstruction(const size_t N, Float const *const x, Float const *const *const u, Float * const * const uleft, 
//		Float * const * const uright, const Float k);

/// Piecewise constant "reconstruction" to compute states at faces
void compute_noReconstruction(const Grid* const grid, Euler1d *const sim);

/// Computes states at the faces by MUSCL reconstruction with Van Albada limiter
void compute_MUSCLReconstruction(const Grid *const grid, Euler1d *const sim);

/// Computes face states by MUSCL reconstruction without limiter
void compute_MUSCLReconstruction_unlimited(const Grid *const grid, Euler1d *const sim);

#endif
