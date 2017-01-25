/** \file reconstruction.h
 * \brief Schemes for computing the derivatives by reconstructio, and for computing face values from cell-centred values.
 * \author Aditya Kashi
 * \date Jan 2017
 */

#ifndef __RECONSTRUCTION_H

#ifndef __DEFINITIONS_H
#include "definitions.h"
#endif

#ifndef __LIMITERS_H
#include "limiters.h"
#endif

/// (Kernal function) Computes states at the faces
void compute_MUSCLReconstruction(const size_t N, Float const *const x, Float const *const *const u, Float * const * const uleft, 
		Float * const * const uright, const Float k);

#endif
