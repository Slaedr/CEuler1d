/** \file inviscidflux.h
 * \brief Approximate Riemann solvers for Euler equations
 * \author Aditya Kashi
 * \date Jan 2017
 */

#ifndef __INVISCIDFLUX_H
#define __INVISCIDFLUX_H

#ifndef __DEFINITIONS_H
#include "definitions.h"
#endif

#pragma acc routine seq
void compute_llfflux_prim(Float const *const uleft, Float const *const uright, Float *const flux, const Float g);

#pragma acc routine seq
void compute_vanleerflux_prim(Float const *const uleft, Float const *const uright, Float *const flux, const Float g);

#endif

