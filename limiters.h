/** \file limiters.h
 * \brief Limiters for TVD/MUSCL schemes
 * \author Aditya Kashi
 * \date Jan 2017
 */

#ifndef __LIMITERS_H

#define __LIMITERS_H

#ifndef __DEFINITIONS_H
#include "definitions.h"
#endif

#pragma acc routine seq
Float vanalbada_limiter_function(Float r);

#endif
