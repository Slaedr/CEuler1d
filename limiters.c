/** \file limiters.c
 * \brief Implementation of slope limiters
 */

#include "limiters.h"

inline Float vanalbada_limiter_function(Float r)
{
	return 2.0*r/(1.0+r*r);
}
