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

/// Interface for computing face values from cell-centred values
class FaceReconstruction
{
protected:
	const int N;										///< Number of cells
	Float const *const x;						///< Positions of cell centres
	Float const *const *const u;			///< Cell-centred variables
	Float const *const *const dudx;		///< Cell-centred slopes
	Float * const * const uleft;			///< Left value at each face
	Float * const * const uright;			///< Right value at each face
public:
	FaceReconstruction(const int _N, Float const *const x, Float const *const *const _u, Float const *const *const _dudx, Float * const * const uleft,
			Float * const * const uright);
	virtual ~FaceReconstruction();
	virtual void compute_face_values() = 0;
};

class MUSCLReconstruction : public FaceReconstruction
{
	const Float k;										///< Controls order of reconstruction; people generally use 1/3
	std::string limiter;								///< String describing the limiter to use
	const Limiter* lim;									///< Slope limiter to use
public:
	MUSCLReconstruction(const int _N, Float const *const x, Float const *const *const _u, Float const *const *const _dudx, Float * const * const uleft,
			Float * const * const uright, std::string _limiter, const Float _k);
	~MUSCLReconstruction();
	void compute_face_values();
};
	
void compute_MUSCLReconstruction(const int N, Float const *const x, Float const *const *const u, Float const *const *const _dudx, Float * const * const uleft, 
		Float * const * const uright, const Float k);

#endif
