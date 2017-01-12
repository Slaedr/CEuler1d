/** \file 1deuler.h
 * \brief Code to control the solution process for time-accurate 1D Euler equations.
 * \author Aditya Kashi
 * \date Jan 2017
 */

#ifndef __1DEULER_H
#define __1DEULER_H

#ifndef __INVISCIDFLUX_H
#include "inviscidflux.h"
#endif

#ifndef __RECONSTRUCTION_H
#include "reconstruction.h"
#endif

/// Base class for Euler solution processes
class Euler1d
{
protected:
	int N;							///< Number of real cells in the grid
	int ncell;						///< total number of cells including ghost cells
	int nface;						///< total number of interfaces to compute fluxes across
	Float* x;						///< Cell centers
	Float* dx;						///< (1D) Size of each cell
	Float* vol;					///< (3D) Volume of each cell
	Float* nodes;					///< Mesh nodes
	Float domlen;					///< Physical length of the domain
	Float* A;						///< Cross-sectional areas at cell centers
	Float* Af;						///< Cross-sectional areas at interfaces
	Float** u;						///< Conserved variables - u[i][0] is density of cell i and so on
	Float** prim;					///< Primitive variables - u[i][0] is density of cell i and so on
	Float** uleft;					///< Left state of each face
	Float** uright;				///< Right state at each face
	Float** prleft;				///< Left state of each face in terms of primitive variables
	Float** prright;				///< Right state at each face in terms of primitive variables
	Float** dudx;					///< Slope of variables in each cell
	Float** res;					///< residual
	int bcL;									///< left BC type
	int bcR;									///< right BC type
	Float bcvalL[NVARS];						///< left boundary value
	Float bcvalR[NVARS];						///< right boundary value
	InviscidFlux* flux;							///< Inviscid flux computation context
	SlopeReconstruction* cslope;				///< Slope reconstruction context
	FaceReconstruction* rec;					///< Context responsible for computation of face values of flow variables from their cell-centred values
	Float cfl;									///< CFL number

public:
	Euler1d(int num_cells, Float length, int leftBCflag, int rightBCflag, std::vector<Float> leftBVs, std::vector<Float> rightBVs, Float cfl, 
			std::string inviscid_flux, std::string slope_scheme, std::string face_extrap_scheme, std::string limiter);

	virtual ~Euler1d();

	/// Generates a grid depending on 
	/** \param type If type == 0, a uniform grid is generated and used. If type == 1, then grid points need to passed to the function in
	 * \param pointlist an array of positions of mesh points.
	 *
	 * Note that the domain is assumed to start at x=0.
	 */
	void generate_mesh(int type, const std::vector<Float>& pointlist);

	/// Set cross-sectional areas
	void set_area(int type, std::vector<Float>& cellCenteredAreas);

	void compute_slopes();

	void compute_face_values();

	void compute_inviscid_fluxes(Float** prleft, Float** prright, Float** res, Float* Af);
	
	void compute_inviscid_fluxes_cellwise(Float** prleft, Float** prright, Float** res, Float* Af);

	void compute_source_term(Float** u, Float** res, Float* Af);

	/// Find new ghost cell values
	void apply_boundary_conditions();
	
	/// Find new values of left boundary face external state
	/** Note that interior states at boundary faces should already be computed.
	 */
	void apply_boundary_conditions_at_left_boundary(std::vector<Float>& ul, const std::vector<Float>& ur);
	
	/// Find new values of right boundary face external state
	/** Note that interior states at boundary faces should already be computed.
	 */
	void apply_boundary_conditions_at_right_boundary(const std::vector<Float>& ul, std::vector<Float>& ur);
};

/// Explicit RK solver for time-dependent 1D Euler equations
class Euler1dExplicit : public Euler1d
{
	Float ftime;								///< Physical time for which to simulate
	Float mws;									///< for computing time steps
	Float a;									///< temp variable
	int temporalOrder;							///< desired temporal order of accuracy
	Float** RKCoeffs;							///< Low-storage multi-stage TVD RK coefficients

public:
	Euler1dExplicit(int num_cells, Float length, int leftBCflag, int rightBCflag, std::vector<Float> leftBVs, std::vector<Float> rightBVs, Float cfl, std::string inviscidFlux,
			std::string slope_scheme, std::string face_rec_scheme, std::string limiter, Float fTime, int temporal_order, std::string RKfile);

	~Euler1dExplicit();

	void run();
	
	void postprocess(std::string outfilename);
};

/// Explicit RK solver for steady-state 1D Euler
class Euler1dSteadyExplicit : public Euler1d
{
	Float tol;
	int maxiter;
	Float mws;									///< for computing time steps

public:
	Euler1dSteadyExplicit(int num_cells, Float length, int leftBCflag, int rightBCflag, std::vector<Float> leftBVs, std::vector<Float> rightBVs, Float cfl, std::string inviscidFlux,
			std::string slope_scheme, std::string face_rec_scheme, std::string limiter, Float toler, int max_iter);

	void run();
	
	void postprocess(std::string outfilename);
};

#endif
