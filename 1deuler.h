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

typedef struct
{
	size_t N;							///< Number of real cells in the grid
	size_t ncell;						///< total number of cells including ghost cells
	size_t nface;						///< total number of interfaces to compute fluxes across
	Float* x;							///< Cell centers
	Float* dx;							///< (1D) Size of each cell
	Float* nodes;						///< Mesh nodes
} Grid;

/// Base class for Euler solution processes
typedef struct
{
	Float domlen;					///< Physical length of the domain
	Float* vol;						///< (3D) Volume of each cell
	Float* A;						///< Cross-sectional areas at cell centers
	Float* Af;						///< Cross-sectional areas at interfaces
	Float** u;						///< Conserved variables - u[i][0] is density of cell i and so on
	Float** prim;					///< Primitive variables - u[i][0] is density of cell i and so on
	Float** uleft;					///< Left state of each face
	Float** uright;					///< Right state at each face
	Float** prleft;					///< Left state of each face in terms of primitive variables
	Float** prright;				///< Right state at each face in terms of primitive variables
	Float** dudx;					///< Slope of variables in each cell
	Float** fluxes;					///< Fluxes across the faces
	Float** res;					///< residuals
	int bcL;						///< left BC type
	int bcR;						///< right BC type
	Float bcvalL[NVARS];			///< left boundary value
	Float bcvalR[NVARS];			///< right boundary value
	char* fluxstr;					///< Inviscid flux computation to be used
	char* cslopestr;				///< Slope reconstruction to be used (irrelevant)
	char* recstr;					///< Method for computation of face values of flow variables from their cell-centred values (MUSCL is used always)
	Float cfl;						///< CFL number
	Float g;						///< Adiabatic index
} Euler1d;

/// Setup variables
/** Allocates memory on host for all arrays and copies inputs to variables where needed.
 */
void setup(Grid *const grid, Euler1d *const sim, const size_t num_cells, const int bcleft, const int bcright, const Float bcvalleft[NVARS], const Float bcvalright[NVARS], const Float domain_length,
		const Float _cfl, const char *const _flux);

/// Frees arrays allocated in setup()
void finalize(Grid *const grid, Euler1d *const sim);

/// Generates a grid depending on 
/** \param type If type == 0, a uniform grid is generated and used. If type == 1, then grid points need to passed to the function in
 * \param pointlist an array of positions of mesh points.
 *
 * Note that the domain is assumed to start at x=0.
 */
void generate_mesh(int type, const Float *const pointlist, Grid *const grid);

/// Set cross-sectional areas
void set_area(int type, const Float *const cellCenteredAreas, const Grid *const grid, Euler1d *const sim);

//void compute_slopes();

/// Computes Van Leer flux across all faces into a global flux array
void compute_inviscid_fluxes_vanleer(const Float *const *const prleft, const Float *const *const prright, const Float *const Af, Float *const *const flux);

/// Updates cell residuals from computed fluxes
void update_residual(const Float *const *const flux, Float *const *const res);

//void compute_inviscid_fluxes_cellwise(const Float *const *const prleft, const Float *const *const prright, Float *const *const res, const Float *const Af);

void compute_source_term(const Float *const *const u, Float *const *const res, const Float *const Af);

/// Find new ghost cell values
void apply_boundary_conditions(const Grid* grid, Euler1d* sim);

/// Find new values of left boundary face external state
/** Note that interior states at boundary faces should already be computed.
 */
void apply_boundary_conditions_at_left_boundary(Float *const ul, const Float *const ur);

/// Find new values of right boundary face external state
/** Note that interior states at boundary faces should already be computed.
 */
void apply_boundary_conditions_at_right_boundary(const Float *const ul, Float *const ur);

/// Explicit RK solver for time-dependent 1D Euler equations
typedef struct
{
	Float ftime;								///< Physical time for which to simulate
	Float mws;									///< for computing time steps
	//Float a;									///< temp variable
	int temporalOrder;							///< desired temporal order of accuracy
	Float** RKCoeffs;							///< Low-storage multi-stage TVD RK coefficients
} Euler1dUnsteadyExplicit;

void unsteady_run(const Grid *const grid, Euler1d *const sim, Euler1dUnsteadyExplicit *const tsim);

/// Explicit RK solver for steady-state 1D Euler
typedef struct 
{
	Float tol;
	int maxiter;
	Float mws;									///< for computing time steps
} Euler1dSteadyExplicit;

void steady_run(const Grid *const grid, Euler1d *const sim, Euler1dSteadyExplicit *const tsim);

void postprocess(const Grid *const grid, const Euler1d *const sim, const char *const outfilename);

#endif
