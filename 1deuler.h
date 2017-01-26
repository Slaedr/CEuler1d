/** \file 1deuler.h
 * \brief Code to control the solution process for 1D Euler equations.
 * \author Aditya Kashi
 * \date Jan 2017
 */

#ifndef __1DEULER_H
#define __1DEULER_H

#ifndef __DEFINITIONS_H
#include "definitions.h"
#endif

#ifndef __INVISCIDFLUX_H
#include "inviscidflux.h"
#endif

typedef struct
{
	size_t N;							///< Number of real cells in the grid
	size_t ncell;						///< total number of cells including ghost cells
	size_t nface;						///< total number of interfaces to compute fluxes across
	Float* x;							///< Cell centers
	Float* dx;							///< (1D) Size of each cell
	Float* nodes;						///< Mesh nodes
	Float domlen;						///< Physical length of the domain
} Grid;

/// Base struct for Euler solution processes
typedef struct
{
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
	int fluxid;						///< Integer denoting which flux to use: 0 -> LLF, 1 -> Van Leer
	char* cslopestr;				///< Slope reconstruction to be used (irrelevant)
	char* recstr;					///< Method for computation of face values of flow variables from their cell-centred values (MUSCL is used always)
	Float cfl;						///< CFL number
	Float Cv;						///< Specific heat at constant volume
	Float g;						///< Adiabatic index
	Float muscl_k;					///< MUSCL reconstruction factor (ideally, 1/3); currently set statically in setup()
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
void compute_inviscid_fluxes_vanleer(const Grid *const grid, Euler1d *const sim);

/// Computes LLF fluxes across all faces into a global flux array
void compute_inviscid_fluxes_llf(const Grid *const grid, Euler1d *const sim);

/// Updates cell residuals from computed fluxes
void update_residual(const Grid *const grid, Euler1d *const sim);

//void compute_inviscid_fluxes_cellwise(const Float *const *const prleft, const Float *const *const prright, Float *const *const res, const Float *const Af);

/// Adds contribution of source terms to the residual
void compute_source_term(const Grid *const grid, Euler1d *const sim);

/// Find new ghost cell values
void apply_boundary_conditions(const Grid *const grid, Euler1d *const sim);

#endif
