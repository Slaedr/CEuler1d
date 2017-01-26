#include "explicitsteady.h"

void setup_data_steady(const size_t num_cells, const int bcleft, const int bcright, const Float bcvalleft[NVARS], const Float bcvalright[NVARS], const Float domain_length,	const char *const _flux, 
		const Float _cfl, Float _tol, int max_iter, Grid *const grid, Euler1d *const sim, Euler1dSteadyExplicit *const tsim)
{
	tsim->tol = _tol;
	tsim->maxiter = max_iter;
	setup(grid, sim, num_cells, bcleft, bcright, bcvalleft, bcvalright, domain_length, _cfl, _flux);
}

/** Does not converge when MUSCL is used with the Van Albada limiter, for some reason. Reaches a limit cycle too soon.
 * When used without a limiter, MUSCL converges for smooth solutions when a continuous initial condition is used.
 */
void run_steady(const Grid *const grid, Euler1d *const sim, Euler1dSteadyExplicit *const tsim)
{
	int step = 0;
	Float resnorm = 1.0, resnorm0 = 1.0;

	Float* dt = (Float*)malloc((grid->N+2)*sizeof(Float));
	Float** uold = (Float**)malloc((grid->N+2)*sizeof(Float*));
	uold[0] = (Float*)malloc((grid->N+2)*NVARS*sizeof(Float));
	for(int i = 0; i < grid->N+2; i++)
		uold[i] = *uold + i*NVARS;

	// initial conditions
	
	Float p_t = sim->bcvalL[0];
	Float T_t = sim->bcvalL[1];
	Float M = sim->bcvalL[2];
	Float term = 1.0 + (sim->g-1.0)*0.5*M*M;
	Float Tin = T_t/term;
	Float pin = p_t*pow(term, -sim->g/(sim->g-1.0));
	Float cin;

	Float pex = sim->bcvalR[0], tex = sim->bcvalR[1], Mex = sim->bcvalR[2], cex, vex;
	
	// set some cells according to inlet condition
	int pn = grid->N/4;
	for(int i = 0; i <= pn; i++)
	{
		sim->u[i][0] = pin/(GAS_CONSTANT*Tin);
		cin = sqrt(sim->g*pin/sim->u[i][0]);
		sim->u[i][1] = sim->u[i][0]*M*cin;
		sim->u[i][2] = pin/(sim->g-1.0)+0.5*sim->u[i][1]*M*cin;
		sim->prim[i][0] = sim->u[i][0];
		sim->prim[i][1] = M*cin;
		sim->prim[i][2] = pin;
	}

	// set last cells according to exit conditions
	for(int i = grid->N+1; i <= grid->N+1; i++)
	{
		sim->u[i][0] = pex/(GAS_CONSTANT*tex);
		cex = sqrt(sim->g*pex/sim->u[i][0]);
		vex = Mex*cex;
		sim->u[i][1] = sim->u[i][0]*vex;
		sim->u[i][2] = pex/(sim->g-1.0) + 0.5*sim->u[i][0]*vex*vex;
		sim->prim[i][0] = sim->u[i][0];
		sim->prim[i][1] = vex;
		sim->prim[i][2] = pex;
	}

	// linearly interpolate cells in the middle
	for(int i = pn; i < grid->N+1; i++)
	{
		for(int j = 0; j < NVARS; j++)
		{
			sim->u[i][j] = (sim->u[grid->N+1][j]-sim->u[pn][j])/(grid->x[grid->N+1]-grid->x[pn])*(grid->x[i]-grid->x[pn]) + sim->u[pn][j];
			sim->prim[i][j] = (sim->prim[grid->N+1][j]-sim->prim[pn][j])/(grid->x[grid->N+1]-grid->x[pn])*(grid->x[i]-grid->x[pn]) + sim->prim[pn][j];
		}
	}

	/*Float** u = sim->u;
	Float** prim = sim->prim;
	Float** dudx = sim->dudx;
	Float** res = sim->res;
	Float** fluxes = sim->fluxes;
	Float** prleft = sim->prleft;
	Float** prright = sim->prright;
	Float* x = grid->x;
	Float* dx = grid->dx;
	Float* A = sim->A;
	Float* Af = sim->Af;
	Float* vol = sim->vol;
	Float* nodes = grid->nodes;
	Float* bcvalL = sim->bcvalL;
	Float* bcvalR = sim->bcvalR;
	int* bcL = &(sim->bcL);
	int* bcR = &(sim->bcR);

	int N = grid->N;
	Float g = sim->g;
	Float cfl = sim->cfl;
	Float* pcfl = &(sim->cfl);	// pointer to CFL
	Float* pg = &(sim->g);		// pointer to g
	*/

	// Start time loop
	while(resnorm/resnorm0 > tsim->tol && step < tsim->maxiter)
	{
		int i,j;
		for(i = 0; i < grid->N+2; i++)
		{
			for(j = 0; j < NVARS; j++)
			{
				uold[i][j] = sim->u[i][j];
				sim->res[i][j] = 0;
			}
		}
		
		/** NOTE: for some reason, this only works with MUSCL reconstruction WITHOUT a limiter.
		 */
			
		compute_MUSCLReconstruction_unlimited(grid, sim);
		//compute_noReconstruction(grid, sim);
		//printf("Euler1dExplicit: run():  Computed face values\n");

		if(sim->fluxid == 0)
			compute_inviscid_fluxes_llf(grid, sim);
		else if(sim->fluxid == 1)
			compute_inviscid_fluxes_vanleer(grid, sim);
		else
			printf("! Euler1dExplicitSteady: Invalid flux scheme!\n");
		//std::cout << "Euler1dExplicit: run():  Computed fluxes" << std::endl;
		
		update_residual(grid, sim);

		compute_source_term(grid, sim);
		//std::cout << "Euler1dExplicit: run():  Computed source terms" << std::endl;

		// find time step as dt[i] = CFL * dx[i]/(|v[i]|+c[i])
		
		for(i = 1; i < grid->N+1; i++)
		{
			Float c = sqrt( sim->g*(sim->g-1.0) * (sim->u[i][2] - 0.5*sim->u[i][1]*sim->u[i][1]/sim->u[i][0]) / sim->u[i][0] );
			dt[i] = sim->cfl * grid->dx[i]/(fabs(sim->u[i][1]) + c);
		}

		resnorm = 0;
		for(i = 1; i <= grid->N; i++)
			resnorm += sim->res[i][0]*sim->res[i][0]*grid->dx[i];
		resnorm = sqrt(resnorm);
		if(step==0)
			resnorm0 = resnorm;

		// RK step
		for(i = 1; i < grid->N+1; i++)
		{
			for(j = 0; j < NVARS; j++)
				sim->u[i][j] = uold[i][j] + dt[i]/sim->vol[i]*sim->res[i][j];

			sim->prim[i][0] = sim->u[i][0];
			sim->prim[i][1] = sim->u[i][1]/sim->u[i][0];
			sim->prim[i][2] = (sim->g-1.0)*(sim->u[i][2] - 0.5*sim->u[i][1]*sim->prim[i][1]);
		}

		// apply BCs
		apply_boundary_conditions(grid,sim);

		if(step % 10 == 0)
		{
			printf("Euler1dSteadyExplicit: run(): Step %d, relative mass flux norm = %f\n", step, resnorm/resnorm0);
		}

		step++;
	}

	printf("Euler1dExplicit: run(): Done. Number of time steps = %d\n", step);

	if(step == tsim->maxiter)
		printf("Euler1dExplicit: run(): Not converged!\n");

	free(dt);
	free(uold[0]);
	free(uold);
}

void postprocess_steady(const Grid *const grid, const Euler1d *const sim, const char *const outfilename)
{
	printf("postprocess_steady: Writing output to file.\n");
	FILE* ofile = fopen(outfilename, "w");
	Float pressure, mach, c;
	for(int i = 1; i < grid->N+1; i++)
	{
		pressure = (sim->g-1)*(sim->u[i][2] - 0.5*sim->u[i][1]*sim->u[i][1]/sim->u[i][0]);
		c = sqrt(sim->g*pressure/sim->u[i][0]);
		mach = (sim->u[i][1]/sim->u[i][0])/c;
		fprintf(ofile, "%.10e %.10e %.10e %.10e %.10e %.10e\n", grid->x[i], sim->u[i][0], mach, pressure/sim->bcvalL[0], sim->u[i][1], c);
	}
	fclose(ofile);
}
