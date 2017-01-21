#include "explicitsteady.h"

Euler1dSteadyExplicit::Euler1dSteadyExplicit(int num_cells, Float length, int leftBCflag, int rightBCflag, std::vector<Float> leftBVs, std::vector<Float> rightBVs, 
		Float CFL, std::string inviscidFlux, std::string slope_scheme, std::string face_extrap_scheme, std::string limiter, Float toler, int max_iter)
	: Euler1d(num_cells,length,leftBCflag,rightBCflag,leftBVs,rightBVs, CFL, inviscidFlux, slope_scheme, face_extrap_scheme, limiter), tol(toler), maxiter(max_iter)
{
}

void setup_data_steady(const size_t num_cells, const int bcleft, const int bcright, const Float bcvalleft[NVARS], const Float bcvalright[NVARS], const Float domain_length,	const char *const _flux, 
		const Float _cfl, Float _tol, int max_iter, Grid *const grid, Euler1d *const sim, Euler1dSteadyExplicit *const tsim)
{
	tsim->tol = _tol;
	tsim->maxiter = max_iter;
	set_data(grid, sim, num_cells, bcleft, bcright, bcvalleft, bcvalright, domain_length, _cfl, _flux);
}

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
	Float term = 1.0 + (g-1.0)*0.5*M*M;
	Float Tin = T_t/term;
	Float pin = p_t*pow(term, -sim->g/(sim->g-1.0));
	Float cin;

	Float pex = sim->bcvalR[0], tex = sim->bcvalR[1], Mex = sim->bcvalR[2], cex, vex;
	
	// set some cells according to inlet condition
	int pn = grid->N/4;
	for(int i = 0; i <= pn; i++)
	{
		sim->u[i][0] = pin/(R*Tin);
		cin = sqrt(g*pin/sim->u[i][0]);
		sim->u[i][1] = sim->u[i][0]*M*cin;
		sim->u[i][2] = pin/(sim->g-1.0)+0.5*sim->u[i][1]*M*cin;
		sim->prim[i][0] = sim->u[i][0];
		sim->prim[i][1] = M*cin;
		sim->prim[i][2] = pin;
	}

	// set last cells according to exit conditions
	for(int i = grid->N+1; i <= grid->N+1; i++)
	{
		sim->u[i][0] = pex/(R*tex);
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

	Float** u = sim->u;
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
	
	// Start time loop
	while(resnorm/resnorm0 > tol && step < maxiter)
	{
		int i,j;
		for(i = 0; i < N+2; i++)
		{
			for(j = 0; j < NVARS; j++)
			{
				uold[i][j] = u[i][j];
				res[i][j] = 0;
			}
		}
			
		compute_MUSCLReconstruction(N, x, u, prleft, prright, k);
		//printf("Euler1dExplicit: run():  Computed face values\n");

		compute_inviscid_fluxes_vanleer(prleft, prright, Af, fluxes, g);
		//std::cout << "Euler1dExplicit: run():  Computed fluxes" << std::endl;
		
		update_residual(fluxes, res);

		compute_source_term(u,res,Af);
		//std::cout << "Euler1dExplicit: run():  Computed source terms" << std::endl;

		// find time step as dt = CFL * min{ dx[i]/(|v[i]|+c[i]) }
		
		for(i = 1; i < N+1; i++)
		{
			Float c = sqrt( g*(g-1.0) * (u[i][2] - 0.5*u[i][1]*u[i][1]/u[i][0]) / u[i][0] );
			dt[i] = cfl * dx[i]/(fabs(u[i][1]) + c);
		}

		resnorm = 0;
		for(i = 1; i <= N; i++)
			resnorm += res[i][0]*res[i][0]*dx[i];
		resnorm = sqrt(resnorm);
		if(step==0)
			resnorm0 = resnorm;

		// RK step
		for(i = 1; i < N+1; i++)
		{
			for(j = 0; j < NVARS; j++)
				u[i][j] = uold[i][j] + dt[i]/vol[i]*res[i][j];

			prim[i][0] = u[i][0];
			prim[i][1] = u[i][1]/u[i][0];
			prim[i][2] = (g-1.0)*(u[i][2] - 0.5*u[i][1]*prim[i][1]);
		}

		// apply BCs
		apply_boundary_conditions();

		if(step % 10 == 0)
		{
			printf("Euler1dSteadyExplicit: run(): Step %d, relative mass flux norm = %f\n", step, resnorm/resnorm0);
		}

		step++;
	}

	printf("Euler1dExplicit: run(): Done. Number of time steps = %d\n", step);

	if(step == maxiter)
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
	for(int i = 1; i < sim->N+1; i++)
	{
		pressure = (sim->g-1)*(sim->u[i][2] - 0.5*sim->u[i][1]*sim->u[i][1]/sim->u[i][0]);
		c = sqrt(sim->g*pressure/sim->u[i][0]);
		mach = (sim->u[i][1]/sim->u[i][0])/c;
		fprintf(ofile, "%.10e %.10e %.10e %.10e %.10e %.10e\n", grid->x[i], sim->u[i][0], mach, pressure/sim->bcvalL[0], sim->u[i][1], c);
	}
	fclose(ofile);
}
