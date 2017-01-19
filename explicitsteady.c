#include "explicitsteady.h"

Euler1dSteadyExplicit::Euler1dSteadyExplicit(int num_cells, Float length, int leftBCflag, int rightBCflag, std::vector<Float> leftBVs, std::vector<Float> rightBVs, 
		Float CFL, std::string inviscidFlux, std::string slope_scheme, std::string face_extrap_scheme, std::string limiter, Float toler, int max_iter)
	: Euler1d(num_cells,length,leftBCflag,rightBCflag,leftBVs,rightBVs, CFL, inviscidFlux, slope_scheme, face_extrap_scheme, limiter), tol(toler), maxiter(max_iter)
{
}

void run_steady()
{
	int step = 0;
	Float resnorm = 1.0, resnorm0 = 1.0;
	std::vector<Float> dt(N+2);

	std::vector<Float> c(N+2);
	std::vector<std::vector<Float>> uold;
	uold.resize(N+2);
	for(int i = 0; i < N+2; i++)
		uold[i].resize(NVARS);

	// initial conditions
	
	Float p_t = bcvalL[0];
	Float T_t = bcvalL[1];
	Float M = bcvalL[2];
	Float term = 1.0 + (g-1.0)*0.5*M*M;
	Float Tin = T_t/term;
	Float pin = p_t*pow(term, -g/(g-1.0));
	Float cin;

	Float pex = bcvalR[0], tex = bcvalR[1], Mex = bcvalR[2], cex, vex;
	
	// set some cells according to inlet condition
	int pn = N/4;
	for(int i = 0; i <= pn; i++)
	{
		u[i][0] = pin/(R*Tin);
		cin = sqrt(g*pin/u[i][0]);
		u[i][1] = u[i][0]*M*cin;
		u[i][2] = pin/(g-1.0)+0.5*u[i][1]*M*cin;
		prim[i][0] = u[i][0];
		prim[i][1] = M*cin;
		prim[i][2] = pin;
	}

	// set last cells according to exit conditions
	for(int i = N+1; i <= N+1; i++)
	{
		u[i][0] = pex/(R*tex);
		cex = sqrt(g*pex/u[i][0]);
		vex = Mex*cex;
		u[i][1] = u[i][0]*vex;
		u[i][2] = pex/(g-1.0) + 0.5*u[i][0]*vex*vex;
		prim[i][0] = u[i][0];
		prim[i][1] = vex;
		prim[i][2] = pex;
	}

	// linearly interpolate cells in the middle
	for(int i = pn; i < N+1; i++)
	{
		for(int j = 0; j < NVARS; j++)
		{
			u[i][j] = (u[N+1][j]-u[pn][j])/(x[N+1]-x[pn])*(x[i]-x[pn]) + u[pn][j];
			prim[i][j] = (prim[N+1][j]-prim[pn][j])/(x[N+1]-x[pn])*(x[i]-x[pn]) + prim[pn][j];
		}
	}

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

		//cslope->compute_slopes();
		rec->compute_face_values();

		compute_inviscid_fluxes(prleft,prright,res,Af);
		compute_source_term(u,res,Af);

		// find time step as dt = CFL * min{ dx[i]/(|v[i]|+c[i]) }
		
		for(i = 1; i < N+1; i++)
		{
			c[i] = sqrt( g*(g-1.0) * (u[i][2] - 0.5*u[i][1]*u[i][1]/u[i][0]) / u[i][0] );
		}

		for(i = 1; i < N+1; i++)
		{
			dt[i] = cfl * dx[i]/(fabs(u[i][1]) + c[i]);
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
			std::cout << "Euler1dSteadyExplicit: run(): Step " << step << ", relative mass flux norm = " << resnorm/resnorm0 << std::endl;
		}

		step++;
	}

	std::cout << "Euler1dExplicit: run(): Done. Number of time steps = " << step << std::endl;

	/*for(int j = 0; j < NVARS; j++)
	{
		for(int i = 0; i <= N+1; i++)
			std::cout << dudx[i][j] << " ";
		std::cout << std::endl;
	}
	std::cout << std::endl;*/

	if(step == maxiter)
		std::cout << "Euler1dExplicit: run(): Not converged!" << std::endl;
}

void postprocess_steady(const Grid *const grid, const Euler1d *const sim, const char *const outfilename)
{
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
