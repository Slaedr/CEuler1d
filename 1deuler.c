#include "1deuler.h"

void setup(Grid *const grid, Euler1d *const sim, const size_t num_cells, const int bcleft, const int bcright, const Float leftBVs[NVARS], const Float rightBVs[NVARS], const Float domain_length, 
		const Float _cfl, const char *const _flux)
{
	grid->N = num_cells;
	grid->ncell = grid->N+2;
	grid->nface = grid->N+1;
	
	grid->x = (Float*)malloc((grid->N+2)*sizeof(Float));
	grid->dx = (Float*)malloc((grid->N+2)*sizeof(Float));
	grid->nodes = (Float*)malloc((grid->N+1)*sizeof(Float));

	sim->A = (Float*)malloc((grid->N+2)*sizeof(Float));
	sim->Af = (Float*)malloc((grid->N+1)*sizeof(Float));
	sim->vol = (Float*)malloc((grid->N+2)*sizeof(Float));
	sim->u = (Float**)malloc((grid->N+2)*sizeof(Float*));
	sim->u[0] = (Float*)malloc((grid->N+2)*NVARS*sizeof(Float));
	sim->prim = (Float**)malloc((grid->N+2)*sizeof(Float*));
	sim->prim[0] = (Float*)malloc((grid->N+2)*NVARS*sizeof(Float));
	sim->dudx = (Float**)malloc((grid->N+2)*sizeof(Float*));
	sim->dudx[0] = (Float*)malloc((grid->N+2)*NVARS*sizeof(Float));
	sim->res = (Float**)malloc((grid->N+2)*sizeof(Float*));
	sim->res[0] = (Float*)malloc((grid->N+2)*NVARS*sizeof(Float));

	sim->fluxes = (Float**)malloc((grid->N+1)*sizeof(Float*));
	sim->fluxes[0] = (Float*)malloc(NVARS*(grid->N+1)*sizeof(Float));
	sim->prleft = (Float**)malloc((grid->N+1)*sizeof(Float*));
	sim->prleft[0] = (Float*)malloc((grid->N+1)*NVARS*sizeof(Float));
	sim->prright = (Float**)malloc((grid->N+1)*sizeof(Float*));
	sim->prright[0] = (Float*)malloc((grid->N+1)*NVARS*sizeof(Float));

	sim->cfl = _cfl;
	sim->domlen = domain_length;
	sim->bcL = bcleft;
	sim->bcR = bcright;

	for(int i = 0; i < NVARS; i++)
	{
		sim->bcvalL[i] = leftBVs[i];
		sim->bcvalR[i] = rightBVs[i];
	}

	/*for(int i = 0; i < ncell; i++)
		u[i] = (Float*)malloc(NVARS*sizeof(Float));*/

	for(int i = 1; i < grid->N+2; i++)
	{
		sim->u[i] = *sim->u + i*NVARS;
		sim->prim[i] = *sim->prim + i*NVARS;
		sim->dudx[i] = *sim->dudx + i*NVARS;
		sim->res[i] = *sim->res + i*NVARS;
	}
	for(int i = 0; i < grid->N+1; i++)
	{
		sim->fluxes[i] = *(sim->fluxes) + i*NVARS;
		sim->prleft[i] = *(sim->prleft) + i*NVARS;
		sim->prright[i] = *(sim->prright) + i*NVARS;
	}

	sim->fluxstr = (char*)malloc(10*sizeof(char));
	strcpy(sim->fluxstr, _flux);
}

void finalize(Grid *const grid, Euler1d *const sim);
{
	free(sim->fluxstr);
	
	free(grid->x);			
	free(grid->dx);			
	free(grid->nodes);

	free(sim->A);
	free(sim->Af);
	free(sim->vol);		
	free(sim->u[0]);			
	free(sim->prim[0]);		
	free(sim->fluxes[0]);
	free(sim->prleft[0]);	
	free(sim->prright[0]);	
	free(sim->dudx[0]);		
	free(sim->res[0]);

	/*for(int i = 0; i < N+2; i++)
		free(u[i]);*/

	free(sim->u);
	free(sim->prim);		
	free(sim->fluxes);
	free(sim->prleft);	
	free(sim->prright);	
	free(sim->dudx);		
	free(sim->res);
}

void generate_mesh(int type, const Float *const pointlist, Grid *const grid)
{
	if(type == 0)
	{
		grid->nodes[0] = 0.0;
		Float delx = sim->domlen/grid->N;
		grid->x[0] = -delx/2.0;
		grid->dx[0] = delx;
		for(int i = 1; i < grid->N+2; i++)
		{
			grid->dx[i] = delx;
			grid->x[i] = i*delx - delx/2.0;
			if(i < grid->N+1)
				grid->nodes[i] = i*delx;
		}
	}
	else
	{
		for(int i = 0; i < grid->N+1; i++)
			grid->nodes[i] = pointlist[i];
		for(int i = 1; i < grid->N+1; i++)
		{
			grid->x[i] = (grid->nodes[i]+grid->nodes[i-1])/2.0;
			grid->dx[i] = grid->nodes[i]-grid->nodes[i-1];
		}

		// ghost cells
		grid->x[0] = -grid->x[1]; 
		grid->dx[0] = grid->dx[1];
		grid->x[grid->N+1] = grid->nodes[grid->N] + grid->dx[grid->N]/2.0;
		grid->dx[grid->N+1] = grid->dx[grid->N];
	}
}

void set_area(int type, const Float *const cellCenteredAreas, const Grid *const grid, Euler1d *const sim)
{
	if(type == 0)
		for(int i = 0; i < grid->N+2; i++)
			sim->A[i] = cellCenteredAreas[0];
	else
	{
		for(int i = 1; i < grid->N+1; i++)
			sim->A[i] = cellCenteredAreas[i-1];
		/*Float h = 0.15, t1 = 0.8, t2 = 3.0;
		for(int i = 1; i < N+1; i++)
			A[i] = 1.0 - h*pow(sin(PI*pow(x[i],t1)),t2);*/

		// maybe assign ghost cell areas by linear extrapolation?
		sim->A[0] = sim->A[1];
		sim->A[N+1] = sim->A[grid->N];
	}

	/** Get interface areas as inverse-distance weighted averages of cell-centered areas so that they are exact for linear profiles.
	 */
	for(int i = 0; i <= grid->N; i++)
		sim->Af[i] = (sim->A[i]*grid->dx[i+1] + sim->A[i+1]*grid->dx[i])/(grid->dx[i]+grid->dx[i+1]);

	for(int i = 0; i < grid->N+2; i++)
	{
		sim->vol[i] = grid->dx[i]*sim->A[i];
	}
}

void compute_inviscid_fluxes_vanleer(const Grid *const grid, Euler1d *const sim)
{
	#pragma acc kernels present( grid, sim) 
	// sim->fluxes, sim->prleft, sim->prright, sim->Af, sim->res, sim->g)
	{
		// iterate over interfaces
		#pragma acc loop independent gang worker device_type(nvidia) vector(NVIDIA_VECTOR_LENGTH)
		for(int i = 0; i < grid->N+1; i++)
		{
			compute_vanleerflux_prim(sim->prleft[i], sim->prright[i], sim->fluxes[i], sim->g);

			// update residual
			for(int j = 0; j < NVARS; j++)
				sim->fluxes[i][j] *= sim->Af[i];
		}
	}
}

void compute_inviscid_fluxes_llf(const Grid *const grid, Euler1d *const sim)
{
	#pragma acc kernels present(grid, sim) 
	//fluxes, prleft, prright, Af, res, g)
	{
		// iterate over interfaces
		#pragma acc loop independent gang worker device_type(nvidia) vector(NVIDIA_VECTOR_LENGTH)
		for(int i = 0; i < grid->N+1; i++)
		{
			compute_llfflux_prim(sim->prleft[i], sim->prright[i], sim->fluxes[i], sim->g);

			// update residual
			for(int j = 0; j < NVARS; j++)
				sim->fluxes[i][j] *= sim->Af[i];
		}
	}
}

void update_residual(const Grid *const grid, Euler1d *const sim)
{
#pragma acc kernels present(grid, sim)
	{
#pragma acc loop independent gang worker device_type(nvidia) vector(NVIDIA_VECTOR_LENGTH)
		for(int i = 1; i < grid->N+1; i++)
		{
			for(int j = 0; j < NVARS; j++)
				sim->res[i][j] += sim->flux[i-1][j] - sim->flux[i][j];
		}
	}
}

/*void compute_inviscid_fluxes_cellwise(const Float *const *const prleft, const Float *const *const prright, Float *const *const res, const Float *const Af, const Grid *const gr)
{
	// NOTE: Do we really need to allocate this much?
	Float** fluxes = (Float**)malloc((gr->nface)*sizeof(Float*));
	fluxes[0] = (Float*)malloc(NVARS*gr->nface*sizeof(Float));
	for(int i = 0; i < gr->nface; i++)
		fluxes[i] = *fluxes + i*NVARS;

	#pragma acc kernels present( prleft, prright, Af, res, g, gr[0:1]) create(fluxes[:N+1][:NVARS])
	{
		// iterate over interfaces
		#pragma acc loop independent gang worker device_type(nvidia) vector(NVIDIA_VECTOR_LENGTH)
		for(int i = 1; i < gr->nface; i++)
		{
			compute_vanleerflux_prim(prleft[i-1], prright[i-1], fluxes[i], g);

			for(int j = 0; j < NVARS; j++)
			{
				fluxes[i][j] *= Af[i-1];
				res[i][j] += fluxes[i][j];
			}
			
			compute_vanleerflux_prim(prleft[i], prright[i], fluxes[i], g);

			for(int j = 0; j < NVARS; j++)
			{
				fluxes[i][j] *= Af[i];
				res[i][j] -= fluxes[i][j];
			}
		}
	}
	
	free(fluxes[0]);
	free(fluxes);
}*/

void compute_source_term(const Grid *const grid, Euler1d *const sim)
{
	#pragma acc kernels present(grid, sim)
	{
		#pragma acc loop independent gang worker device_type(nvidia) vector(NVIDIA_VECTOR_LENGTH)
		for(int i = 1; i <= grid->N; i++)
		{
			Float p = (sim->g-1.0)*(sim->u[i][2] - 0.5*sim->u[i][1]*sim->u[i][1]/sim->u[i][0]);
			sim->res[i][1] += p*(sim->Af[i] - sim->Af[i-1]);
		}
	}
}

void apply_boundary_conditions(Euler1d *const sim)
{
	Float** u = sim->u;
	Float** prim = sim->prim;
	Float* bcvalL = sim->bcvalL;
	Float* bcvalR = sim->bcvalR;
	int* bcL = &(sim->bcL);
	int* bcR = &(sim->bcR);
	#pragma acc parallel present(sim, prim[:N+2][:NVARS], u[:N+2][:NVARS], bcvalL[:NVARS], bcvalR[:NVARS], bcL[:1], bcR[:1]) num_gangs(1)
	{
		if(*bcL == 0)
		{
			// homogeneous Dirichlet - wall
			u[0][0] = u[1][0];
			u[0][1] = -u[1][1];
			u[0][2] = u[1][2];
			prim[0][0] = prim[1][0];
			prim[0][1] = -prim[1][1];
			prim[0][2] = prim[1][2];
		}
		else if(*bcL == 1)
		{
			Float M_in, c_in, v_in, p_in;
			v_in = u[0][1]/u[0][0];
			p_in = (sim->g-1.0)*(u[0][2] - 0.5*u[0][0]*v_in*v_in);
			c_in = sqrt(sim->g*p_in/u[0][0]);
			M_in = v_in/c_in;

			if(M_in >= 1.0)
			{
				// supersonic inflow
				// get conserved variables from pt, Tt and M
				//Float astar = 2*g*(g-1.0)/(g+1.0)*Cv*bcvalL[1];
				Float T = bcvalL[1]/(1 + (sim->g-1.0)/2.0*bcvalL[2]*bcvalL[2]);
				Float c = sqrt(g*R*T);
				Float v = bcvalL[2]*c;
				Float p = bcvalL[0]*pow( 1+(sim->g-1.0)/2.0*bcvalL[2]*bcvalL[2], -sim->g/(sim->g-1.0) );
				Float rho = p/(R*T);
				Float E = p/(sim->g-1.0) + 0.5*rho*v*v;
				/*u[0][0] = 2*rho - u[1][0];
				u[0][1] = 2*rho*v - u[1][1];
				u[0][2] = 2*E - u[1][2];*/
				u[0][0] = rho;
				u[0][1] = rho*v;
				u[0][2] = E;
				prim[0][0] = rho;
				prim[0][1] = v;
				prim[0][2] = p;
			}
			else if(M_in >= 0)
			{
				// subsonic inflow
				// get conserved variables from pt and Tt specified in bcvalL[0] and bcvalL[1] respectively
				Float vold, pold, cold, vold1, pold1, cold1, astar, dpdu, dt0, lambda, du, v, T, p;
				Float* uold0 = u[0];
				vold = uold0[1]/uold0[0];
				pold = (sim->g-1)*(uold0[2] - 0.5*uold0[1]*uold0[1]/uold0[0]);
				cold = sqrt( sim->g*pold/uold0[0] );
				vold1 = u[1][1]/u[1][0];
				pold1 = (sim->g-1)*(u[1][2] - 0.5*u[1][1]*u[1][1]/u[1][0]);
				cold1 = sqrt( sim->g*pold1/u[1][0] );

				astar = 2*sim->g*(sim->g-1.0)/(sim->g+1.0)*Cv*bcvalL[1];
				dpdu = bcvalL[0]*sim->g/(sim->g-1.0)*pow(1.0-(sim->g-1)/(sim->g+1.0)*vold*vold/astar, 1.0/(sim->g-1.0)) * (-2.0)*(sim->g-1)/(sim->g+1.0)*vold/astar;
				dt0 = sim->cfl*dx[0]/(fabs(vold)+cold);
				lambda = (vold1+vold - cold1-cold)*0.5*dt0/dx[0];
				du = -lambda * (pold1-pold-uold0[0]*cold*(vold1-vold)) / (dpdu-uold0[0]*cold);

				v = vold + du;
				T = bcvalL[1]*(1.0 - (sim->g-1)/(sim->g+1.0)*vold*vold/astar);
				p = bcvalL[0]*pow(T/bcvalL[1],sim-> g/(sim->g-1.0));
				u[0][0] = p/(R*T);
				u[0][1] = u[0][0]*v;
				u[0][2] = p/(sim->g-1.0) + 0.5*u[0][0]*v*v;
				prim[0][0] = u[0][0];
				prim[0][1] = v;
				prim[0][2] = p;

				/*c = sqrt(g*p/u[0][0]);
				M = v/c;
				std::cout << "  apply_boundary_conditions(): Inlet ghost cell mach number = " << M << std::endl;*/
			}
#ifndef _OPENACC
			else
				std::cout << "! Euler1d: apply_boundary_conditions(): Error! Inlet is becoming outlet!" << std::endl;
#endif
		}
		//else std::cout << "! Euler1D: apply_boundary_conditions(): BC type not recognized!" << std::endl;

		if(*bcR == 0)
		{
			u[N+1][0] = u[N][0];
			u[N+1][1] = -u[N][1];
			u[N+1][2] = u[N][2];
			prim[N+1][0] = prim[N][0];
			prim[N+1][1] = -prim[N][1];
			prim[N+1][2] = prim[N][2];
		}
		else if(*bcR == 3)
		{
			// outflow
			Float l1, l2, l3, cold, cold1, pold, pold1, vold, vold1, dt0, r1, r2, r3, Mold, dp, drho, dv, p;
			Float* uold = u[N+1];
			vold = uold[1]/uold[0];
			pold = (g-1.0)*(uold[2]-0.5*uold[0]*vold*vold);
			cold = sqrt(sim->g*pold/uold[0]);
			vold1 = u[N][1]/u[N][0];
			pold1 = (sim->g-1.0)*(u[N][2]-0.5*u[N][0]*vold1*vold1);
			cold1 = sqrt(sim->g*pold1/u[N][0]);

			dt0 = sim->cfl*dx[N+1]/(fabs(vold)+cold);
			l1 = (vold+vold1)*0.5*dt0/dx[N+1];
			l2 = (vold+vold1 + cold+cold1)*0.5*dt0/dx[N+1];
			l3 = (vold+vold1 - cold-cold1)*0.5*dt0/dx[N+1];

			r1 = -l1*( u[N+1][0] - u[N][0] - 1.0/(cold*cold)*(pold - pold1));
			r2 = -l2*( pold - pold1 + u[N+1][0]*cold*(vold - vold1));
			r3 = -l3*( pold - pold1 - u[N+1][0]*cold*(vold - vold1));
			Mold = (vold+vold1)/(cold+cold1);

			// check whether supersonic or subsonic
			if(Mold > 1)
				dp = 0.5*(r2+r3);
			else
				dp = 0;

			drho = r1 + dp/(cold*cold);
			dv = (r2-dp)/(u[N+1][0]*cold);

			u[N+1][0] += drho;
			u[N+1][1] = u[N+1][0]*(vold + dv);

			if(Mold > 1)
				p = pold + dp;
			else
				p = bcvalR[0];

			u[N+1][2] = p/(sim->g-1.0) + 0.5*u[N+1][0]*(vold+dv)*(vold+dv);
			prim[N+1][0] = u[N+1][0];
			prim[N+1][1] = vold+dv;
			prim[N+1][2] = p;
		}
	//else std::cout << "! Euler1D: apply_boundary_conditions(): BC type not recognized!" << std::endl;
	}
}

