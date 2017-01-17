#include "1deuler.h"

void setup(Grid *const grid, Euler1d *const sim, const size_t num_cells, const int bcleft, const int bcright, const Float leftBVs[NVARS], const Float rightBVs[NVARS], const Float domain_length, 
		const Float _cfl, const char *const _flux)
{
	grid->N = num_cells;
	grid->ncell = N+2;
	grid->nface = N+1;
	
	grid->x = (Float*)malloc((N+2)*sizeof(Float));
	grid->dx = (Float*)malloc((N+2)*sizeof(Float));
	grid->nodes = (Float*)malloc((N+1)*sizeof(Float));

	sim->A = (Float*)malloc((N+2)*sizeof(Float));
	sim->Af = (Float*)malloc((N+1)*sizeof(Float));
	sim->vol = (Float*)malloc((N+2)*sizeof(Float));
	sim->u = (Float**)malloc((N+2)*sizeof(Float*));
	sim->u[0] = (Float*)malloc((N+2)*NVARS*sizeof(Float));
	sim->prim = (Float**)malloc((N+2)*sizeof(Float*));
	sim->prim[0] = (Float*)malloc((N+2)*NVARS*sizeof(Float));
	sim->dudx = (Float**)malloc((N+2)*sizeof(Float*));
	sim->dudx[0] = (Float*)malloc((N+2)*NVARS*sizeof(Float));
	sim->res = (Float**)malloc((N+2)*sizeof(Float*));
	sim->res[0] = (Float*)malloc((N+2)*NVARS*sizeof(Float));

	sim->fluxes = (Float**)malloc((N+1)*sizeof(Float*));
	sim->fluxes[0] = (Float*)malloc(NVARS*(N+1)*sizeof(Float));
	sim->prleft = (Float**)malloc((N+1)*sizeof(Float*));
	sim->prleft[0] = (Float*)malloc((N+1)*NVARS*sizeof(Float));
	sim->prright = (Float**)malloc((N+1)*sizeof(Float*));
	sim->prright[0] = (Float*)malloc((N+1)*NVARS*sizeof(Float));

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

	for(int i = 1; i < N+2; i++)
	{
		sim->u[i] = *sim->u + i*NVARS;
		sim->prim[i] = *sim->prim + i*NVARS;
		sim->dudx[i] = *sim->dudx + i*NVARS;
		sim->res[i] = *sim->res + i*NVARS;
	}
	for(int i = 0; i < N+1; i++)
	{
		fluxes[i] = *fluxes + i*NVARS;
		prleft[i] = *prleft + i*NVARS;
		prright[i] = *prright + i*NVARS;
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
		Float delx = domlen/grid->N;
		grid->x[0] = -delx/2.0;
		grid->dx[0] = delx;
		for(int i = 1; i < grid->N+2; i++)
		{
			dx[i] = delx;
			x[i] = i*delx - delx/2.0;
			if(i < N+1)
				nodes[i] = i*delx;
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
		grid->x[N+1] = grid->nodes[grid->N] + grid->dx[grid->N]/2.0;
		grid->dx[N+1] = grid->dx[N];
	}
}

void set_area(int type, const Float *const cellCenteredAreas, const Grid *const grid, Euler1d *const sim);
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

void compute_inviscid_fluxes_vanleer(const Float *const *const prleft, const Float *const *const prright, const Float *const Af, Float *const *const flux)
{
	#pragma acc kernels present( flux, prleft, prright, Af, res)
	{
		// iterate over interfaces
		#pragma acc loop independent gang worker device_type(nvidia) vector(NVIDIA_VECTOR_LENGTH)
		for(int i = 0; i < N+1; i++)
		{
			compute_vanleerflux_prim(prleft[i], prright[i], flux[i], g);

			// update residual
			for(int j = 0; j < NVARS; j++)
				flux[i][j] *= Af[i];
		}
	}
}

void update_residual(const Float *const *const flux, Float *const *const res)
{
#pragma acc kernels present(flux, res)
	{
#pragma acc loop independent gang worker device_type(nvidia) vector(NVIDIA_VECTOR_LENGTH)
		for(int i = 1; i < N+1; i++)
		{
			for(int j = 0; j < NVARS; j++)
				res[i][j] += flux[i-1][j] - flux[i][j];
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

void compute_source_term(const Float *const *const u, Float *const *const res, const Float *const Af);
{
	#pragma acc kernels present(u, Af, res, g)
	{
		#pragma acc loop independent gang worker device_type(nvidia) vector(NVIDIA_VECTOR_LENGTH)
		for(int i = 1; i <= N; i++)
		{
			Float p = (g-1.0)*(u[i][2] - 0.5*u[i][1]*u[i][1]/u[i][0]);
			res[i][1] += p*(Af[i] - Af[i-1]);
		}
	}
}

void Euler1d::apply_boundary_conditions()
{
	Float** u = this->u;
	Float** prim = this->prim;
	Float* bcvalL = this->bcvalL;
	Float* bcvalR = this->bcvalR;
	int* bcL = &(this->bcL);
	int* bcR = &(this->bcR);
	#pragma acc parallel present(prim[:N+2][:NVARS], u[:N+2][:NVARS], bcvalL[:NVARS], bcvalR[:NVARS], bcL[:1], bcR[:1], g, this) num_gangs(1)
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
			p_in = (g-1.0)*(u[0][2] - 0.5*u[0][0]*v_in*v_in);
			c_in = sqrt(g*p_in/u[0][0]);
			M_in = v_in/c_in;

			if(M_in >= 1.0)
			{
				// supersonic inflow
				// get conserved variables from pt, Tt and M
				//Float astar = 2*g*(g-1.0)/(g+1.0)*Cv*bcvalL[1];
				Float T = bcvalL[1]/(1 + (g-1.0)/2.0*bcvalL[2]*bcvalL[2]);
				Float c = sqrt(g*R*T);
				Float v = bcvalL[2]*c;
				Float p = bcvalL[0]*pow( 1+(g-1.0)/2.0*bcvalL[2]*bcvalL[2], -g/(g-1.0) );
				Float rho = p/(R*T);
				Float E = p/(g-1.0) + 0.5*rho*v*v;
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
				pold = (g-1)*(uold0[2] - 0.5*uold0[1]*uold0[1]/uold0[0]);
				cold = sqrt( g*pold/uold0[0] );
				vold1 = u[1][1]/u[1][0];
				pold1 = (g-1)*(u[1][2] - 0.5*u[1][1]*u[1][1]/u[1][0]);
				cold1 = sqrt( g*pold1/u[1][0] );

				astar = 2*g*(g-1.0)/(g+1.0)*Cv*bcvalL[1];
				dpdu = bcvalL[0]*g/(g-1.0)*pow(1.0-(g-1)/(g+1.0)*vold*vold/astar, 1.0/(g-1.0)) * (-2.0)*(g-1)/(g+1.0)*vold/astar;
				dt0 = cfl*dx[0]/(fabs(vold)+cold);
				lambda = (vold1+vold - cold1-cold)*0.5*dt0/dx[0];
				du = -lambda * (pold1-pold-uold0[0]*cold*(vold1-vold)) / (dpdu-uold0[0]*cold);

				v = vold + du;
				T = bcvalL[1]*(1.0 - (g-1)/(g+1.0)*vold*vold/astar);
				p = bcvalL[0]*pow(T/bcvalL[1], g/(g-1.0));
				u[0][0] = p/(R*T);
				u[0][1] = u[0][0]*v;
				u[0][2] = p/(g-1.0) + 0.5*u[0][0]*v*v;
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
			cold = sqrt(g*pold/uold[0]);
			vold1 = u[N][1]/u[N][0];
			pold1 = (g-1.0)*(u[N][2]-0.5*u[N][0]*vold1*vold1);
			cold1 = sqrt(g*pold1/u[N][0]);

			dt0 = cfl*dx[N+1]/(fabs(vold)+cold);
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

			u[N+1][2] = p/(g-1.0) + 0.5*u[N+1][0]*(vold+dv)*(vold+dv);
			prim[N+1][0] = u[N+1][0];
			prim[N+1][1] = vold+dv;
			prim[N+1][2] = p;
		}
	//else std::cout << "! Euler1D: apply_boundary_conditions(): BC type not recognized!" << std::endl;
	}
}

Euler1dExplicit::Euler1dExplicit(int num_cells, Float length, int leftBCflag, int rightBCflag, std::vector<Float> leftBVs, std::vector<Float> rightBVs,
		Float CFL, std::string inviscid_flux, std::string slope_scheme, std::string face_extrap_scheme, std::string limiter, Float f_time, int temporal_order, std::string RKfile)
	: Euler1d(num_cells,length,leftBCflag,rightBCflag,leftBVs,rightBVs, CFL, inviscid_flux, slope_scheme, face_extrap_scheme, limiter), ftime(f_time), temporalOrder(temporal_order)
{
	RKCoeffs = (Float**)malloc(temporalOrder*sizeof(Float*));
	RKCoeffs[0] = (Float*)malloc(temporalOrder*3*sizeof(Float));
	for(int i = 0; i < temporalOrder; i++)
		RKCoeffs[i] = *RKCoeffs + i*3;

	std::ifstream rkfile(RKfile);
	for(int i = 0; i < temporalOrder; i++)
		for(int j = 0; j < 3; j++)
			rkfile >> RKCoeffs[i][j];
	rkfile.close();

	std::cout << "Euler1dExplicit: Using " << temporalOrder << "-stage TVD RK scheme; loaded coefficients.\n";
}

Euler1dExplicit::~Euler1dExplicit()
{
	free(RKCoeffs[0]);
	free(RKCoeffs);
}

void Euler1dExplicit::run()
{
	int step = 0, istage;
	Float dt = 1.0, time = 0;

	// IC for Sod shock tube
	for(int i = 0; i < N+2; i++)
		if(x[i] <= 0.5)
		{
			u[i][0] = 1.0;
			u[i][1] = 0;
			u[i][2] = 2.5;
			prim[i][0] = 1.0;
			prim[i][1] = 0;
			prim[i][2] = 1.0;
		}
		else
		{
			u[i][0] = 0.125;
			u[i][1] = 0;
			u[i][2] = 0.25;
			prim[i][0] = 0.125;
			prim[i][1] = 0;
			prim[i][2] = 0.1;
		}

	Float* c = (Float*)malloc((N+2)*sizeof(Float));
	Float** uold = (Float**)malloc((N+2)*sizeof(Float*));
	//uold[0] = (Float*)malloc((N+2)*NVARS*sizeof(Float));
	Float** ustage = (Float**)malloc((N+2)*sizeof(Float*));
	ustage[0] = (Float*)malloc((N+2)*NVARS*sizeof(Float));
	for(int i = 0; i < N+2; i++)
	{
		//uold[i] = *uold + i*NVARS;
		ustage[i] = *ustage + i*NVARS;
	}
	
	for(int i = 0; i < ncell; i++)
		uold[i] = (Float*)malloc(NVARS*sizeof(Float));

	Float** u = this->u;
	Float** prim = this->prim;
	Float** RKCoeffs = this->RKCoeffs;
	Float** dudx = this->dudx;
	Float** res = this->res;
	Float** prleft = this->prleft;
	Float** prright = this->prright;
	Float* x = this->x;
	Float* dx = this->dx;
	Float* A = this->A;
	Float* Af = this->Af;
	Float* vol = this->vol;
	Float* nodes = this->nodes;
	Float* bcvalL = this->bcvalL;
	Float* bcvalR = this->bcvalR;
	int* bcL = &(this->bcL);
	int* bcR = &(this->bcR);
	//int* N = &(this->N);
	//Float* cfl = &(this->cfl);
	
	#pragma acc enter data copyin(this)
	#pragma acc enter data copyin(u[:N+2][:NVARS], prim[:N+2][:NVARS], x[:N+2], dx[:N+2], A[:N+2], vol[:N+2], Af[:N+1], nodes[:N+1], RKCoeffs[:temporalOrder][:3], bcvalL[:NVARS], bcvalR[:NVARS], bcL[:1], bcR[:1], g, cfl, dt)
	#pragma acc enter data create(dudx[:N+2][:NVARS], res[:N+2][:NVARS], prleft[:N+1][:NVARS], prright[:N+1][:NVARS], istage, c[:N+2], uold[:N+2][:NVARS], ustage[:N+2][:NVARS])

	while(time < ftime)
	{
		//std::cout << "Euler1dExplicit: run(): Started time loop" << std::endl;

		#pragma acc update self(u[:N+2][:NVARS])
		
		//std::cout << "Euler1dExplicit: run(): Updated self" << std::endl;

		#pragma acc parallel loop present(u[:N+2][:NVARS], uold[:N+2][:NVARS], this) gang worker vector device_type(nvidia) vector_length(NVIDIA_VECTOR_LENGTH)
		for(int i = 0; i < N+2; i++)
		{
			for(int j = 0; j < NVARS; j++)
			{
				uold[i][j] = u[i][j];
			}
		}
		//std::cout << "Euler1dExplicit: run(): Set uold" << std::endl;
		
		// find time step as dt = CFL * min{ dx[i]/(|v[i]|+c[i]) }
		
		#pragma acc parallel loop present(u, c, cfl, dt, g, this) reduction(min:dt) gang worker vector device_type(nvidia) vector_length(NVIDIA_VECTOR_LENGTH)
		for(int i = 1; i < N+1; i++)
		{
			c[i] = cfl*dx[i]/(fabs(u[i][1]) + sqrt(g*(g-1.0) * (u[i][2] - 0.5*u[i][1]*u[i][1]/u[i][0]) / u[i][0]));
			dt = fmin(dt, c[i]);
		}

		/*for(int i = 2; i < N+1; i++)
		{
			a = dx[i]/(fabs(u[i][1]) + c[i]);
			if(a < mws) {
				mws = a;
			}
		}

		dt = cfl*mws;*/

		//std::cout << "Euler1dExplicit: run(): Computed dt" << std::endl;

		// NOTE: moved apply_boundary_conditions() to the top of the inner loop
		for(istage = 0; istage < temporalOrder; istage++)
		{
			// apply BCs
			{
				apply_boundary_conditions();
			}
			
			//std::cout << "Euler1dExplicit: run():  Applied BCs" << std::endl;

			#pragma acc parallel loop present(u, ustage, res, this) gang worker vector device_type(nvidia) vector_length(NVIDIA_VECTOR_LENGTH)
			for(int i = 0; i < N+2; i++)
			{
				for(int j = 0; j < NVARS; j++)
				{
					ustage[i][j] = u[i][j];
					res[i][j] = 0;
				}
			}
			
			//std::cout << "Euler1dExplicit: run():  Set ustage and res" << std::endl;

			cslope->compute_slopes();

			rec->compute_face_values();
			//std::cout << "Euler1dExplicit: run():  Computed face values" << std::endl;

			compute_inviscid_fluxes_cellwise(prleft,prright,res,Af);
			//std::cout << "Euler1dExplicit: run():  Computed fluxes" << std::endl;

			compute_source_term(u,res,Af);
			//std::cout << "Euler1dExplicit: run():  Computed source terms" << std::endl;

			// RK stage
			#pragma acc parallel loop present(prim, u, uold, ustage, res, vol, dt, RKCoeffs, g, this) gang worker vector device_type(nvidia) vector_length(NVIDIA_VECTOR_LENGTH)
			for(int i = 1; i < N+1; i++)
			{
				for(int j = 0; j < NVARS; j++)
					u[i][j] = RKCoeffs[istage][0]*uold[i][j] + RKCoeffs[istage][1]*ustage[i][j] + RKCoeffs[istage][2]*dt/vol[i]*res[i][j];

				prim[i][0] = u[i][0];
				prim[i][1] = u[i][1]/u[i][0];
				prim[i][2] = (g-1.0)*(u[i][2] - 0.5*u[i][1]*prim[i][1]);
			}

		}

		if(step % 10 == 0)
			std::cout << "Euler1dExplicit: run(): Step " << step << " - Time = " << time << std::endl;

#pragma acc update self(dt)
		time += dt;
		step++;
	}

	#pragma update self(u[:N+2][:NVARS], prim[:N+2][:NVARS])
	
	#pragma acc exit data delete(u[:N+2][:NVARS], prim[:N+2][:NVARS], x[:N+2], dx[:N+2], A[:N+2], vol[:N+2], Af[:N+1], nodes[:N+1], RKCoeffs[:temporalOrder][:3], bcvalL[:NVARS], bcvalR[:NVARS], bcL, bcR, g, N, cfl, dudx[:N+2][:NVARS], res[:N+2][:NVARS], prleft[:N+1][:NVARS], prright[:N+1][:NVARS], dt, mws, c[:N+2], uold[:N+2][:NVARS], ustage[:N+2][:NVARS])
	#pragma acc exit data delete(this)

	free(c);
	for(int i = 0; i < ncell; i++)
		free(uold[i]);
	//free(uold[0]);
	free(uold);
	free(ustage[0]);
	free(ustage);

	std::cout << "Euler1dExplicit: run(): Done. Number of time steps = " << step << ", final time = " << time << std::endl;
}

void Euler1dExplicit::postprocess(std::string outfilename)
{
	std::ofstream ofile(outfilename);
	Float pressure, mach, c;
	for(int i = 1; i < N+1; i++)
	{
		pressure = (g-1)*(u[i][2] - 0.5*u[i][1]*u[i][1]/u[i][0]);
		c = sqrt(g*pressure/u[i][0]);
		mach = (u[i][1]/u[i][0])/c;
		ofile << x[i] << " " << u[i][0] << " " << mach << " " << pressure << " " << u[i][1] << " " << c << '\n';
	}
	ofile.close();
}


Euler1dSteadyExplicit::Euler1dSteadyExplicit(int num_cells, Float length, int leftBCflag, int rightBCflag, std::vector<Float> leftBVs, std::vector<Float> rightBVs, 
		Float CFL, std::string inviscidFlux, std::string slope_scheme, std::string face_extrap_scheme, std::string limiter, Float toler, int max_iter)
	: Euler1d(num_cells,length,leftBCflag,rightBCflag,leftBVs,rightBVs, CFL, inviscidFlux, slope_scheme, face_extrap_scheme, limiter), tol(toler), maxiter(max_iter)
{
}

void Euler1dSteadyExplicit::run()
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

void Euler1dSteadyExplicit::postprocess(std::string outfilename)
{
	std::ofstream ofile(outfilename);
	Float pressure, mach, c;
	for(int i = 1; i < N+1; i++)
	{
		pressure = (g-1)*(u[i][2] - 0.5*u[i][1]*u[i][1]/u[i][0]);
		c = sqrt(g*pressure/u[i][0]);
		mach = (u[i][1]/u[i][0])/c;
		ofile << x[i] << " " << u[i][0] << " " << mach << " " << pressure/bcvalL[0] << " " << u[i][1] << " " << c << '\n';
	}
	ofile.close();
}
