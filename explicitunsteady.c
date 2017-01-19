#include "explicitunsteady.h"

void set_data_unsteady(const size_t num_cells, const int bcleft, const int bcright, const Float bcvalleft[NVARS], const Float bcvalright[NVARS], const Float domain_length, 
		Float f_time, int temporal_order, const char* RKfile, Grid *const grid, Euler1d *const sim, Euler1dUnsteadyExplicit *const tsim)
{
	setup(grid, sim, num_cells, bcleft, bcright, bcvalleft, bcvalright, domain_length, _cfl, _flux);

	RKCoeffs = (Float**)malloc(temporalOrder*sizeof(Float*));
	RKCoeffs[0] = (Float*)malloc(temporalOrder*3*sizeof(Float));
	for(int i = 0; i < temporalOrder; i++)
		RKCoeffs[i] = *RKCoeffs + i*3;

	FILE* rkfile = fopen(RKfile, "r");
	for(int i = 0; i < temporalOrder; i++)
		for(int j = 0; j < 3; j++)
			fscanf(rkfile, "%f", &RKCoeffs[i][j]);
	fclose(rkfile);

	printf("Euler1dExplicit: Using %d-stage TVD RK scheme; loaded coefficients.\n", temporalOrder);
}

void run_unsteady(const Grid *const grid, Euler1d *const sim, Euler1dUnsteadyExplicit *const tsim)
{
	int step = 0, istage;
	Float dt = 1.0, time = 0;

	// IC for Sod shock tube
	for(int i = 0; i < grid->grid->N+2; i++)
		if(grid->x[i] <= 0.5)
		{
			sim->u[i][0] = 1.0;
			sim->u[i][1] = 0;
			sim->u[i][2] = 2.5;
			sim->prim[i][0] = 1.0;
			sim->prim[i][1] = 0;
			sim->prim[i][2] = 1.0;
		}
		else
		{
			sim->u[i][0] = 0.125;
			sim->u[i][1] = 0;
			sim->u[i][2] = 0.25;
			sim->prim[i][0] = 0.125;
			sim->prim[i][1] = 0;
			sim->prim[i][2] = 0.1;
		}

	Float* c = (Float*)malloc((grid->N+2)*sizeof(Float));
	Float** uold = (Float**)malloc((grid->N+2)*sizeof(Float*));
	uold[0] = (Float*)malloc((grid->N+2)*NVARS*sizeof(Float));
	Float** ustage = (Float**)malloc((grid->N+2)*sizeof(Float*));
	ustage[0] = (Float*)malloc((grid->N+2)*NVARS*sizeof(Float));
	for(int i = 0; i < grid->N+2; i++)
	{
		uold[i] = *uold + i*NVARS;
		ustage[i] = *ustage + i*NVARS;
	}
	
	/*for(int i = 0; i < grid->ncell; i++)
		uold[i] = (Float*)malloc(NVARS*sizeof(Float));*/

	Float** u = sim->u;
	Float** prim = sim->prim;
	Float** RKCoeffs = tsim->RKCoeffs;
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

	int N = sim->N;
	Float g = sim->g;
	Float cfl = sim->cfl;
	//int* N = &(grid->N);
	//Float* cfl = &(sim->cfl);
	
	#pragma acc enter data copyin(u[:N+2][:NVARS], prim[:N+2][:NVARS], x[:N+2], dx[:N+2], A[:N+2], vol[:N+2], Af[:N+1], nodes[:N+1], RKCoeffs[:temporalOrder][:3], bcvalL[:NVARS], bcvalR[:NVARS], bcL[:1], bcR[:1], g, cfl, dt)
	#pragma acc enter data create(dudx[:N+2][:NVARS], res[:N+2][:NVARS], fluxes[:N+1][:NVARS], prleft[:N+1][:NVARS], prright[:N+1][:NVARS], istage, c[:N+2], uold[:N+2][:NVARS], ustage[:N+2][:NVARS])

	while(time < tsim->ftime)
	{
		//printf("Euler1dExplicit: run(): Started time loop\n");

		//#pragma acc update self(u[:N+2][:NVARS])
		
		//std::cout << "Euler1dExplicit: run(): Updated self" << std::endl;

		#pragma acc parallel loop present(u[:N+2][:NVARS], uold[:N+2][:NVARS]) gang worker vector device_type(nvidia) vector_length(NVIDIA_VECTOR_LENGTH)
		for(int i = 0; i < N+2; i++)
		{
			for(int j = 0; j < NVARS; j++)
			{
				uold[i][j] = u[i][j];
			}
		}
		//std::cout << "Euler1dExplicit: run(): Set uold" << std::endl;
		
		// find time step as dt = CFL * min{ dx[i]/(|v[i]|+c[i]) }
		
		#pragma acc parallel loop present(u, c, cfl, dt, g) reduction(min:dt) gang worker vector device_type(nvidia) vector_length(NVIDIA_VECTOR_LENGTH)
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

			//cslope->compute_slopes();

			rec->compute_face_values();
			//std::cout << "Euler1dExplicit: run():  Computed face values" << std::endl;

			compute_inviscid_fluxes_vanleer(prleft, prright, Af, fluxes, g);
			//std::cout << "Euler1dExplicit: run():  Computed fluxes" << std::endl;
			
			update_residual(fluxes, res);

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
			printf("Euler1dExplicit: run(): Step %d - Time = %f\n" step, time);

#pragma acc update self(dt)
		time += dt;
		step++;
	}

	#pragma update self(u[:N+2][:NVARS], prim[:N+2][:NVARS])
	
	#pragma acc exit data delete(u[:N+2][:NVARS], prim[:N+2][:NVARS], x[:N+2], dx[:N+2], A[:N+2], vol[:N+2], Af[:N+1], nodes[:N+1], RKCoeffs[:temporalOrder][:3], bcvalL[:NVARS], bcvalR[:NVARS], bcL, bcR, g, N, cfl, dudx[:N+2][:NVARS], res[:N+2][:NVARS], fluxes[:N+1][:NVARS], prleft[:N+1][:NVARS], prright[:N+1][:NVARS], dt, mws, c[:N+2], uold[:N+2][:NVARS], ustage[:N+2][:NVARS])
	#pragma acc exit data delete(this)

	free(c);
	/*for(int i = 0; i < ncell; i++)
		free(uold[i]);*/
	free(uold[0]);
	free(uold);
	free(ustage[0]);
	free(ustage);

	std::cout << "Euler1dExplicit: run(): Done. Number of time steps = " << step << ", final time = " << time << std::endl;
}

void postprocess_unsteady(const Grid *const grid, const Euler1d *const sim, const char *const outfilename)
{
	FILE* ofile = fopen(outfilename, "w");
	Float pressure, mach, c;
	for(int i = 1; i < sim->N+1; i++)
	{
		pressure = (sim->g-1)*(sim->u[i][2] - 0.5*sim->u[i][1]*sim->u[i][1]/sim->u[i][0]);
		c = sqrt(sim->g*pressure/sim->u[i][0]);
		mach = (sim->u[i][1]/sim->u[i][0])/c;
		fprintf(ofile, "%.10e %.10e %.10e %.10e %.10e %.10e\n", grid->x[i], sim->u[i][0], mach, pressure, sim->u[i][1], c);
	}
	fclose(ofile);
}

