#include "explicitunsteady.h"

void setup_data_unsteady(const size_t num_cells, const int bcleft, const int bcright, const Float bcvalleft[NVARS], const Float bcvalright[NVARS], const Float domain_length,
		const char *const _flux, const Float _cfl, Float f_time, int temporal_order, const char* RKfile, Grid *const grid, Euler1d *const sim, Euler1dUnsteadyExplicit *const tsim)
{
	setup(grid, sim, num_cells, bcleft, bcright, bcvalleft, bcvalright, domain_length, _cfl, _flux);
	
	tsim->temporalOrder = temporal_order;
	tsim->ftime = f_time;
	tsim->RKCoeffs = (Float**)malloc(tsim->temporalOrder*sizeof(Float*));
	tsim->RKCoeffs[0] = (Float*)malloc(tsim->temporalOrder*3*sizeof(Float));
	for(int i = 0; i < tsim->temporalOrder; i++)
		tsim->RKCoeffs[i] = *(tsim->RKCoeffs) + i*3;

	FILE* rkfile = fopen(RKfile, "r");
	for(int i = 0; i < tsim->temporalOrder; i++)
		for(int j = 0; j < 3; j++)
			fscanf(rkfile, "%lf", &(tsim->RKCoeffs[i][j]));
	fclose(rkfile);

	printf("Euler1dExplicit: Using %d-stage TVD RK scheme; loaded coefficients.\n", tsim->temporalOrder);
}

void finalize_unsteady(Euler1dUnsteadyExplicit *const tsim)
{
	free(tsim->RKCoeffs[0]);
	free(tsim->RKCoeffs);
}

void run_unsteady(const Grid *const grid, Euler1d *const sim, Euler1dUnsteadyExplicit *const tsim)
{
	int step = 0, istage;
	Float dt = 1.0, time = 0;

	// IC for Sod shock tube
	for(int i = 0; i < grid->N+2; i++)
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

	/*Float** u = sim->u;
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

	int N = grid->N;
	Float g = sim->g;
	Float cfl = sim->cfl;
	Float* pcfl = &(sim->cfl);	// pointer to CFL
	Float* pg = &(sim->g);		// pointer to g
	//int* N = &(grid->N);
	//Float* cfl = &(sim->cfl);
	*/
	
	//printf("Initial: %f, %f, %f\n", sim->u[1][0], sim->u[1][1], sim->u[1][2]);
	
#pragma acc enter data copyin(grid[:1], sim[:1], tsim[:1])
#pragma acc enter data copyin(sim->u[:grid->N+2][:NVARS], sim->prim[:grid->N+2][:NVARS], grid->x[:grid->N+2], grid->dx[:grid->N+2], sim->A[:grid->N+2], sim->vol[:grid->N+2], sim->Af[:grid->N+1], grid->nodes[:grid->N+1])
#pragma acc enter data copyin(sim->bcvalL[:NVARS], sim->bcvalR[:NVARS], dt, tsim->RKCoeffs[:tsim->temporalOrder][:3])
#pragma acc enter data create(sim->dudx[:grid->N+2][:NVARS], sim->res[:grid->N+2][:NVARS], sim->fluxes[:grid->N+1][:NVARS], sim->prleft[:grid->N+1][:NVARS], sim->prright[:grid->N+1][:NVARS], istage, c[:grid->N+2], uold[:grid->N+2][:NVARS], ustage[:grid->N+2][:NVARS])

	while(time < tsim->ftime)
	{
		//printf("Euler1dExplicit: run(): Started time loop\n");

		//#pragma acc update self(u[:N+2][:NVARS])
		
		//std::cout << "Euler1dExplicit: run(): Updated self" << std::endl;

		#pragma acc parallel loop present(sim, uold) gang worker vector device_type(nvidia) vector_length(NVIDIA_VECTOR_LENGTH) //num_workers(NVIDIA_WORKERS_PER_GANG)
		for(int i = 0; i < grid->N+2; i++)
		{
			for(int j = 0; j < NVARS; j++)
			{
				uold[i][j] = sim->u[i][j];
			}
		}
		//std::cout << "Euler1dExplicit: run(): Set uold" << std::endl;
		
		// find time step as dt = CFL * min{ dx[i]/(|v[i]|+c[i]) }
		
		#pragma acc parallel loop present(sim, grid, c, dt) reduction(min:dt) gang worker vector device_type(nvidia) vector_length(NVIDIA_VECTOR_LENGTH) //num_workers(NVIDIA_WORKERS_PER_GANG)
		for(int i = 1; i < grid->N+1; i++)
		{
			c[i] = sim->cfl*grid->dx[i]/(fabs(sim->u[i][1]) + sqrt(sim->g*(sim->g-1.0) * (sim->u[i][2] - 0.5*sim->u[i][1]*sim->u[i][1]/sim->u[i][0]) / sim->u[i][0]));
			dt = fmin(dt, c[i]);
		}

		//std::cout << "Euler1dExplicit: run(): Computed dt" << std::endl;

		// NOTE: moved apply_boundary_conditions() to the top of the inner loop
		for(istage = 0; istage < tsim->temporalOrder; istage++)
		{
			// apply BCs
			apply_boundary_conditions(grid,sim);
			
			//std::cout << "Euler1dExplicit: run():  Applied BCs" << std::endl;

			#pragma acc parallel loop present(sim, ustage) gang worker vector device_type(nvidia) vector_length(NVIDIA_VECTOR_LENGTH) //num_workers(NVIDIA_WORKERS_PER_GANG)
			for(int i = 0; i < grid->N+2; i++)
			{
				for(int j = 0; j < NVARS; j++)
				{
					ustage[i][j] = sim->u[i][j];
					sim->res[i][j] = 0;
				}
			}
			
			//std::cout << "Euler1dExplicit: run():  Set ustage and res" << std::endl;

			//cslope->compute_slopes();

			// compute face values of primitive variables
			//compute_MUSCLReconstruction(N, x, u, prleft, prright, k);
			compute_MUSCLReconstruction(grid, sim);
			//compute_noReconstruction(grid, sim);
			//printf("Euler1dExplicit: run():  Computed face values\n");

			if(sim->fluxid == 0)
				compute_inviscid_fluxes_llf(grid, sim);
			else if(sim->fluxid == 1)
				compute_inviscid_fluxes_vanleer(grid, sim);
			else
				printf("! Euler1dExplicitUnsteady: Invalid flux scheme!\n");
			//std::cout << "Euler1dExplicit: run():  Computed fluxes" << std::endl;
			
			update_residual(grid, sim);

			compute_source_term(grid, sim);
			//std::cout << "Euler1dExplicit: run():  Computed source terms" << std::endl;

			// RK stage
			//#pragma acc parallel loop present(prim, u, uold, ustage, res, vol, dt, RKCoeffs, g) gang worker vector device_type(nvidia) vector_length(NVIDIA_VECTOR_LENGTH)
			#pragma acc parallel loop present(sim, tsim, dt, ustage, uold) gang worker vector device_type(nvidia) vector_length(NVIDIA_VECTOR_LENGTH) //num_workers(NVIDIA_WORKERS_PER_GANG)
			for(int i = 1; i < grid->N+1; i++)
			{
				for(int j = 0; j < NVARS; j++)
					sim->u[i][j] = tsim->RKCoeffs[istage][0]*uold[i][j] + tsim->RKCoeffs[istage][1]*ustage[i][j] + tsim->RKCoeffs[istage][2]*dt/sim->vol[i]*sim->res[i][j];

				sim->prim[i][0] = sim->u[i][0];
				sim->prim[i][1] = sim->u[i][1]/sim->u[i][0];
				sim->prim[i][2] = (sim->g-1.0)*(sim->u[i][2] - 0.5*sim->u[i][1]*sim->prim[i][1]);
			}

		}

		if(step % 10 == 0)
			printf("Euler1dExplicit: run(): Step %d - Time = %f\n", step, time);

#pragma acc update self(dt)
		time += dt;
		step++;
	}

	#pragma update self(u[:grid->N+2][:NVARS], prim[:grid->N+2][:NVARS])
	
#pragma acc exit data delete(sim->u[:grid->N+2][:NVARS], sim->prim[:grid->N+2][:NVARS], grid->x[:grid->N+2], grid->dx[:grid->N+2], sim->A[:grid->N+2], sim->vol[:grid->N+2], sim->Af[:grid->N+1], grid->nodes[:grid->N+1], tsim->RKCoeffs[:tsim->temporalOrder][:3], sim->bcvalL[:NVARS], sim->bcvalR[:NVARS])
#pragma acc exit data delete(sim->dudx[:grid->N+2][:NVARS], sim->res[:grid->N+2][:NVARS], sim->fluxes[:grid->N+1][:NVARS], sim->prleft[:grid->N+1][:NVARS], sim->prright[:grid->N+1][:NVARS],dt,c[:grid->N+2])
#pragma acc exit data delete (uold[:grid->N+2][:NVARS], ustage[:grid->N+2][:NVARS])
#pragma acc exit data delete (sim[:1], tsim[:1], grid[:1])

	free(c);
	/*for(int i = 0; i < ncell; i++)
		free(uold[i]);*/
	free(uold[0]);
	free(uold);
	free(ustage[0]);
	free(ustage);

	printf("Euler1dExplicit: run(): Done. Number of time steps = %d, final time = %f\n", step, time);
}

void postprocess_unsteady(const Grid *const grid, const Euler1d *const sim, const char *const outfilename)
{
	FILE* ofile = fopen(outfilename, "w");
	Float pressure, mach, c;
	for(int i = 1; i < grid->N+1; i++)
	{
		pressure = (sim->g-1)*(sim->u[i][2] - 0.5*sim->u[i][1]*sim->u[i][1]/sim->u[i][0]);
		c = sqrt(sim->g*pressure/sim->u[i][0]);
		mach = (sim->u[i][1]/sim->u[i][0])/c;
		fprintf(ofile, "%.10e %.10e %.10e %.10e %.10e %.10e\n", grid->x[i], sim->u[i][0], mach, pressure, sim->u[i][1], c);
	}
	fclose(ofile);
}

