#include "explicitunsteady.h"
#include "explicitsteady.h"

int main(int argc, char* argv[])
{
	if(argc < 2)
	{
		printf("Euler1d needs a configuration file name. Quitting.\n");
		return -1;
	}

	char* confile = argv[1];
	FILE* conf = fopen(confile, "r");

	int leftbc, rightbc, temporal_order, N, areatype, maxiter;
	Float leftbv[NVARS], rightbv[NVARS], *areas;
	char inv_flux[20], areafile[50], outputfile[20], simtype[20], slope_scheme[20], rec_scheme[20], limiter[20], rkfile[20], dum[50];
	Float cfl, f_time, L, tol;

	int ierr = fscanf(conf, "%s", dum); TASSERT(ierr);
	ierr = fscanf(conf, "%s",simtype); TASSERT(ierr);
	if(strcmp(simtype,"unsteady")==0)
	{
		ierr = fscanf(conf, "%s", dum); ierr = fscanf(conf, "%s",outputfile);
		ierr = fscanf(conf, "%s", dum); ierr = fscanf(conf, "%d",&N);
		ierr = fscanf(conf, "%s", dum); ierr = fscanf(conf, "%lf",&L);
		ierr = fscanf(conf, "%s", dum); ierr = fscanf(conf, "%d",&leftbc);
		ierr = fscanf(conf, "%s", dum); ierr = fscanf(conf, "%d",&rightbc);
		ierr = fscanf(conf, "%s", dum);
		for(int i = 0; i < NVARS; i++)
			ierr = fscanf(conf, "%lf",&leftbv[i]);
		ierr = fscanf(conf, "%s", dum);
		for(int i = 0; i < NVARS; i++)
			ierr = fscanf(conf, "%lf",&rightbv[i]);
		ierr = fscanf(conf, "%s", dum); ierr = fscanf(conf, "%lf",&cfl);
		ierr = fscanf(conf, "%s", dum); ierr = fscanf(conf, "%s",inv_flux);
		ierr = fscanf(conf, "%s", dum); ierr = fscanf(conf, "%s",slope_scheme);
		ierr = fscanf(conf, "%s", dum); ierr = fscanf(conf, "%s",rec_scheme);
		ierr = fscanf(conf, "%s", dum); ierr = fscanf(conf, "%s",limiter);
		ierr = fscanf(conf, "%s", dum); ierr = fscanf(conf, "%lf",&f_time);
		ierr = fscanf(conf, "%s", dum); ierr = fscanf(conf, "%d",&temporal_order);
		ierr = fscanf(conf, "%s", dum); ierr = fscanf(conf, "%s",rkfile);
		ierr = fscanf(conf, "%s", dum); ierr = fscanf(conf, "%d",&areatype);
		if(areatype == 1)
		{
			areas = (Float*)malloc(N*sizeof(Float));
			ierr = fscanf(conf, "%s", dum); ierr = fscanf(conf, "%s",areafile);
			FILE* areaf = fopen(areafile,"r");
			for(int i = 0; i < N; i++)
			{
				ierr = fscanf(areaf, "%lf", &areas[i]);
			}
			fclose(areaf);
		}
		else
		{
			areas = (Float*)malloc(sizeof(Float));
			areas[0] = 1.0;
		}

		TASSERT(ierr);
	
		Float *plist = NULL;

		printf("Unsteady: %d, %f, %d, %d, %s, %f, %f, %d\n", N, L, leftbc, rightbc, inv_flux, cfl, f_time, temporal_order);

		Grid* grid = (Grid*)malloc(sizeof(Grid));
		Euler1dUnsteadyExplicit* tsim = (Euler1dUnsteadyExplicit*)malloc(sizeof(Euler1dUnsteadyExplicit));
		Euler1d* sim = (Euler1d*)malloc(sizeof(Euler1d));

		setup_data_unsteady(N, leftbc, rightbc, leftbv, rightbv, L, inv_flux,
		                    cfl, f_time, temporal_order, rkfile, grid, sim, tsim);
		generate_mesh(0,plist,grid);  // plist is not used here
		set_area(0,areas,grid,sim);
		printf("Setup done.\n");

		run_unsteady(grid, sim, tsim);

		postprocess_unsteady(grid, sim, outputfile);

		finalize(grid,sim);
		finalize_unsteady(tsim);
		free(grid); free(sim); free(tsim);
		free(areas);
	}
	else
	{
		ierr = fscanf(conf, "%s", dum); ierr = fscanf(conf, "%s",outputfile);
		ierr = fscanf(conf, "%s", dum); ierr = fscanf(conf, "%d",&N);
		ierr = fscanf(conf, "%s", dum); ierr = fscanf(conf, "%lf",&L);
		ierr = fscanf(conf, "%s", dum); ierr = fscanf(conf, "%d",&leftbc);
		ierr = fscanf(conf, "%s", dum); ierr = fscanf(conf, "%d",&rightbc);
		ierr = fscanf(conf, "%s", dum);
		for(int i = 0; i < NVARS; i++)
			ierr = fscanf(conf, "%lf",&leftbv[i]);
		ierr = fscanf(conf, "%s", dum);
		for(int i = 0; i < NVARS; i++)
			ierr = fscanf(conf, "%lf",&rightbv[i]);
		ierr = fscanf(conf, "%s", dum); ierr = fscanf(conf, "%lf",&cfl);
		ierr = fscanf(conf, "%s", dum); ierr = fscanf(conf, "%s",inv_flux);
		ierr = fscanf(conf, "%s", dum); ierr = fscanf(conf, "%s",slope_scheme);
		ierr = fscanf(conf, "%s", dum); ierr = fscanf(conf, "%s",rec_scheme);
		ierr = fscanf(conf, "%s", dum); ierr = fscanf(conf, "%s",limiter);
		ierr = fscanf(conf, "%s", dum); ierr = fscanf(conf, "%lf",&tol);
		ierr = fscanf(conf, "%s", dum); ierr = fscanf(conf, "%d",&maxiter);
		ierr = fscanf(conf, "%s", dum); ierr = fscanf(conf, "%d",&areatype);
		if(areatype != 0)
		{
			areas = (Float*)malloc(N*sizeof(Float));
			ierr = fscanf(conf, "%s", dum); ierr = fscanf(conf, "%s",areafile);
			FILE* areaf = fopen(areafile, "r");
			for(int i = 0; i < N; i++)
			{
				ierr = fscanf(areaf, "%lf", &areas[i]);
			}
			fclose(areaf);
		}
		else
		{
			areas = (Float*)malloc(sizeof(Float));
			areas[0] = 1.0;
		}

		TASSERT(ierr);
		
		Float *plist = NULL;

		Grid* grid = (Grid*)malloc(sizeof(Grid));
		Euler1dSteadyExplicit* tsim = (Euler1dSteadyExplicit*)malloc(sizeof(Euler1dUnsteadyExplicit));
		Euler1d* sim = (Euler1d*)malloc(sizeof(Euler1d));

		setup_data_steady(N, leftbc, rightbc, leftbv, rightbv, L, inv_flux, cfl, tol, maxiter, grid, sim, tsim);
		generate_mesh(0,plist,grid);
		set_area(areatype,areas,grid,sim);
		printf("Simulation parameters: %d, %f, %d, %d, %s, %f, %f, %d\n", N, L, leftbc, rightbc, inv_flux, cfl, tol, maxiter);

		run_steady(grid, sim, tsim);

		postprocess_steady(grid, sim, outputfile);
		
		finalize(grid,sim);
		free(grid); free(sim); free(tsim);
		free(areas);
	}

	fclose(conf);

	return 0;
}
