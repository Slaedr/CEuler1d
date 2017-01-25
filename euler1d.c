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
	printf("%s\n", confile);
	FILE* conf = fopen(confile, "r");

	int leftbc, rightbc, temporal_order, N, areatype, maxiter;
	Float leftbv[NVARS], rightbv[NVARS], *areas;
	char inv_flux[20], areafile[50], outputfile[20], simtype[20], slope_scheme[20], rec_scheme[20], limiter[20], rkfile[20], dum[50];
	Float cfl, f_time, L, tol;

	fscanf(conf, "%s", dum);
	fscanf(conf, "%s",simtype);
	if(strcmp(simtype,"unsteady"))
	{
		fscanf(conf, "%s", dum); fscanf(conf, "%s",outputfile);
		fscanf(conf, "%s", dum); fscanf(conf, "%d",&N);
		fscanf(conf, "%s", dum); fscanf(conf, "%f",&L);
		fscanf(conf, "%s", dum); fscanf(conf, "%d",&leftbc);
		fscanf(conf, "%s", dum); fscanf(conf, "%d",&rightbc);
		fscanf(conf, "%s", dum);
		for(int i = 0; i < NVARS; i++)
			fscanf(conf, "%f",&leftbv[i]);
		fscanf(conf, "%s", dum);
		for(int i = 0; i < NVARS; i++)
			fscanf(conf, "%f",&rightbv[i]);
		fscanf(conf, "%s", dum); fscanf(conf, "%f",&cfl);
		fscanf(conf, "%s", dum); fscanf(conf, "%s",inv_flux);
		fscanf(conf, "%s", dum); fscanf(conf, "%s",slope_scheme);
		fscanf(conf, "%s", dum); fscanf(conf, "%s",rec_scheme);
		fscanf(conf, "%s", dum); fscanf(conf, "%s",limiter);
		fscanf(conf, "%s", dum); fscanf(conf, "%f",&f_time);
		fscanf(conf, "%s", dum); fscanf(conf, "%d",&temporal_order);
		fscanf(conf, "%s", dum); fscanf(conf, "%s",rkfile);
		fscanf(conf, "%s", dum); fscanf(conf, "%d",&areatype);
		if(areatype == 1)
		{
			areas = (Float*)malloc(N*sizeof(Float));
			fscanf(conf, "%s", dum); fscanf(conf, "%s",areafile);
			FILE* areaf = fopen(areafile,"r");
			for(int i = 0; i < N; i++)
			{
				fscanf(areaf, "%f", &areas[i]);
			}
			fclose(areaf);
		}
		else
		{
			areas = (Float*)malloc(sizeof(Float));
			areas[0] = 1.0;
		}
	
		Float *plist;

		printf("Unsteady: %d, %f, %d, %d, %s, %f, %f, %d\n", N, L, leftbc, rightbc, inv_flux, cfl, f_time, temporal_order);

		Grid* grid = (Grid*)malloc(sizeof(Grid));
		Euler1dUnsteadyExplicit* tsim = (Euler1dUnsteadyExplicit*)malloc(sizeof(Euler1dUnsteadyExplicit));
		Euler1d* sim = (Euler1d*)malloc(sizeof(Euler1d));

		setup_data_unsteady(N, leftbc, rightbc, leftbv, rightbv, L, inv_flux, cfl, f_time, temporal_order, rkfile, grid, sim, tsim);
		generate_mesh(0,plist,grid);
		set_area(0,areas,grid,sim);

		run_unsteady(grid, sim, tsim);

		postprocess_unsteady(grid, sim, outputfile);

		finalize(grid,sim);
		finalize_unsteady(tsim);
		free(grid); free(sim); free(tsim);
		free(areas);
	}
	else
	{
		fscanf(conf, "%s", dum); fscanf(conf, "%s",outputfile);
		fscanf(conf, "%s", dum); fscanf(conf, "%d",&N);
		fscanf(conf, "%s", dum); fscanf(conf, "%f",&L);
		fscanf(conf, "%s", dum); fscanf(conf, "%d",&leftbc);
		fscanf(conf, "%s", dum); fscanf(conf, "%d",&rightbc);
		fscanf(conf, "%s", dum);
		for(int i = 0; i < NVARS; i++)
			fscanf(conf, "%f",&leftbv[i]);
		fscanf(conf, "%s", dum);
		for(int i = 0; i < NVARS; i++)
			fscanf(conf, "%f",&rightbv[i]);
		fscanf(conf, "%s", dum); fscanf(conf, "%f",&cfl);
		fscanf(conf, "%s", dum); fscanf(conf, "%s",inv_flux);
		fscanf(conf, "%s", dum); fscanf(conf, "%s",slope_scheme);
		fscanf(conf, "%s", dum); fscanf(conf, "%s",rec_scheme);
		fscanf(conf, "%s", dum); fscanf(conf, "%s",limiter);
		fscanf(conf, "%s", dum); fscanf(conf, "%f",&tol);
		fscanf(conf, "%s", dum); fscanf(conf, "%d",&maxiter);
		fscanf(conf, "%s", dum); fscanf(conf, "%d",&areatype);
		if(areatype != 0)
		{
			areas = (Float*)malloc(N*sizeof(Float));
			fscanf(conf, "%s", dum); fscanf(conf, "%s",areafile);
			FILE* areaf = fopen(areafile, "r");
			for(int i = 0; i < N; i++)
			{
				fscanf(areaf, "%f", &areas[i]);
			}
			fclose(areaf);
		}
		else
		{
			areas = (Float*)malloc(sizeof(Float));
			areas[0] = 1.0;
		}
		
		Float *plist;

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
