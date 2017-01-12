#include "1deuler.hpp"

int main(int argc, char* argv[])
{
	if(argc < 2)
	{
		std::cout << "Euler1d needs a configuration file name. Quitting.\n";
		return -1;
	}

	std::string confile(argv[1]);
	std::ifstream conf(confile);

	int leftbc, rightbc, temporal_order, N, areatype, maxiter;
	std::vector<double> leftbv(NVARS,0), rightbv(NVARS,0), areas;
	std::string inv_flux, areafile, outputfile, simtype, slope_scheme, rec_scheme, limiter, rkfile, dum;
	double cfl, f_time, L, tol;

	conf >> dum; conf >> simtype;
	if(simtype == "unsteady")
	{
		conf >> dum; conf >> outputfile;
		conf >> dum; conf >> N;
		conf >> dum; conf >> L;
		conf >> dum; conf >> leftbc;
		conf >> dum; conf >> rightbc;
		conf >> dum;
		for(int i = 0; i < NVARS; i++)
			conf >> leftbv[i];
		conf >> dum;
		for(int i = 0; i < NVARS; i++)
			conf >> rightbv[i];
		conf >> dum; conf >> cfl;
		conf >> dum; conf >> inv_flux;
		conf >> dum; conf >> slope_scheme;
		conf >> dum; conf >> rec_scheme;
		conf >> dum; conf >> limiter;
		conf >> dum; conf >> f_time;
		conf >> dum; conf >> temporal_order;
		conf >> dum; conf >> rkfile;
		conf >> dum; conf >> areatype;
		if(areatype == 1)
		{
			areas.resize(N);
			conf >> dum; conf >> areafile;
			std::ifstream areaf(areafile);
			for(int i = 0; i < N; i++)
			{
				areaf >> areas[i];
			}
			areaf.close();
		}
		else
		{
			areas.resize(1);
			areas[0] = 1.0;
		}
	
		std::vector<double> plist;

		std::cout << "Unsteady: " << N << " " << L << " " << leftbc << " " << rightbc << " " << inv_flux << " " << cfl << " " << f_time << " " << temporal_order << std::endl;

		Euler1dExplicit prob(N, L, leftbc, rightbc, leftbv, rightbv, cfl, inv_flux, slope_scheme, rec_scheme, limiter, f_time, temporal_order, rkfile);
		prob.generate_mesh(0,plist);
		prob.set_area(0,areas);

		prob.run();

		prob.postprocess(outputfile);
	}
	else
	{
		conf >> dum; conf >> outputfile;
		conf >> dum; conf >> N;
		conf >> dum; conf >> L;
		conf >> dum; conf >> leftbc;
		conf >> dum; conf >> rightbc;
		conf >> dum;
		for(int i = 0; i < NVARS; i++)
			conf >> leftbv[i];
		conf >> dum;
		for(int i = 0; i < NVARS; i++)
			conf >> rightbv[i];
		conf >> dum; conf >> cfl;
		conf >> dum; conf >> inv_flux;
		conf >> dum; conf >> slope_scheme;
		conf >> dum; conf >> rec_scheme;
		conf >> dum; conf >> limiter;
		conf >> dum; conf >> tol;
		conf >> dum; conf >> maxiter;
		conf >> dum; conf >> areatype;
		if(areatype != 0)
		{
			areas.resize(N);
			conf >> dum; conf >> areafile;
			std::ifstream areaf(areafile);
			for(int i = 0; i < N; i++)
			{
				areaf >> areas[i];
			}
			areaf.close();
		}
		
		std::vector<double> plist;

		Euler1dSteadyExplicit prob(N, L, leftbc, rightbc, leftbv, rightbv, cfl, inv_flux, slope_scheme, rec_scheme, limiter, tol, maxiter);
		prob.generate_mesh(0,plist);
		prob.set_area(areatype,areas);
		std::cout << N << " " << L << " " << leftbc << " " << rightbc << " " << inv_flux << " " << cfl << " " << tol << " " << maxiter << std::endl;

		prob.run();

		prob.postprocess(outputfile);
	}
	conf.close();

	return 0;
}
