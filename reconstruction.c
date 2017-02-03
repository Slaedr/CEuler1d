/** \file reconstruction.c
 * \brief Implementation of reconstruction schemes.
 * \author Aditya Kashi
 * \date Jan 2017
 */

#include "reconstruction.h"

void compute_noReconstruction(const Grid* const grid, Euler1d *const sim)
{
#pragma acc parallel loop gang worker vector present(grid, sim) device_type(nvidia) vector_length(NVIDIA_VECTOR_LENGTH)
	for(size_t i = 0; i <= grid->N; i++)
	{
#pragma acc loop seq
		for(size_t j = 0; j < NVARS; j++)
		{
			sim->prleft[i][j] = sim->prim[i][j];
			sim->prright[i][j] = sim->prim[i+1][j];
		}
	}
}

void compute_MUSCLReconstruction(const Grid *const grid, Euler1d *const sim)
{
	// interior faces
	#pragma acc parallel present(grid, sim) device_type(nvidia) vector_length(NVIDIA_VECTOR_LENGTH) //num_workers(NVIDIA_WORKERS_PER_GANG)
	{
		#pragma acc loop gang worker vector 
		for(size_t i = 1; i <= grid->N-1; i++)
		{
			Float denL, denR, num, rL, rR;
			
			#pragma acc loop seq
			for(size_t j = 0; j < NVARS; j++)
			{
				denL = sim->prim[i][j]-sim->prim[i-1][j];
				denR = sim->prim[i+2][j]-sim->prim[i+1][j];
				num = sim->prim[i+1][j]-sim->prim[i][j];

				if(fabs(denL) > ZERO_TOL*10)
				{
					rL = num/denL;
					Float ss = vanalbada_limiter_function(rL);
					sim->prleft[i][j] = sim->prim[i][j] + 0.25*vanalbada_limiter_function(rL) * ((1.0+(sim->muscl_k)*ss)*num + (1.0-(sim->muscl_k*ss))*denL);
				}
				else
					sim->prleft[i][j] = sim->prim[i][j];

				if(fabs(denR) > ZERO_TOL*10)
				{
					rR = num/denR;
					Float ss = vanalbada_limiter_function(rR);
					sim->prright[i][j] = sim->prim[i+1][j] - 0.25*vanalbada_limiter_function(rR) * ((1.0+(sim->muscl_k)*ss)*num + (1.0-(sim->muscl_k*ss))*denR);
				}
				else
					sim->prright[i][j] = sim->prim[i+1][j];
			}
		}
	}
	
	// boundaries
	#pragma acc parallel present(grid, sim) num_gangs(1)
	{
		Float denL, denR, num, rL, rR;
		Float slope, cc, um0, xm0, umN, xmN;

		#pragma acc loop seq
		for(size_t j = 0; j < NVARS; j++)
		{
			// extrapolate variables

			slope = (sim->prim[1][j] - sim->prim[0][j])/(grid->x[1]-grid->x[0]);
			cc = sim->prim[0][j] - (sim->prim[1][j]-sim->prim[0][j])/(grid->x[1]-grid->x[0])*grid->x[0];
			xm0 = grid->x[0] - (grid->x[1]-grid->x[0]);
			um0 = slope*xm0 + cc;

			slope = (sim->prim[grid->N+1][j] - sim->prim[grid->N][j])/(grid->x[grid->N+1]-grid->x[(grid->N)]);
			cc = sim->prim[grid->N][j] - (sim->prim[grid->N+1][j]-sim->prim[grid->N][j])/(grid->x[grid->N+1]-grid->x[grid->N])*grid->x[grid->N];
			xmN = grid->x[grid->N+1] + grid->x[grid->N]-grid->x[grid->N-1];
			umN = slope*xmN + cc;

			// left
			denL = sim->prim[0][j] - um0;
			denR = sim->prim[2][j]-sim->prim[1][j];
			num = sim->prim[1][j]-sim->prim[0][j];

			if(fabs(denL) > ZERO_TOL*10)
			{
				rL = num/denL;
				Float ss = vanalbada_limiter_function(rL);
				sim->prleft[0][j] = sim->prim[0][j] + 0.25*vanalbada_limiter_function(rL) * ((1.0+(sim->muscl_k)*ss)*num + (1.0-(sim->muscl_k*ss))*denL);
			}
			else
				sim->prleft[0][j] = sim->prim[0][j];

			if(fabs(denR) > ZERO_TOL*10)
			{
				rR = num/denR;
				Float ss = vanalbada_limiter_function(rR);
				sim->prright[0][j] = sim->prim[1][j] - 0.25*vanalbada_limiter_function(rR) * ((1.0+(sim->muscl_k*ss))*num + (1.0-(sim->muscl_k*ss))*denR);
			}
			else
				sim->prright[0][j] = sim->prim[1][j];
			
			// right
			denL = sim->prim[grid->N][j]-sim->prim[grid->N-1][j];
			denR = umN - sim->prim[grid->N+1][j];
			num = sim->prim[grid->N+1][j]-sim->prim[grid->N][j];

			if(fabs(denL) > ZERO_TOL*10)
			{
				rL = num/denL;
				Float ss = vanalbada_limiter_function(rL);
				sim->prleft[grid->N][j] = sim->prim[grid->N][j] + 0.25*vanalbada_limiter_function(rL) * ((1.0+sim->muscl_k*ss)*num + (1.0-sim->muscl_k*ss)*denL);
			}
			else
				sim->prleft[grid->N][j] = sim->prim[grid->N][j];

			if(fabs(denR) > ZERO_TOL*10)
			{
				rR = num/denR;
				Float ss = vanalbada_limiter_function(rR);
				sim->prright[grid->N][j] = sim->prim[grid->N+1][j] - 0.25*vanalbada_limiter_function(rR) * ((1.0+sim->muscl_k*ss)*num + (1.0-sim->muscl_k*ss)*denR);
			}
			else
				sim->prright[grid->N][j] = sim->prim[grid->N+1][j];
		}
	}
}

void compute_MUSCLReconstruction_unlimited(const Grid *const grid, Euler1d *const sim)
{
	// interior faces
	#pragma acc parallel present(grid, sim) device_type(nvidia) vector_length(NVIDIA_VECTOR_LENGTH) //num_workers(NVIDIA_WORKERS_PER_GANG)
	{
		#pragma acc loop gang worker vector 
		for(size_t i = 1; i <= grid->N-1; i++)
		{
			Float denL, denR, num, rL, rR;
			
			#pragma acc loop seq
			for(size_t j = 0; j < NVARS; j++)
			{
				denL = sim->prim[i][j]-sim->prim[i-1][j];
				denR = sim->prim[i+2][j]-sim->prim[i+1][j];
				num = sim->prim[i+1][j]-sim->prim[i][j];

				if(fabs(denL) > ZERO_TOL*10)
				{
					rL = num/denL;
					sim->prleft[i][j] = sim->prim[i][j] + 0.25 * ((1.0+(sim->muscl_k))*num + (1.0-(sim->muscl_k))*denL);
				}
				else
					sim->prleft[i][j] = sim->prim[i][j];

				if(fabs(denR) > ZERO_TOL*10)
				{
					rR = num/denR;
					sim->prright[i][j] = sim->prim[i+1][j] - 0.25* ((1.0+(sim->muscl_k))*num + (1.0-(sim->muscl_k))*denR);
				}
				else
					sim->prright[i][j] = sim->prim[i+1][j];
			}
		}
	}
	
	// boundaries
	#pragma acc parallel present(grid, sim) num_gangs(1)
	{
		Float denL, denR, num, rL, rR;
		Float slope, cc, um0, xm0, umN, xmN;

		#pragma acc loop seq
		for(size_t j = 0; j < NVARS; j++)
		{
			// extrapolate variables

			slope = (sim->prim[1][j] - sim->prim[0][j])/(grid->x[1]-grid->x[0]);
			cc = sim->prim[0][j] - (sim->prim[1][j]-sim->prim[0][j])/(grid->x[1]-grid->x[0])*grid->x[0];
			xm0 = grid->x[0] - (grid->x[1]-grid->x[0]);
			um0 = slope*xm0 + cc;

			slope = (sim->prim[grid->N+1][j] - sim->prim[grid->N][j])/(grid->x[grid->N+1]-grid->x[(grid->N)]);
			cc = sim->prim[grid->N][j] - (sim->prim[grid->N+1][j]-sim->prim[grid->N][j])/(grid->x[grid->N+1]-grid->x[grid->N])*grid->x[grid->N];
			xmN = grid->x[grid->N+1] + grid->x[grid->N]-grid->x[grid->N-1];
			umN = slope*xmN + cc;

			// left
			denL = sim->prim[0][j] - um0;
			denR = sim->prim[2][j]-sim->prim[1][j];
			num = sim->prim[1][j]-sim->prim[0][j];

			if(fabs(denL) > ZERO_TOL*10)
			{
				rL = num/denL;
				sim->prleft[0][j] = sim->prim[0][j] + 0.25 * ((1.0+(sim->muscl_k))*num + (1.0-(sim->muscl_k))*denL);
			}
			else
				sim->prleft[0][j] = sim->prim[0][j];

			if(fabs(denR) > ZERO_TOL*10)
			{
				rR = num/denR;
				sim->prright[0][j] = sim->prim[1][j] - 0.25 * ((1.0+(sim->muscl_k))*num + (1.0-(sim->muscl_k))*denR);
			}
			else
				sim->prright[0][j] = sim->prim[1][j];
			
			// right
			denL = sim->prim[grid->N][j]-sim->prim[grid->N-1][j];
			denR = umN - sim->prim[grid->N+1][j];
			num = sim->prim[grid->N+1][j]-sim->prim[grid->N][j];

			if(fabs(denL) > ZERO_TOL*10)
			{
				rL = num/denL;
				sim->prleft[grid->N][j] = sim->prim[grid->N][j] + 0.25 * ((1.0+sim->muscl_k)*num + (1.0-sim->muscl_k)*denL);
			}
			else
				sim->prleft[grid->N][j] = sim->prim[grid->N][j];

			if(fabs(denR) > ZERO_TOL*10)
			{
				rR = num/denR;
				sim->prright[grid->N][j] = sim->prim[grid->N+1][j] - 0.25 * ((1.0+sim->muscl_k)*num + (1.0-sim->muscl_k)*denR);
			}
			else
				sim->prright[grid->N][j] = sim->prim[grid->N+1][j];
		}
	}
}
