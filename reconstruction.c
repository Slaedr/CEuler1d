/** \file reconstruction.c
 * \brief Implementation of reconstruction schemes.
 * \author Aditya Kashi
 * \date Jan 2017
 */

#include "reconstruction.h"

void compute_MUSCLReconstruction(const size_t N, Float const *const x, Float const *const *const u, Float const *const *const _dudx, Float * const * const uleft, 
		Float * const * const uright, const Float k)

{
	// interior faces
	#pragma acc parallel present(uleft, uright, u, x, lim, k, N) device_type(nvidia) vector_length(NVIDIA_VECTOR_LENGTH)
	{
		#pragma acc loop gang worker vector 
		for(size_t i = 1; i <= N-1; i++)
		{
			Float denL, denR, num, rL, rR;
			
			#pragma acc loop seq
			for(size_t j = 0; j < NVARS; j++)
			{
				denL = u[i][j]-u[i-1][j];
				denR = u[i+2][j]-u[i+1][j];
				num = u[i+1][j]-u[i][j];

				if(fabs(denL) > ZERO_TOL*10)
				{
					rL = num/denL;
					uleft[i][j] = u[i][j] + 0.25*vanalbada_limiter_function(rL) * ((1.0+(*k))*num + (1.0-(*k))*denL);
				}
				else
					uleft[i][j] = u[i][j];

				if(fabs(denR) > ZERO_TOL*10)
				{
					rR = num/denR;
					uright[i][j] = u[i+1][j] - 0.25*vanalbada_limiter_function(rR) * ((1.0+(*k))*num + (1.0-(*k))*denR);
				}
				else
					uright[i][j] = u[i+1][j];
			}
		}
	}
	
	// boundaries
	#pragma acc parallel present(uleft, uright, u, x, lim, k, N) num_gangs(1)
	{
		Float denL, denR, num, rL, rR;
		Float slope, cc, um0, xm0, umN, xmN;

		#pragma acc loop seq
		for(size_t j = 0; j < NVARS; j++)
		{
			// extrapolate variables

			slope = (u[1][j] - u[0][j])/(x[1]-x[0]);
			cc = u[0][j] - (u[1][j]-u[0][j])/(x[1]-x[0])*x[0];
			xm0 = x[0] - (x[1]-x[0]);
			um0 = slope*xm0 + cc;

			slope = (u[N+1][j] - u[*N][j])/(x[N+1]-x[(*N)]);
			cc = u[N][j] - (u[N+1][j]-u[N][j])/(x[N+1]-x[N])*x[N];
			xmN = x[N+1] + x[N]-x[N-1];
			umN = slope*xmN + cc;

			// left
			denL = u[0][j] - um0;
			denR = u[2][j]-u[1][j];
			num = u[1][j]-u[0][j];

			if(fabs(denL) > ZERO_TOL*10)
			{
				rL = num/denL;
				uleft[0][j] = u[0][j] + 0.25*vanalbada_limiter_function(rL) * ((1.0+(*k))*num + (1.0-(*k))*denL);
			}
			else
				uleft[0][j] = u[0][j];

			if(fabs(denR) > ZERO_TOL*10)
			{
				rR = num/denR;
				uright[0][j] = u[1][j] - 0.25*vanalbada_limiter_function(rR) * ((1.0+(*k))*num + (1.0-(*k))*denR);
			}
			else
				uright[0][j] = u[1][j];
			
			// right
			denL = u[N][j]-u[N-1][j];
			denR = umN - u[N+1][j];
			num = u[N+1][j]-u[N][j];

			if(fabs(denL) > ZERO_TOL*10)
			{
				rL = num/denL;
				uleft[N][j] = u[N][j] + 0.25*vanalbada_limiter_function(rL) * ((1.0+k)*num + (1.0-k)*denL);
			}
			else
				uleft[N][j] = u[N][j];

			if(fabs(denR) > ZERO_TOL*10)
			{
				rR = num/denR;
				uright[N][j] = u[N+1][j] - 0.25*vanalbada_limiter_function(rR) * ((1.0+k)*num + (1.0-k)*denR);
			}
			else
				uright[N][j] = u[N+1][j];
		}
	}
}
