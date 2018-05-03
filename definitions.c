#include "definitions.h"

void matprint(Float const *const *const mat, size_t m, size_t n)
{
	printf("\n");
	for(int i = 0; i < m; i++)
	{
		for(int j = 0; j < n; j++)
			printf("%f\t", mat[i][j]);
		printf("\n");
	}
	printf("\n");
}

