#ifndef CEULER1D_DEFINITIONS_H
#define CEULER1D_DEFINITIONS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define ZERO_TOL 2.2e-16
#define SMALL_NUMBER 1e-10

#define NVARS 3

typedef double Float;

const Float PI = 3.14159265359;

/// adiabatic index
//#define gamma 1.4
const Float gamma = 1.4;

// MUSCL factor (not needed)
//const Float k = 1.0/3.0;

/// Specific gas constant of air
//const Float R = 287.1; // metric units
#define GAS_CONSTANT 1716.0

/// Write out a dense matrix to console
void matprint(Float const *const *const mat, size_t m, size_t n);

/// OpenACC constants
#define NVIDIA_VECTOR_LENGTH 32

/// Hard assertion
#define TASSERT(x) if(!(x)) printf("Assertion failed!\n"); exit(-1)

#endif
