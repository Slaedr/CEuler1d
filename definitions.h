#ifndef __DEFITIONS_H
#define __DEFINITIONS_H

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
#define gamma 1.4

/// MUSCL factor
const Float k = 1.0/3.0;

/// Specific gas constant of air
//const Float R = 287.1;
#define R 1716.0

#define Cv R/(g-1.0)

/// Constant for MUSCL reconstruction
const Float muscl_k = 1.0/3.0;

void matprint(Float const *const *const mat);

// OpenACC constants
#define NVIDIA_VECTOR_LENGTH 32

#endif
