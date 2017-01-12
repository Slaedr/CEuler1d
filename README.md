Quasi-1D Euler
==============

Quasi 1D Euler solver for steady and unsteady problems using 2nd order MUSCL scheme with various limiters and two Riemann fluxes - local Lax-Freidrichs and Van Leer flux vector splitting. Explicit time stepping.

Two variants of MUSCL have been implemented - MUSCLReconstructionG based on least-squares slopes, and MUSCLReconstruction based on Blazek's method. The latter works better at the moment. With the former, note that only Van Albada limiter works.

Accelerator-enabled
-------------------
The code is being designed to take advantage of accelerator devices via OpenACC. For this to work, array storage was converted from `std::vector` to C-style arrays using `malloc()` and `free()`.
Notes:
- Currently, C-style arrays are being used. Try the Array1d and Array2d classes which implement single-pointer flattened storage. `arr[i][j]` type access is enabled by operator overloading. This would ensure that arrays remain contiguous in GPU memory, which I think does not happen now. It would also use less storage.
