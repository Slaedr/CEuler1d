Quasi-1D Euler
==============

Quasi 1D Euler solver for steady and unsteady problems using 2nd order MUSCL scheme with various limiters and two Riemann fluxes - local Lax-Freidrichs (LLF) and Van Leer flux vector splitting. Explicit time stepping. Currently does not handle steady-state calculations with shocks.

Accelerator-enabled
-------------------
The code is being designed to take advantage of accelerator devices via OpenACC.

Build scripts
-------------
The scripts "build_cpu.sh", "build_multicore.sh" and "build_gpu.sh" setup the compilation environment on the Guillimin cluster. 
Set the DEBUG variable to compile a debug version. 
A directory called "build" needs to be made before running make for the first time.

Running
-------
Make sure to run the executable from within the runs/ directory, as the locations of the cross-sectional area file etc are specified relative to this.

Visualization
-------------
Python, along with Numpy and Matplotlib, is required to run "plot.py" (for steady flow results) and "plot-unsteady.py" (for unsteady cases). These scripts require as arguments (1) the file containing data output by the code, and (2) a name for a title of the plot which is also used as the filename of the plot file that is automatically saved.
