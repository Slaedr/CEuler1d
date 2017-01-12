#! /bin/bash

# set up environment for compiling OpenACC program

module load CUDA/7.5.18 
module load PGI/16.3-GCC-4.9.3-2.25

export CXX=pgc++
unset DEBUG
export BUILD_WITH_ACC=1
export BUILD_FOR_MULTICORE=1

echo "CUDA loaded, PGI loaded, CXX set to pgc++, BUILD_WITH_ACC set to 1, BUILD_FOR_MULTICORE set to 1 and DEBUG set to 0"
