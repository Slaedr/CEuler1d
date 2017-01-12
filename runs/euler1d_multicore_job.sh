#! /bin/bash
#PBS -l nodes=1:ppn=16:sandybridge
#PBS -l walltime=00:20:00
#PBS -A rck-371-aa
#PBS -o sodst-100k_acc_cpus-16threads-cellflux.log
#PBS -e sodst_acc_cpus.err
#PBS -N sodst_acc_cpus

# to be executed from Euler1d-acc/runs/ or equivalent

# ppn=16:sandybridge will always be run on a SW2 sandy bridge node, according to http://www.hpc.mcgill.ca/index.php/starthere/81-doc-pages/91-guillimin-job-submit

module load CUDA/7.5.18
module load PGI/16.3-GCC-4.9.3-2.25
#export PGI_ACC_TIME=1
export ACC_NUM_CORES=16

cd $PBS_O_WORKDIR
cat /proc/cpuinfo | grep 'model name'
../build/euler1d-acc ../cases/sodst.control
