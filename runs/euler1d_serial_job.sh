#! /bin/bash
#PBS -l nodes=1:ppn=1:sandybridge
#PBS -l walltime=00:05:00
#PBS -A rck-371-aa
#PBS -o nozzle_serial.log
#PBS -e nozzle_serial.err
#PBS -N nozzle_serial

# to be executed from Euler1d_root_dir/runs/ or equivalent

# ppn=16:sandybridge will always be run on a SW2 sandy bridge node, according to http://www.hpc.mcgill.ca/index.php/starthere/81-doc-pages/91-guillimin-job-submit

cd $PBS_O_WORKDIR
cat /proc/cpuinfo | grep 'model name'
../build/euler1d-acc ../cases/nozzle.control
