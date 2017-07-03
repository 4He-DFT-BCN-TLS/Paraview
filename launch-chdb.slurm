#!/bin/bash
#SBATCH --nodes=9
#SBATCH --ntasks=180
#SBATCH --ntasks-per-node=20
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task=1
#SBATCH --time 0:10:0
#SBATCH -J pvb-chdb

module purge
module load intel/16.1.3 intelmpi/5.1.3.210
module load paraview/pvbatch-4.4.0

export I_MPI_WAIT_MODE=1
export OMP_NUM_THREADS=1
export WORKDIR=$(pwd)
cd ${WORKDIR}

mpirun chdb --in-dir=../density-chdb.in --in-files ../density-chdb.in/filter.txt --out-dir pvbatch-out.${SLURM_JOBID} --in-type dat --command-line './draw-slurm-chdb.sh %name% %in-dir% %out-dir%' --out-files %name% --on-error=errors.${SLURM_JOBID} --report=report.${SLURM_JOBID}

infoincidentjob
jobinfo ${SLURM_JOB_ID}
