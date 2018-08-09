#!/bin/bash

### Script assumes that all the wave-function files are located in '../he-wfs/' relative to this script.
### CALL:	./update-filterlist.sh <START ID w/ leading zeros> <STOP ID w/ leading zeros>
### EXMPL:	./update-filterlist.sh 0005 0234

START=${1}
FINISH=${2}

rm -v ../he-wfs/filter.txt
for ID in  $(seq -w ${START} 1 ${FINISH})
do
	echo "density.${ID}.dat" >> ../he-wfs/filter.txt
done
echo "New filter list created"
N=$(cat ../he-wfs/filter.txt | wc -l)
M=$((N+1)) # We need one extra MPI task for the master program that controls the MPI slaves
echo "${N} entries in filter list + $((M-N)) task for the MPI master task"
if ((M<180))
then
	if ((M%20>0))
	then
		sed -i "/^#SBATCH --nodes=/c\#SBATCH --nodes=$((M/20+1))" launch-chdb.slurm
	else
		sed -i "/^#SBATCH --nodes=/c\#SBATCH --nodes=$((M/20))" launch-chdb.slurm
	fi
	sed -i "/^#SBATCH --ntasks=/c\#SBATCH --ntasks=${M}" launch-chdb.slurm
else
	sed -i "/^#SBATCH --nodes=/c\#SBATCH --nodes=9" launch-chdb.slurm
	sed -i "/^#SBATCH --ntasks=/c\#SBATCH --ntasks=180" launch-chdb.slurm
fi
echo "'launch-chdb.slurm' updated."

