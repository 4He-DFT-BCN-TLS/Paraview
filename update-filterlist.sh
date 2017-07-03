#!/bin/bash

START=${1}
FINISH=${2}

rm -v ../density-chdb.in/filter.txt
for ID in  $(seq -w ${START} 1 ${FINISH})
do
	echo "density.${ID}.dat" >> ../density-chdb.in/filter.txt
done
echo "New filter list created"
N=$(cat ../density-chdb.in/filter.txt | wc -l)
echo "${N} entries in filter list"
if ((N<180))
then
	if ((N%20>0))
	then
		sed -i "/^#SBATCH --nodes=/c\#SBATCH --nodes=$((N/20+1))" launch-chdb.slurm
	else
		sed -i "/^#SBATCH --nodes=/c\#SBATCH --nodes=$((N/20))" launch-chdb.slurm
	fi
	sed -i "/^#SBATCH --ntasks=/c\#SBATCH --ntasks=${N}" launch-chdb.slurm
else
	sed -i "/^#SBATCH --nodes=/c\#SBATCH --nodes=9" launch-chdb.slurm
	sed -i "/^#SBATCH --ntasks=/c\#SBATCH --ntasks=180" launch-chdb.slurm
fi
echo "'launch-chdb.slurm' updated."