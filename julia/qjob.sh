
#!/bin/sh
#QSUB -queue i18cpu
#QSUB -node 1
#QSUB -mpi 24
#QSUB -omp 1
#QSUB -place pack
#QSUB -over false
#PBS -l walltime=0:30:00
#PBS -N srvo3
cd $PBS_O_WORKDIR

. /etc/profile.d/modules.sh
module remove mpt
module load intel-mpi

date > output
mpijob julia ~/repos/lecture/julia/main.jl 2d.ini >> output
date >> output
