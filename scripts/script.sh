#!/bin/bash
#PBS -q pdlab
#PBS -N mainProgram16n2
#PBS -j oe
#PBS -l nodes=2:ppn=8
#PBS -l walltime=01:00:00

module load mpich-x86_64

cd $PBS_O_WORKDIR


NHOSTS = 'cat $PBS_NODEFILE|wc -l'
echo "nodes -> $NHOSTS" > log.log

mpdboot -n $NHOSTS -f $PBS_NODEFILE -r ssh

mpiexec -np 16 ../../mainProgram 22 2 5 > log.txt 2>&1

mpdallexit

