#!/bin/bash

#PBS -P RDS-FSC-Scape-RW
#PBS -l select=1:ncpus=1:mem=1GB
#PBS -l walltime=48:00:00

WORKDIR=/scratch/Scape/fred/2008_manual
NFDIR=/home/fjay0039/mutscape

cd ${WORKDIR}
module load nextflow/20.04.1

nextflow run ${NFDIR}/main.nf \
	--mode haplotypeCallerUnrecal \
	-profile pbs \
	-resume
