#!/bin/bash

#PBS -P RDS-FSC-Scape-RW
#PBS -l select=1:ncpus=1:mem=1GB
#PBS -l walltime=12:00:00

WORKDIR=/scratch/Scape/fred/nf_tests
NFDIR=/home/fjay0039/mutscape/tests

cd ${WORKDIR}
module load nextflow/20.04.1

nextflow run ${NFDIR}/tests.nf \
	--mode hapCallSingle
	
