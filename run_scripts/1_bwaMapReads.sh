#!/bin/bash

#PBS -P RDS-FSC-Scape-RW
#PBS -l select=1:ncpus=1:mem=1GB
#PBS -l walltime=120:00:00

WORKDIR=/scratch/Scape/fred/2008_manual
NFDIR=/home/fjay0039/mutscape

cd ${WORKDIR}
module load nextflow/20.04.1

nextflow run ${NFDIR}/main.nf \
	--mode bwaMapReads \
	-profile pbs \
	--ref "/project/Scape/Trimmomatic/Paired/Fdrone_{R1,R2}_paired.fastq.gz"
