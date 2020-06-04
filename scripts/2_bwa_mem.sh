#!/bin/bash

#PBS -P RDS-FSC-Scape-RW 
#PBS -l select=1:ncpus=1:mem=16GB
#PBS -l walltime=12:00:00

SCRATCH_PATH="/scratch/Scape/fred/"
READ_PATH="/project/Scape/Trimmomatic/Paired/"

module load bwa
cd ${SCRATCH_PATH}1_bwa/

bwa mem \
	GCA_003254395.2_Amel_HAv3.1_genomic.fna.gz \
	${READ_PATH}Larv01_R1_paired.fastq.gz \
	${READ_PATH}Larv01_R2_paired.fastq.gz \
	> Larv01_pe.sam
