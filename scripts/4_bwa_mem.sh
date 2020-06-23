#!/bin/bash

#PBS -P RDS-FSC-Scape-RW 
#PBS -l select=1:ncpus=1:mem=16GB
#PBS -l walltime=6:00:00

# could increase cores, mem and time

SCRATCH_PATH="/scratch/Scape/fred"
READ_PATH="/project/Scape/Trimmomatic/Paired"

module load bwa/0.7.17

mkdir -p ${SCRATCH_PATH}/4_bwa_mem
cd ${SCRATCH_PATH}/4_bwa_mem

bwa mem \
	${SCRATCH_PATH}/1_ref_idx/GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz \
	${READ_PATH}/Larv09_R1_paired.fastq.gz \
	${READ_PATH}/Larv09_R2_paired.fastq.gz \
	> Larv09_pe.sam
