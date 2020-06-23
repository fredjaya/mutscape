#!/bin/bash

#PBS -P RDS-FSC-Scape-RW 
#PBS -l select=1:ncpus=2:mem=16GB
#PBS -l walltime=1:00:00

SCRATCH_PATH="/scratch/Scape/fred"

module load samtools/1.9

mkdir -p ${SCRATCH_PATH}/5_sorted_sam
cd ${SCRATCH_PATH}/5_sorted_sam

samtools sort -o Larv09_pe_sorted.bam -O bam ${SCRATCH_PATH}/4_bwa_mem/Larv09_pe.sam
