#!/bin/bash

#PBS -P RDS-FSC-Scape-RW 
#PBS -l select=1:ncpus=2:mem=32GB
#PBS -l walltime=2:00:00

SCRATCH_PATH="/scratch/Scape/fred"

module load picard/2.7.1

mkdir -p ${SCRATCH_PATH}/6_markdups
cd ${SCRATCH_PATH}/6_markdups

picard MarkDuplicates \
	I=${SCRATCH_PATH}/5_sorted_sam/Larv09_pe_sorted.bam \
	O=Larv09_pe_marked_dups.bam \
	M=marked_dup_metrics.txt \
