#!/bin/bash

#PBS -P RDS-FSC-Scape-RW 
#PBS -l select=1:ncpus=2:mem=32GB
#PBS -l walltime=00:45:00

SCRATCH_PATH="/scratch/Scape/fred"

module load picard/2.7.1

cd ${SCRATCH_PATH}/8_rtc_idx

picard MarkDuplicates \
	I=${SCRATCH_PATH}/5_sorted_sam/Larv09_pe_sorted.bam \
	O=Larv09_pe_marked_dups.bam \
	M=marked_dup_metrics.txt \
	CREATE_INDEX=true

