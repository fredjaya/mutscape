#!/bin/bash

#PBS -P RDS-FSC-Scape-RW 
#PBS -l select=1:ncpus=2:mem=64GB
#PBS -l walltime=6:00:00

SCRATCH_PATH="/scratch/Scape/fred/"

module load picard
cd ${SCRATCH_PATH}

picard MarkDuplicates \
	I=1_bwa/Larv01_pe.sam \
	O=3_markdups/Larv01_pe_marked_duplicates.bam \
	M=3_markdups/marked_dup_metrics.txt \
	ASO=queryname
