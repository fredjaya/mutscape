#!/bin/bash

#PBS -P RDS-FSC-Scape-RW 
#PBS -l select=1:ncpus=2:mem=32GB
#PBS -l walltime=00:45:00

DIR=/scratch/Scape/fred/rtc_idx

cd $DIR
module load picard/2.7.1

picard MarkDuplicates \
	I=Larv09_pe_sorted.bam \
	O=Larv09_marked_dups.bam \
	M=marked_dup_metrics.txt \
	CREATE_INDEX=true

