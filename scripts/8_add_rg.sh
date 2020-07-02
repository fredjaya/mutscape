#!/bin/bash

#PBS -P RDS-FSC-Scape-RW 
#PBS -l select=1:ncpus=1:mem=16GB
#PBS -l walltime=00:30:00

SCRATCH_PATH="/scratch/Scape/fred"

module load picard/2.7.1

cd ${SCRATCH_PATH}/8_rtc_idx

picard AddOrReplaceReadGroups \
	I=Larv09_pe_marked_dups.bam \
	O=Larv09_marked_dups_rg.bam \
	RGID=Larv09 \
	RGLB=lib1 \
	RGPL=ILLUMINA \
	RGPU=unit1 \
	RGSM=Larv09
