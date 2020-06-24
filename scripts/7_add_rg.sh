#!/bin/bash

#PBS -P RDS-FSC-Scape-RW 
#PBS -l select=1:ncpus=2:mem=16GB
#PBS -l walltime=1:00:00

SCRATCH_PATH="/scratch/Scape/fred"

module load picard/2.7.1

mkdir -p ${SCRATCH_PATH}/7_add_rg
cd ${SCRATCH_PATH}

picard AddOrReplaceReadGroups \
	I=6_markdups/Larv09_pe_marked_dups.bam \
	O=7_add_rg/Larv09_marked_dups_rg.bam \
	RGID=Larv09 \
	RGLB=lib1 \
	RGPL=ILLUMINA \
	RGPU=unit1 \
	RGSM=Larv09
