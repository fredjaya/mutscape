#!/bin/bash

#PBS -P RDS-FSC-Scape-RW 
#PBS -l select=1:ncpus=2:mem=32GB
#PBS -l walltime=2:00:00

SCRATCH_PATH="/scratch/Scape/fred"
REF="GCF_003254395.2_Amel_HAv3.1_genomic"
BAM="Larv09_marked_dups_rg.bam"

module load gatk/3.8.1

cd ${SCRATCH_PATH}/8_rtc_idx

gatk -T RealignerTargetCreator \
	-R ${REF}.fna \
	-I ${BAM} \
	-o Larv09_target_intervals.list
