#!/bin/bash

#PBS -P RDS-FSC-Scape-RW 
#PBS -l select=1:ncpus=2:mem=32GB
#PBS -l walltime=2:00:00

DIR=/scratch/Scape/fred/rtc_idx
REF="GCF_003254395.2_Amel_HAv3.1_genomic"
BAM="Larv09_marked_dups_rg.bam"

cd $DIR
module load gatk/3.8.1

gatk -T RealignerTargetCreator \
	-R ${REF}.fna \
	-I ${BAM} \
	-o Larv09_target_intervals.list
