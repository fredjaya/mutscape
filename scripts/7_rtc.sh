#!/bin/bash

#PBS -P RDS-FSC-Scape-RW 
#PBS -l select=1:ncpus=1:mem=8GB
#PBS -l walltime=00:30:00

DIR=/scratch/Scape/fred/rtc_idx
REF="GCF_003254395.2_Amel_HAv3.1_genomic"
BAM="Larv09_marked_dups.bam"

cd $DIR
module load gatk/3.8.1

gatk -T RealignerTargetCreator \
	-R ${REF}.fna \
	-I ${BAM} \
	-o Larv09_target_intervals.list
