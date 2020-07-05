#!/bin/bash

#PBS -P RDS-FSC-Scape-RW 
#PBS -l select=1:ncpus=2:mem=16GB
#PBS -l walltime=05:00:00

DIR=/scratch/Scape/fred/rtc_idx
REF="GCF_003254395.2_Amel_HAv3.1_genomic"
BAM="Larv09_marked_dups.bam"

cd $DIR
module load gatk/3.8.1

gatk -T IndelRealigner \
	-R ${REF}.fna \
	-I ${BAM} \
	-targetIntervals Larv09_target_intervals.list \
	-o Larv09_realigned_reads.bam
