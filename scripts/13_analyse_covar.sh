#!/bin/bash

#PBS -P RDS-FSC-Scape-RW 
#PBS -l select=1:ncpus=1:mem=1GB
#PBS -l walltime=00:01:00

DIR=/scratch/Scape/fred/rtc_idx
REF="GCF_003254395.2_Amel_HAv3.1_genomic"
PRE="Larv09_pre_recal.table"
POST="Larv09_post_recal.table"

cd $DIR
module load gatk/3.8.1
module load R/3.6.0

gatk -T AnalyzeCovariates \
	-R ${REF}.fna \
	-before ${PRE} \
	-after ${POST} \
	-plots Larv09_bqsr_report.pdf

# 11_post_bsqr.sh may not be neccessary, and can be incorporated 
