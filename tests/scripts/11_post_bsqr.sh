#!/bin/bash

#PBS -P RDS-FSC-Scape-RW 
#PBS -l select=1:ncpus=2:mem=10GB
#PBS -l walltime=4:00:00

DIR=/scratch/Scape/fred/rtc_idx
REF="GCF_003254395.2_Amel_HAv3.1_genomic"
BAM="Larv09_realigned_reads.bam"
RECAL="Larv09_pre_recal.table"
VCF="Larv09_raw_variants.vcf "

cd $DIR
module load gatk/3.8.1

gatk -T BaseRecalibrator \
	-R ${REF}.fna \
	-I ${BAM} \
	-knownSites ${VCF} \
	-BQSR ${RECAL} \
	-o Larv09_post_recal.table \
	-nct 24 # Heldenbrand (2019) == 40; HPC limit == 24
