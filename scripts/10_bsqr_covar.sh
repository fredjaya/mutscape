#!/bin/bash

#PBS -P RDS-FSC-Scape-RW 
#PBS -l select=1:ncpus=2:mem=12GB
#PBS -l walltime=4:00:00

DIR=/scratch/Scape/fred/rtc_idx
REF="GCF_003254395.2_Amel_HAv3.1_genomic"
BAM="Larv09_realigned_reads.bam"
VCF="Larv09_raw_variants_conf_sites.vcf"
RECAL="Larv09_recal_covar.table"

cd $DIR
module load gatk/3.8.1

gatk -T BaseRecalibrator \
	-R ${REF}.fna \
	-I ${BAM} \
	-knownSites ${VCF} \
	-o Larv09_pre_recal.table
# add -nt -nct 
