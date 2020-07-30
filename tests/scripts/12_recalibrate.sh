#!/bin/bash

#PBS -P RDS-FSC-Scape-RW 
#PBS -l select=1:ncpus=3:mem=30GB
#PBS -l walltime=3:00:00

# Check how much memory is required
DIR=/scratch/Scape/fred/rtc_idx
REF="GCF_003254395.2_Amel_HAv3.1_genomic"
BAM="Larv09_realigned_reads.bam"
VCF="Larv09_raw_variants_conf_sites.vcf"

cd $DIR
module load gatk/3.8.1

gatk -T PrintReads \
	-R ${REF}.fna \
	-I ${BAM} \
	-BQSR Larv09_pre_recal.table \
	-o Larv09_recal_reads.bam \
	-nct 3 # Heldenbrand (2019) 

