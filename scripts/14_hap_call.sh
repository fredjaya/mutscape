#!/bin/bash

#PBS -P RDS-FSC-Scape-RW 
#PBS -l select=1:ncpus=4:mem=12GB
#PBS -l walltime=15:00:00

DIR=/scratch/Scape/fred/rtc_idx
REF="GCF_003254395.2_Amel_HAv3.1_genomic"
RECAL_BAM="Larv09_recal_reads.bam"

cd $DIR
module load gatk/3.8.1

gatk -T HaplotypeCaller \
	-R ${REF}.fna \
	-I ${RECAL_BAM} \
	--genotyping_mode DISCOVERY \
	--output_mode EMIT_ALL_SITES \
	-stand_call_conf 10 \
	-o Larv09_raw_variants_recal.vcf \
	-nt 1 -nct 39 # Heldenbrand et al. (2019)	 
