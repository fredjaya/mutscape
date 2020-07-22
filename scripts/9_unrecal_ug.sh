#!/bin/bash

#PBS -P RDS-FSC-Scape-RW 
#PBS -l select=1:ncpus=4:mem=20GB
#PBS -l walltime=01:30:00

DIR=/scratch/Scape/fred/rtc_idx
REF="GCF_003254395.2_Amel_HAv3.1_genomic"
BAM="Larv09_realigned_reads.bam"

cd $DIR
module load gatk/3.8.1

gatk -T UnifiedGenotyper \
	-R ${REF}.fna \
	-I ${BAM} \
	-ploidy 1 \
	--genotyping_mode DISCOVERY \
	--output_mode EMIT_ALL_CONFIDENT_SITES \
	-stand_call_conf 30 \
	-o Larv09_ug_variants_conf_sites.vcf \
	-nt 1 -nct 39 # Heldenbrand et al. (2019)	 
