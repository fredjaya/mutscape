#!/bin/bash

#PBS -P RDS-FSC-Scape-RW 
#PBS -l select=1:ncpus=2:mem=16GB
#PBS -l walltime=05:00:00

DIR=/scratch/Scape/fred/rtc_idx
REF="GCF_003254395.2_Amel_HAv3.1_genomic"
BAM="Larv09_realigned_reads.bam"

cd $DIR
module load gatk/3.8.1

gatk -T HaplotypeCaller \
	-R ${REF}.fna \
	-I ${BAM} \
	--genotyping_mode DISCOVERY \
	--output_mode EMIT_ALL_SITES \
	-stand_call_conf 10 \
	-o Larv09_raw_variants.vcf \
	-nt 1 -nct 39 # Heldenbrand et al. (2019)	 
