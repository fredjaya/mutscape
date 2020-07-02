#!/bin/bash

#PBS -P RDS-FSC-Scape-RW 
#PBS -l select=1:ncpus=1:mem=1GB
#PBS -l walltime=00:01:00

SCRATCH_PATH="/scratch/Scape/fred"

module load picard/2.18.23

cd ${SCRATCH_PATH}/8_rtc_idx

picard CreateSequenceDictionary \
	R=GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz \
	O=GCF_003254395.2_Amel_HAv3.1_genomic.dict
