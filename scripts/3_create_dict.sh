#!/bin/bash

#PBS -P RDS-FSC-Scape-RW 
#PBS -l select=1:ncpus=1:mem=1GB
#PBS -l walltime=00:05:00

SCRATCH_PATH="/scratch/Scape/fred"

module load picard/2.18.23

mkdir -p ${SCRATCH_PATH}/3_create_dict
cd ${SCRATCH_PATH}/3_create_dict

picard CreateSequenceDictionary \
	R=${SCRATCH_PATH}/GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz \
	O=${SCRATCH_PATH}/3_create_dict/GCF_003254395.2_Amel_HAv3.1_genomic.dict
