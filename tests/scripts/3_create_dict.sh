#!/bin/bash

#PBS -P RDS-FSC-Scape-RW 
#PBS -l select=1:ncpus=1:mem=1GB
#PBS -l walltime=00:01:00

DIR='/scratch/Scape/fred/rtc_idx'

module load picard/2.18.23

cd $DIR

picard CreateSequenceDictionary \
	R=GCF_003254395.2_Amel_HAv3.1_genomic.fna \
	O=GCF_003254395.2_Amel_HAv3.1_genomic.dict
