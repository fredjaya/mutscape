#!/bin/bash

#PBS -P RDS-FSC-Scape-RW 
#PBS -l select=1:ncpus=1:mem=1GB
#PBS -l walltime=00:30:00

SCRATCH_PATH='/scratch/Scape/fred'

module load bwa/0.7.17

mkdir -p ${SCRATCH_PATH}/1_ref_idx
cd ${SCRATCH_PATH}/1_ref_idx

bwa index ${SCRATCH_PATH}/GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz
