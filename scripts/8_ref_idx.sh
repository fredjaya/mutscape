#!/bin/bash

#PBS -P RDS-FSC-Scape-RW 
#PBS -l select=1:ncpus=1:mem=1GB
#PBS -l walltime=00:5:00

SCRATCH_PATH='/scratch/Scape/fred'

module load bwa/0.7.17

mkdir -p ${SCRATCH_PATH}/8_rtc_idx
cd ${SCRATCH_PATH}/8_rtc_idx

bwa index GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz
