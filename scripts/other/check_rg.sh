#!/bin/bash

# Check if read groups vary across samples

module load samtools/1.9

for i in /project/Scape/Picard/*sorted_add_rg.bam
do
    samtools view -H $i | grep '@RG'
done
