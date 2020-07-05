When running RealignTargetCreator, an error with no bam indexing persists

Aim:
- [ ] Re-run pipe from the beginning
- [ ] Decompress, and re-compress REF with bgzip
- [ ] Save all outputs in the same folder
- [ ] Add @RG during bwa-mem, so .bai output from MarkDuplicates is the direct input to RTC


```
# Decompress reference to avoid bgzip indexing issues
cd /scratch/Scape/fred/rtc_idx
bgzip -d GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz

# Generate BWA index
qsub scripts/1_bwa_index.sh

# Generate fasta index 
qsub scripts/2_ref_index.sh

# Create dictionary
qsub scripts/3_create_dict.sh

# Map reads
qsub scripts/4_bwa_mem.sh

# Sort and convert to bam
qsub scripts/5_sortsam.sh

# Mark duplicates
qsub scripts/6_markdups.sh

# Realign target creator
qsub scripts/7_rtc.sh
```
