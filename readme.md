# mutscape

Taking an iterative approach to pipeline development i.e. running nextflow processes individually to generate results, over developing a comprehensive pipeline with singularity.

## Protocol

### Prior to nextflow

Download reference genome, generate bwa and .fasta indexes
```
# Download the reference genome for _Apis mellifera_
cd /scratch/Scape/fred/
rsync --copy-links --times --verbose \
rsync://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/254/395/GCF_003254395.2_Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz \ 
/scratch/Scape/fred/ # Artemis HPC

# Decompress reference to avoid bgzip indexing issues
cd /scratch/Scape/fred/rtc_idx
bgzip -d GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz

# Generate BWA index
qsub scripts/1_bwa_index.sh

# Generate fasta index 
qsub scripts/2_ref_index.sh

# Create dictionary
qsub scripts/3_create_dict.sh
```

### Preparing sequence data (`.fastq` to `.bam`)

#### 1. Map trimmed reads to reference genome
```
qsub run_scripts/1_bwaMapReads.sh   
```

#### 2. Convert .sam into .sorted .bam
```
qsub run_scripts/2_samtoolsSort.sh
```

#### 3. Mark duplicates
```
qsub run_scripts/3_markDuplicates.sh
```
First run: a few files aborted due to files exceeding walltime

Remaining files run separately using (`sampleId_md_ch`). Nextflow might be mixing up the sample names. Renaming `*_marked_dups.bam` files based on @RG .bam headers (filename -> @RG ):
- Larv01 -> Larv04
- Larv02 -> Worker
- Larv03 -> Larv09
- Larv04 -> Larv08
- Larv05 -> Larv01 
- Larv06 -> Larv03
- Larv07 -> Larv07 *
- Larv08 -> Larv13
- Larv09 -> Larv05
- Larv10 -> Larv14
- Larv11 -> Larv06
- Larv12 -> Fdrone
- Larv13 -> Larv10
- 
- Larv16 -> Worker
- Larv15 -> Larv02
- Worker -> Larv08
- Fdrone -> Larv09

Ended up rerunning MarkDuplicates without `${sampleId}` channels. 

- [ ] Delete list above in future commit

#### 4. Mark regions for realignment around indels 
```
qsub run_scripts/4_realignTargetCreator.sh
```
