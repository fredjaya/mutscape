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

#### 4. Mark regions for realignment around indels 
```
qsub run_scripts/4_realignerTargetCreator.sh
```
Added tests here for inputting bam/bai files.

#### 5. Realign regions around indels
```
qsub run_scripts/5_indelRealigner
```

### 6. Call confident sites for BSQR
```
qsub run_scripts/6_haplotypeCallerUnrecal
```
