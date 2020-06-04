# mutscape  

Proof of concept pipeline to estimate _S. capensis_ mutation rates

## Aims
- To familiarise self with the processes and filtering involved
- Produce a filtered .vcf file of `Larv01` for comparison with `Fdrone` 

## Map trimmed reads to reference

### Indexing  

#### Download the reference genome for _Apis mellifera_:
```
rsync --copy-links --times --verbose \ 
rsync://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/254/395/GCA_003254395.2_Amel_HAv3.1/GCA_003254395.2_Amel_HAv3.1_genomic.fna.gz \
/scratch/Scape/fred/ # Artemis HPC
```

Output: `GCA_003254395.2_Amel_HAv3.1_genomic.fna.gz`

#### 1. Index reference:
```
qsub scripts/1_bwa_index.sh
```

Input: `GCA_003254395.2_Amel_HAv3.1_genomic.fna.gz`
Outputs:
```
225251116 Jun  4 17:09 GCA_003254395.2_Amel_HAv3.1_genomic.fna.gz.bwt
112625560 Jun  4 17:10 GCA_003254395.2_Amel_HAv3.1_genomic.fna.gz.sa
 56312755 Jun  4 17:09 GCA_003254395.2_Amel_HAv3.1_genomic.fna.gz.pac
    21064 Jun  4 17:09 GCA_003254395.2_Amel_HAv3.1_genomic.fna.gz.ann
      859 Jun  4 17:09 GCA_003254395.2_Amel_HAv3.1_genomic.fna.gz.amb
```

### 2. BWA-MEM
Suggested use for Illumina paired-end reads longer than ~70:
```
qsub scripts/2_bwa_mem.sh
```
Inputs:
```
GCA_003254395.2_Amel_HAv3.1_genomic.fna.gz.bwt
GCA_003254395.2_Amel_HAv3.1_genomic.fna.gz.sa
GCA_003254395.2_Amel_HAv3.1_genomic.fna.gz.pac
GCA_003254395.2_Amel_HAv3.1_genomic.fna.gz.ann
GCA_003254395.2_Amel_HAv3.1_genomic.fna.gz.amb
```

Output: `Larv01-pe.sam`
