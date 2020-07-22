#!/usr/bin/env nextflow

/*
 * M U T S C A P E 
 *
 * Pipeline to estimate mutation rates in cape honey bees
 * 
 * Fred Jaya
 *
 */

// PIPELINE SETTINGS

//// Define default parameters

params.ref    = "/scratch/Scape/fred/GCF_003254395.2_Amel_HAv3.1_genomic.fna"
params.reads  = "/project/Scape/Trimmomatic/Paired/*_{R1,R2}_paired.fastq.gz"
params.outdir = "$baseDir"

log.info """\
M U T S C A P E
===============
ref    : ${params.ref}
reads  : ${params.reads}
outdir : ${params.outdir}
"""

// Define read channels

reads_ch = Channel.fromFilePairs(params.reads)

// 1. PREPARE SEQUENCE DATA 

//// Map trimmed reads (paired) to reference genome

process '1_bwa_mem' {
  
  input: 
    path ref from params.ref 
    tuple val(sampleId), path(reads) from reads_ch
    
  output:
    tuple val(sampleId), path("${sampleId}_pe_sorted.bam" into markdups_ch 

  script:
    readGroup = "@RG\\tID:${sampleId}\\tSM:${sampleId}\\tPL:ILLUMINA\\tLB:lib1\\tPU:unit1"
    
    '''
    module load bwa/0.7.17
    module load samtools/1.9
 
    bwa mem -R \"${readGroup}\" -t 24 ${ref} ${reads} | \
    samtools sort -@24 -o ${sampleId}_pe_sorted.bam -O bam
    ''' 

}


/*
process {
  
  input: 


  output:


  script:
    """
    
    """

}
*/
