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
//params.reads  = "/project/Scape/Trimmomatic/Paired/*_{R1,R2}_paired.fastq.gz"
params.reads  = "/scratch/Scape/fred/nf_test/*_{R1,R2}_paired.fastq.gz"
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

process bwa {
  
  input: 
    path ref from params.ref 
    tuple val(sampleId), path(reads) from reads_ch
    
  output:
    tuple val(sampleId), path("${sampleId}_pe_sorted.bam") into markdups_ch 

  script:
    readGroup = "@RG\\tID:${sampleId}\\tSM:${sampleId}\\tPL:ILLUMINA\\tLB:lib1\\tPU:unit1"
    
    """
    module load bwa/0.7.17
    module load samtools/1.9
 
    bwa mem -R \"${readGroup}\" -t 24 ${ref} ${reads} | \
    samtools sort -@24 -o ${sampleId}_pe_sorted.bam -O bam
    """
}


process markDuplicates {
  
  input: 
    tuple val(sampleId), path(bam) from markdups_ch

  output:
    tuple \
      val(sampleId) \
      path("${sampleId}_marked_dups.bam") \
      path("${sampleId}_marked_dups.bai") into marked_bam_ch    
    tuple val(sampleId), path("${sampleId}_marked_dup_metrics.txt")
     
  script:
    """ 
    module load picard/2.7.1
  
    picard MarkDuplicates \
           I=${bam} \
           O=${sampleId}_marked_dups.bam \
           M=${sampleId}_marked_dup_metrics.txt \
           CREATE_INDEX=true
    """
}

process realignerTargetCreator {
  
  input: 
    path ref from params.ref 
    tuple val(sampleId), path(bam), path(bai) from marked_bam_ch 

  output:
    tuple \
      val(sampleId) \
      path("${sampleId}"_target_intervals.list) into intervals_ch 

  script:
    """
    module load gatk/3.8.1

    gatk -T RealignerTargerCreator \
         -R ${ref} \
         -I ${bam} \
         -o ${sampleId}_target_intervals.list
    """

}

process indelRealigner {
  
  input: 
    path ref from params.ref
    tuple val(sampleId), path(bam), path(bai) from marked_bam_ch
    tuple val(sampleId), path(intervals) from intervals_ch

  output:
    tuple \
      val(sampleId) \
      path("${sampleId}_realigned_reads.bam") \
      path("${sampleId}_realigned_reads.bai") into realigned_bam_ch
       
  script:
    """
    gatk -T IndelRealigner \
         -R ${ref} \
         -I ${bam} \
         -targetIntervals ${intervals} \
         -o ${sampleId}_realigned_reads.bam 
    """

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
