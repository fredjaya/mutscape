#!/usr/bin/env nextflow

// Define global parameters 
params.scpath = "/scratch/Scape/fred/2008_manual" 
params.mode   = false
params.ref    = "${params.scpath}/GCF_003254395.2_Amel_HAv3.1_genomic.fna" 
params.outdir = "/scratch/Scape/fred/nf_tests"

if (params.mode == 'bamCh') {
/*
 * TEST: realignTargetCreator input channels 
 * ${bam} returns *.bam file from FilePair channel correctly 
 * ${sampleId} returns correctly
 */

    params.markedbam = "${params.scpath}/*.sam_sorted.bam_marked_dups.{bai,bam}"
    markedbam_ch = Channel.fromFilePairs(params.markedbam)

    process testBamCh {
    
        input:
            tuple val(sampleId), path(bamfiles) from markedbam_ch
    
        script:
            def bam = bamfiles.findAll{ it.toString() =~ /.bam$/ }.join('')
            """
            echo ${sampleId}
            echo ${bam}
            """ 
    }
}

