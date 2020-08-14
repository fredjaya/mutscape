#!/usr/bin/env nextflow

// Define global parameters 
params.scpath    = "/scratch/Scape/fred/2008_manual" 
params.mode      = false
params.ref       = "${params.scpath}/GCF_003254395.2_Amel_HAv3.1_genomic.fna" 
params.fai       = "${params.ref}.fai"
params.dict      = "${params.scpath}/GCF_003254395.2_Amel_HAv3.1_genomic.dict"
params.markedbam = "${params.scpath}/*.sam_sorted.bam_marked_dups.{bam,bai}"
params.outdir    = "/scratch/Scape/fred/nf_tests"

markedbam_ch = Channel.fromFilePairs(params.markedbam)

if (params.mode == 'bamCh') {
/*
 * TEST: realignTargetCreator input channels 
 * ${bam} returns *.bam file from FilePair channel correctly 
 * ${sampleId} returns correctly
 */

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

if (params.mode == 'indelCh') {
/*
 * TEST: IndelRealigner channels input correctly
 * Sucessfully returns 16 tuples with
 * [sampleId, [sampleId.bam, sampleId.bai], sampleId.list] 
 */

    params.intervals = "${params.scpath}/*_target_intervals.list" 
    intervals_ch = Channel
                    .fromPath(params.intervals)
                    .map { file -> 
                            def sampleId = file.name.toString().tokenize('_').get(0)
                            return tuple(sampleId, file)
                    } 
    
    // Combine {bam, bai} and {intervals} channels
    // Aim: [sampleId, [sampleId.bam, sampleId.bai], sampleId.list]
    
    realigntargets_ch = markedbam_ch.join(intervals_ch)
   
    cpus = 2                     
    memory = 14.GB               
    time = '5h'                  
    publishDir "${params.outdir}"
    tag "$sampleId"              
     
    process testIndelCh {
    
        input:
            path ref  from params.ref                            
            path fai  from params.fai                            
            path dict from params.dict                           
            tuple val(sampleId), path(bamfiles), path(intervals) from realigntargets_ch
    
        script:
        def bam = bamfiles.findAll{ it.toString() =~ /.bam$/ }.join('')
        """
        module load gatk/3.8.1
    
        gatk -T IndelRealigner \
             -R ${ref} \
             -I ${bam} \
             -targetIntervals ${intervals}
             -o ${sampleId}_realigned.bam 
        """
    }

}


if (params.mode == 'hapCallSingle') {

/*
 * TEST: HaploptypeCaller works on a single bee sample
 */
    
    params.realignedbam = "${params.scpath}/Fdrone_realigned.{bam,bai}"
    realignedbam_ch = Channel.fromFilePairs(params.realignedbam)   

    process testHapCallSingle {
  
        executor = 'pbs'
        clusterOptions = '-P RDS-FSC-Scape-RW'
        cpus = 4
        memory = 16.GB
        time = '12h'
        publishDir "$params.outdir"
        tag "$sampleId"
 
        input:
            path ref  from params.ref 
            path fai  from params.fai 
            path dict from params.dict
            tuple val(sampleId), path(bamfiles) from realignedbam_ch

        output:
            tuple val(sampleId), path("${sampleId}_unrecal_conf_sites.vcf")

        script:
        def bam = bamfiles.findAll{ it.toString() =~ /.bam$/ }.join('')
        """
        module load gatk/3.8.1

        gatk -T HaplotypeCaller \
             -R ${ref} \
             -I ${bam} \
             --genotyping_mode DISCOVERY \
             --output_mode EMIT_ALL_CONFIDENT_SITES \
             -stand_call_conf 30 \
             -o ${sampleId}_unrecal_conf_sites.vcf \
             -nt 1 -nct 24 
        """
    }
}
