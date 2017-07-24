#!/usr/bin/env nextflow

/*
* MIT License
*
* Copyright (c) 2017 Tobias Neumann
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in all
* copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*/

log.info " GATK basic variant calling pipeline on RNA-Seq data "
log.info "====================================================="
log.info "bams              : ${params.bams}"
log.info "genome            : ${params.genome}"

Channel
    .fromFilePairs( params.bams, size: 1 )
    .ifEmpty { error "Cannot find any reads matching: ${params.bams}" }
    .set { bams }
    
if (params.genome == "human") {
	genome = ""
} else {
	genome = "/groups/zuber/zubarchive/USERS/tobias/mareike/rnaseq/GATK/common/Mus_musculus.GRCm38.dna.primary_assembly.chr.fa"
}

process preprocess {

	tag { name }
	
	module params.samtools
	module params.picard
     
    input:
    set val(name), file(bam) from bams
     
    output:
    set val(name), "${name}.ro.dedup.bam" into preprocessChannel
    
    shell:
    ''' 
    shopt -s expand_aliases
	
	java -jar /biosw/debian7-x86_64/picard-tools/2.6.0/picard.jar AddOrReplaceReadGroups INPUT=!{bam} OUTPUT=!{name}.rg.bam RGID=!{name} RGLB=!{name} RGPL="Illumina" RGPU=!{name} RGSM=!{name}

	java -jar /biosw/debian7-x86_64/picard-tools/2.6.0/picard.jar ReorderSam INPUT=!{name}.rg.bam OUTPUT=!{name}.ro.bam REFERENCE=!{genome} ALLOW_CONTIG_LENGTH_DISCORDANCE=false

	samtools sort -@ 10 -o !{name}.ro.sorted.bam !{name}.ro.bam

	java -jar /biosw/debian7-x86_64/picard-tools/2.6.0/picard.jar MarkDuplicates INPUT=!{name}.ro.sorted.bam OUTPUT=!{name}.ro.dedup.bam METRICS_FILE=!{name}.ro.dedup.metrics.txt

	java -jar /biosw/debian7-x86_64/picard-tools/2.6.0/picard.jar ValidateSamFile I=!{name}.ro.dedup.bam O=!{name}.ro.dedup.val.txt VALIDATION_STRINGENCY=STRICT IGNORE="MISSING_TAG_NM"

	samtools index !{name}.ro.dedup.bam

	module unload picard-tools
	
	
    '''
}
 
 
workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
