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

// Remove warnings by setting default to exact

if (!params.exact) {
	params.exact = "no"
}

log.info " Mappability track creation pipeline "
log.info "====================================="
log.info "genome file          : ${params.genome}"
log.info "kmer size            : ${params.kmerSize}"
log.info "chr sizes            : ${params.chrSizes}"
log.info "threads              : ${params.threads}"
log.info "exact                : ${params.exact}"

gemEnv = "/groups/zuber/zubarchive/USERS/tobias/bin/GEM-binaries-Linux-x86_64-core_i3-20130406-045632/bin/"

Channel
    .fromFilePairs( params.genome, size: 1 ) { file -> file.baseName }
    .ifEmpty { error "Cannot find genome fasta matching: ${params.genome}" }
    .set { genomeFile }
    
chrSizesFile = Channel
    .fromPath( params.chrSizes)
    .ifEmpty { error "Cannot find chr sizes matching: ${params.chrSizes}" }

process index {

    input:
    set val(name), file(genome) from genomeFile
    
    output:
    set name, "${name}.gem" into indexChannel
      
    script:
    """
    export PATH=$PATH:${gemEnv}
    
    gem-indexer -T ${params.threads} -c dna -i ${genome} -o ${name}
    """
}

process mappability {
    input:
    set val(name), file(index) from indexChannel
    
    output:
    set name, file(index), "${name}_kmer${params.kmerSize}.mappability" into mappabilityChannel
      
    script:
    
    if (params.exact != "no")
    	"""
    	export PATH=$PATH:${gemEnv}
    
    	gem-mappability -e 0 -T ${params.threads} -I ${index} -l ${params.kmerSize} -o ${name}_kmer${params.kmerSize}
    	"""
    else
    	"""
    	export PATH=$PATH:${gemEnv}
    
    	gem-mappability -T ${params.threads} -I ${index} -l ${params.kmerSize} -o ${name}_kmer${params.kmerSize}
    	"""
}

process toBigWig {

	module params.ucsc_kent
	
    input:
    file(chrSizesFile)
    set val(name), file(index), file(mappability) from mappabilityChannel
    
    output:
    file("${name}_kmer${params.kmerSize}.bw")
      
    script:
    
    """
    export PATH=$PATH:${gemEnv}
    
    gem-2-wig -I ${index} -i ${mappability} -o ${name}_kmer${params.kmerSize}
    wigToBigWig ${name}_kmer${params.kmerSize}.wig ${chrSizesFile} ${name}_kmer${params.kmerSize}.bw	
    """
}

workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
