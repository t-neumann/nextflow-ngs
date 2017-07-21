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

log.info " SICER peak calling pipeline "
log.info "============================="

masterChannel = Channel
         		.from(params.files)
         		.map { condition, list -> 
            		def files = list.collect{ file(it)} 
            		return tuple(condition, files) 
        		}
        		
SICERPath = "/groups/vbcf-ngs/bin/peaks/SICER_V1.1/SICER/SICER.sh"

process sicer {

	echo true

	tag { name }
	
	module params.bedtools
     
    input:
    set val(name), file(reads) from masterChannel
     
    output:
    file "${name}_SWEMBL.bed" into outChannel
 
    shell:
    '''
    bamToBed -i !{reads[0]} > chip.bed
    bamToBed -i !{reads[1]} > control.bed
    
    !{SICERPath} $PWD chip.bed control.bed . mm9 1 200 600 0.74 600 0.01
    
    echo $\'#chrom\tstart\tend\tname\tfold_change\tstrand\' > !{name}_SICER.bed
    
    sort -n -k 8 !{name}-W200-G600-islands-summary-FDR0.01 | awk \'BEGIN{FS="\t"; OFS="\t"} {print $1,$2,$3,NR,$7,"+"}\' >> !{name}_SICER.bed

    '''
}
 
 
workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
