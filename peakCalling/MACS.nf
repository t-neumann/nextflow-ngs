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

log.info " MACS2 peak calling pipeline "
log.info "============================="

masterChannel = Channel
         		.from(params.files)
         		.map { condition, list -> 
            		def files = list.collect{ file(it)} 
            		return tuple(condition, files) 
        		}
        		
process macs {

	echo true

	tag { name }
	
	module params.macs
     
    input:
    set val(name), file(reads) from masterChannel
     
    output:
    file "${name}_MACS.bed" into outChannel
 
    shell:
    '''
   	macs2 callpeak -t !{reads[0]} -c !{reads[1]} -f AUTO -g mm -n !{name} --nomodel --extsize 150
   	
   	grep -v "^#" !{name}_peaks.xls | grep -v "fold_enrichment" | grep -v "^$" | \\
   	awk \'BEGIN{FS="\\t"; OFS="\\t"} {print $1, $2, $3, $10, $9, "+"}\' > !{name}_MACS.bed

    '''
}
 
 
workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
