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

log.info " MACS2/SICER initiation site and zone peak calling pipeline "
log.info "============================================================"

masterChannel = Channel
         		.from(params.files)
         		.map { condition, list -> 
            		def files = list.collect{ file(it)} 
            		return tuple(condition, files) 
        		}
        		
masterChannel.into{sicerInputChannel; macsInputChannel}

        		        		
SICERPath = "/groups/vbcf-ngs/bin/peaks/SICER_V1.1/SICER/SICER.sh"

process sicer {

	echo true

	tag { name }
	
	module params.bedtools
	module params.python
     
    input:
    set val(name), file(reads) from sicerInputChannel
     
    output:
    set val(name), "${name}_SICER.bed" into sicerChannel
 
    shell:
    '''
    bamToBed -i !{reads[0]} > chip.bed
    bamToBed -i !{reads[1]} > control.bed
    
    !{SICERPath} $PWD chip.bed control.bed . mm9 1 200 600 0.74 600 0.01
    
    echo $\'#chrom\\tstart\\tend\\tname\\tfold_change\\tstrand\' > !{name}_SICER.bed
    
    sort -n -k 8 chip-W200-G600-islands-summary-FDR0.01 | awk \'BEGIN{FS="\\t"; OFS="\\t"} {print $1,$2,$3,NR,$7,"+"}\' >> !{name}_SICER.bed

    '''
}

process macs {

	echo true

	tag { name }
	
	module params.macs
     
    input:
    set val(name), file(reads) from macsInputChannel
     
    output:
    set val(name), "${name}_MACS.bed" into macsChannel
 
    shell:
    '''
   	macs2 callpeak -t !{reads[0]} -c !{reads[1]} -f AUTO -g mm -n !{name} --nomodel --extsize 150
   	
   	grep -v "^#" !{name}_peaks.xls | grep -v "fold_enrichment" | grep -v "^$" | \\
   	awk \'BEGIN{FS="\\t"; OFS="\\t"} {print $1, $2, $3, $10, $9, "+"}\' > !{name}_MACS.bed

    '''
}

process zone_site_calling {

echo true

	tag { name }
	
	module params.bedtools
     
    input:
    set val(name), file(peaks) from sicerChannel.phase(macsChannel).map { left, right -> tuple(left[0], [left[1], right[1]]) }
     
    output:
    file "${name}_IS.bed" into outISChannel
    file "${name}_IZ.bed" into outIZChannel
 
    shell:
    '''
    bedtools intersect -a !{peaks[0]} -b !{peaks[1]} -u -wa > !{name}_IZ.bed
    bedtools intersect -b !{peaks[0]} -a !{peaks[1]} -u -wa > !{name}_IS.bed

    '''
}

 
workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
