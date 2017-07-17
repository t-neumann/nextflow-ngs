Nextflow-NGS
============
[![nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.24.0-brightgreen.svg)](http://nextflow.io)

Generic pipelines for next-generation sequencing data analysis using [Nextflow](https://www.nextflow.io/). These pipelines comprise robust and reasonable tool selection for various
analysis steps.

Preprocessing
-------------
* **Adapter trimming** - Standard adapter trimming and QC using [trim_galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)

Genome tools
------------
* **Mappability** - Genome-wide mappability assessment using the [GEM-library](http://algorithms.cnag.cat/wiki/The_GEM_library)
