# ENTEROBIOTA-Hybrid_assembly
This repository contains the codes used in the hybrid assembly of the ENTEROBIOTA project.
This project aims to characterize the genomes of clinical and probiotic isolates of _Enterococcus faecalis_ and _Enterococcus faecium_ for identification of virulence factors and antimicrobial resistance genes, while studying the mobile genetic elements of those genomes for acquisition and dissemination analysis of the virulence and resistance genes.

The following workflow is divided into 3 main subworkflows:

## Quality control and processing
For the quality control of Illumina data, FastQC and MultiQC are used for quality control of raw data, then Trimmomatic is used for the trimming of the raw reads, and a quality control of trimmed reads is performed also with FastQC and MultiQC.
For Nanopore, a quality control of raw data is performed with NanoPlot, the sequences are trimmed with Filtlong and a following quality control of filtered data is also performed by NanoPlot.

## Hybrid assembly
The workflow has integrated a decision block that depending on the quality of the reads, there are two possibilities of assembly:

### Unicycler
This tool is used when the quality of ONT data is insufficient. The assembly is performed with the short-reads-first strategy using Illumina reads to performe the initial assembly that is polished with Nanopore data.
### Autocycler
This tool is used when the quality of ONT data is sufficient. The assembly is performed with the long-reads-first strategy using ONT reads to perform the initial assembly and Polypolish is used to polish the ONT assembly with Illumina reads. Additionally, Medaka is used to polish the ONT data before the initial assembly.

## Characterization
The quality control of the assembly is performed with BUSCO and QUAST, then a mapping is performed for ONT and Illumina reads against the hybrid assembly, and the error rate of the assembly is calculated.

## Execute the workflow
In order to execute the workflow the `samplesheet.csv` must have this structure following the structure of a CSV file:

sample,illumina_r1,illumina_2,ont

For example:

```
sample1_name,sample1_name_R1.fastq.gz,sample1_name_R2.fastq.gz,sample1_name_ONT.fastq.gz
sample2_name,sample2_name_R1.fastq.gz,sample2_name_R2.fastq.gz,sample2_name_ONT.fastq.gz
sample3_name,sample3_name_R1.fastq.gz,sample3_name_R2.fastq.gz,sample3_name_ONT.fastq.gz
```

And to execute the workflow the following command is the basic structure to execute in the terminal:

```
nextflow run main.nf \
    --input samplesheet.csv
    --species Enterococcus_faecalis \
    --outdir results \
    --adapters adapters.fa
    -with-conda
```
