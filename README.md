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

## Citation
### Quality control
Wingett, S. W., & Andrews, S. (2018). FastQ Screen: A tool for multi-genome mapping and quality control. F1000Research, 7, 1338.
Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19), 3047-3048.
De Coster, W., & Rademakers, R. (2023). NanoPack2: population-scale evaluation of long-read sequencing data. `Bioinformatics`, 39(5), btad311.
Simão, F. A., Waterhouse, R. M., Ioannidis, P., Kriventseva, E. V., & Zdobnov, E. M. (2015). BUSCO: assessing genome assembly and annotation completeness with single-copy orthologs. `Bioinformatics`, 31(19), 3210-3212.
Gurevich, A., Saveliev, V., Vyahhi, N., & Tesler, G. (2013). QUAST: quality assessment tool for genome assemblies. `Bioinformatics`, 29(8), 1072-1075.

### Trimming
Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: a flexible trimmer for Illumina sequence data. `Bioinformatics`, 30(15), 2114-2120.
Furqon, A. D. C., Wick, R. R., Roberts, L. W., & Hall, M. B. (2026). Evaluating Nanopore read trimming and quality filtering in bacterial genomics. `bioRxiv`, 2026-10.

### Assembly
Wick, R. R., Judd, L. M., Gorrie, C. L., & Holt, K. E. (2017). Unicycler: resolving bacterial genome assemblies from short and long sequencing reads. `PLoS computational biology`, 13(6), e1005595.

### Mapping
Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. `Bioinformatics`, 34(18), 3094-3100.
Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., ... & Li, H. (2021). Twelve years of SAMtools and BCFtools. `Gigascience`, 10(2), giab008.
Li, H., & Durbin, R. (2009). Fast and accurate short read alignment with Burrows–Wheeler transform. `Bioinformatics`, 25(14), 1754-1760.