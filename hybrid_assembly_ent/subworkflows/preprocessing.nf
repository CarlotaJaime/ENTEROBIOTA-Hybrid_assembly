/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    This subworkflow does:
        - The initial quality control
        - Filtering of Illumina and ONT data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { raw_fastqc; trimmomatic; trimmed_fastqc; raw_multiqc; trimmed_multiqc } from './../modules/illumina'
include { nanoplot; filtlong } from './../modules/ont'

workflow PREPROCESSING {
    take:
    ch_illumina // [meta, r1, r2]
    ch_ont // [meta, ontfile]

    main:
    // For ONT data
    filtlong(ch_ont)
    ch_nanoplot_input = ch_ont.join(filtlong.out.filtered_ont)
    nanoplot(ch_nanoplot_input)

    // For Illumina data
    raw_fastqc(ch_illumina)
    raw_multiqc(raw_fastqc.out.zip_reports.map{ _meta, file -> file }.collect()) // Added meta in .map operator but is not used (silenced with _)
    trimmomatic(ch_illumina)
    ch_trimmed_paired = trimmomatic.out.reads.map { meta, r1_paired, r2_paired, _r1_unpaired, _r2_unpaired -> [meta, r1_paired, r2_paired] } // Added r1_unpaired and r2_unpaired in .map operator but is not used (silenced with _)
    trimmed_fastqc(ch_trimmed_paired)
    trimmed_multiqc(trimmed_fastqc.out.zip_reports.map{ _meta, file -> file }.collect()) // Added meta in .map operator but is not used (silenced with _)

    emit:
    trimmed_illumina = ch_trimmed_paired
    filtered_ont = filtlong.out.filtered_ont

}