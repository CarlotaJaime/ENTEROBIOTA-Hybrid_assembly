#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Project: ENTEROBIOTA
    Process: Hybrid assembly using Illumina and ONT data
    Species: Enterococcus faecalis and Enterococcus faecium
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import modules and subworkflows
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { raw_fastqc; trimmomatic; trimmed_fastqc; raw_multiqc; trimmed_multiqc } from './modules/illumina'
include { nanoplot; filtlong; nanostats } from './modules/ont'
include { unicycler; autocycler; medaka; polipolish } from './modules/assembly'
include { map_ont; map_illumina } from './modules/mapping'
include { busco; quast } from './modules/qc_assembly'
include { plassembler } from './modules/plasmids'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Define parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
params {
    params.input = null
    params.species = null
    params.outdir = 'results'
    params.adapters = null
    params.run_mapping = true 
    params.run_busco = true
    params.run_quast = true
    params.busco_lineage = 'enterococcus_odb12' 
}

def genome_size_specie(String species) {
    switch (species) {
        case 'Enterococcus_faecalis': return 2990000
        case 'Enterococcus_faecium': return 3050000
        default:
            error "ERROR: Species must be Enterococcus faecalis or Enterococcus faecium. ${species} not recognized"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Define workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow {
    if (!params.input) {
        error "ERROR: No input provided. A samplesheet.csv is required as --input"
    }
    if (!params.species) {
        error "ERROR: No species provided. The species must be indicated using --species"
    }
    
    def genome_size = genome_size_specie(params.species) // Define genome size based on species provided

    /* Create a channel for the input file
        Samplesheet must contain 4 field (Sample, Illumina R1, Illumina R2 and ONT)
    */
    input_ch = channel.fromPath(params.input)
                        .splitCsv(header: true)
                        .map { row -> tuple( [id: row.sample.toString()],
                        file(row.illumina_r1),
                        file(row.illumina_r2),
                        file(row.ont), genome_size)}

    // Quality control Illumina
    raw_fastqc(samples_ch.map { meta, r1, r2, ont, gs -> tuple(meta, r1, r2)})
    raw_multiqc(illumina_fastqc.out.zip_reports.map { meta, report -> report }.collect())
    
    // Trimming Illumina
    trimmomatic(samples_ch.map { meta, r1, r2, u1, u2 -> tuple(meta, r1, r2)})
    trimmed_fastqc(trimmomatic.out.reads.map { meta, r1, r2, u1, u2 -> tuple(meta, r1, r2)})
    multiqc_trimmed(trimmed_fastqc.out.zip_reports.map { meta, report -> report }.collect())

    // Quality control ONT
    ont_qc_filter(samples_ch.map { meta, r1, r2, ontfile, gs -> tuple(meta, ontfile, gs)})

    // Hybrid assembly
    assembly_input = trimmomatic.out.reads.join(ont_qc_filter.out.filtered, by: 0).map { meta, r1, r2, u1, u2, ont_filtered -> tuple(meta, r1, r2, ont_filtered)}
    unicycler(assembly_input)

    // Assembly quality control
    if (params.run_busco) {
        busco(unicycler.out.assembly)
    }
    if (params.run_quast) {
        quast(unicycler.out.assembly.collect())
    }

    // Mapping
    if (params.run_mapping) {
        mapping_input = unicycler.out.assembly.join(ont_qc_filter.out.filtered, by: 0).join(trimmomatic.out.reads, by: 0).map { meta, assembly, ontfile, r1, r2, u1, u2 -> 
        tuple(meta, assembly, ontfile, r1, r2)
        }
        map_ont(mapping_input.map { meta, assembly, ontfile, r1, r2 -> tuple(meta, assembly, ontfile)})
        map_illumina(mapping_input.map { meta, assembly, ontfile, r1, r2 -> tuple(meta, assembly, r1, r2)})
    }

                        

}
