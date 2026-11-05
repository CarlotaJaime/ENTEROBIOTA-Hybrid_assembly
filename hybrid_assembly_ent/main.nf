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
// Include subworkflows
include { PREPROCESSING } from './subworkflows/preprocessing'
include { HYBRID_ASSEMBLY } from './subworkflows/hybrid_assembly'
include { CHARACTERIZATION } from './subworkflows/characterization'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Define parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
params {
    input = null
    species = null
    outdir = 'results'
    adapters = null
    busco_lineage = 'enterococcus_odb12'
}

def genome_size_specie(String species) {
    if (species == 'Enterococcus_faecalis') {2990000}
    else if (species == 'Enterococcus_faecium') {3050000}
    else { error "ERROR: Species must be Enterococcus faecalis or Enterococcus faecium. ${species} not recognized" }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Define workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow {

    // Check inputs
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
    ch_input = channel.fromPath(params.input)
                        .splitCsv(header: true)
                        .map { row -> tuple( [id: row.sample.toString()],
                        file(row.illumina_r1),
                        file(row.illumina_r2),
                        file(row.ont), genome_size)}
    ch_illumina_raw = ch_input.map { meta, r1, r2, ont -> tuple(meta, r1, r2) }
    ch_ont_raw = ch_input.map { meta, r1, r2, ont -> tuple(meta, r1, r2) }

    PREPROCESSING ( ch_illumina_raw, ch_ont_raw)
    HYBRID_ASSEMBLY (PREPROCESSING.out.trimmed_illumina, PREPROCESSING.out.filtered_ont)
    CHARACTERIZATION (HYBRID_ASSEMBLY.out.assembly, PREPROCESSING.out.trimmed_illumina, PREPROCESSING.out.filtered_ont)
}