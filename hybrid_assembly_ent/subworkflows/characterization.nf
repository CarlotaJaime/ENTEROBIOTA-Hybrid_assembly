/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    This subworkflow does the mapping of the reads against
    the hybrid assembly generated
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { map_ont; map_illumina } from '../modules/mapping'

workflow CHARACTERIZATION {
    take:
    ch_assembly
    ch_illumina
    ch_ont

    main:
    // For ONT data
    ch_map_ont_input = ch_assembly.join(ch_ont)
    map_ont(ch_map_ont_input)

    // For Illumina data
    ch_map_illumina_input = ch_assembly.join(ch_illumina)
    map_illumina(ch_map_illumina_input)
}