/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    This subworkflow does hybrid assembly using 2 different
    strategies depending on the quality of the initial data:
        - Unicycler (short reads first)
        - Autocycler (long reads first)        
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { unicycler } from './../modules/assembly'
include { busco; quast } from './../modules/qc_assembly'

workflow HYBRID_ASSEMBLY {
    take:
    ch_illumina
    ch_ont

    main:
    ch_unicycler_input = ch_illumina.join(ch_ont) // Combine [meta, r1, r2] + [meta, onfile] to use Unicycler
    unicycler(ch_unicycler_input)

    // Quality control of the assembly
    busco(unicycler.out.assembly)
    quast(unicycler.out.assembly)

    emit:
    assembly = unicycler.out.assembly // [meta, assembly_fasta]
}