/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    This subworkflow does hybrid assembly using 2 different
    strategies depending on the quality of the initial data:
        - Unicycler (short reads first)
        - Autocycler (long reads first)        
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { unicycler; autocycler; medaka; polypolish } from './../modules/assembly'
include { ont_quality } from './../modules/ont'
include { busco; quast } from './../modules/qc_assembly'

workflow HYBRID_ASSEMBLY {
    take:
    ch_illumina
    ch_ont

    main:
    // Calcule ONT quality metrics
    ont_quality(ch_ont)
    // Join quality metrics with the original ONT reads
    ch_ont_quality = ont_quality.out.metrics.join(ch_ont, by:0)
                                            .map{meta, metrics, ontfile -> 
                                            def values = metrics.txt.readLines().last().split('\t')
                                            def ont_depth = values[1] as double
                                            def ont_mean_q = values[2] as double
                                            def ont_n50 = values[3] as double
                                            def strategy // Select strategy for hybrid assembly
                                            if (
                                                ont_depth >= params.ont_min_depth &&
                                                ont_mean_q >= params.ont_min_mean_q &&
                                                ont_n50 >= params.ont_min_n50
                                            ) {
                                                strategy = 'autocycler'
                                            } else {
                                                strategy = 'unicycler'
                                            }
                                            tuple(meta, ontfile, strategy)

                                        }

    // Separate samples depending on the strategy selected
    ch_unicycler_ont = ch_ont_quality.filter { _meta, _ontfile, strategy -> strategy == 'unicycler'} // Added meta and ontfile in .filter operator but is not used (silenced with _)
                                    .map { meta, ontfile, _strategy -> tuple(meta, ontfile) } // Added strategy in .map operator but is not used (silenced with _)
    ch_autocycler_ont = ch_ont_quality.filter { _meta, _ontfile, strategy -> strategy == 'autocycler'} // Added meta and ontfile in .filter operator but is not used (silenced with _)
                                    .map { meta, ontfile, _strategy -> tuple(meta, ontfile) } // Added strategy in .map operator but is not used (silenced with _)

    // Unicycler strategy for assembly
    ch_unicycler_input = ch_illumina.join(ch_unicycler_ont, by: 0)
                                    .map { meta, r1, r2, ontfile -> tuple(meta, r1, r2, ontfile)} // Combine [meta, r1, r2] + [meta, onfile] to use Unicycler
    unicycler(ch_unicycler_input)

    // Autocycler stratgey for assembly
    ch_autocycler_input = ch_autocycler_ont.map { meta, ontfile -> tuple(meta, ontfile, meta.genome_size)}
    autocycler(ch_autocycler_input)
    ch_medaka_input = autocycler.out.assembly.join(ch_autocycler_ont, by: 0)
                                            .map { meta, assembly, ontfile -> tuple(meta, assembly, ontfile)}
    medaka(ch_medaka_input)
    ch_polypolish_input = medaka.out.assembly.join(ch_illumina, by:0)
                                            .map {meta, assembly, r1, r2 -> tuple(meta, assembly, r1, r2)}
    polypolish(ch_polypolish_input)

    // Save the assembly in a channel
    ch_final_assembly = unicycler.out.assembly.mix(polypolish.out.assembly) // No crec que sigui mix no? Només volem guardar el assembly final en un channel

    // Quality control of the assembly
    busco(ch_final_assembly)
    quast(ch_final_assembly)

    emit:
    assembly = ch_final_assembly // [meta, assembly_fasta]
}