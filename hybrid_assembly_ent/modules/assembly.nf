process unicycler {
    tag "$meta.id"
    conda 'bioconda::unicycler=0.5.1'
    publishDir "${params.outdir}/assemblies", mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2), path(ontfile)

    output:
    tuple val(meta), path("${meta.id}_final.fasta"), emit: assembly

    script:
    """
    unicycler -1 "${r1}" -2 "${r2}" -l "${ontfile}" -o "unicycler_output" --mode normal --threads ${task.cpus}
    """
}