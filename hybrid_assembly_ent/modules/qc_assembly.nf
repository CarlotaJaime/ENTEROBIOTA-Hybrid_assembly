process busco {
    tag "$meta.id"
    conda 'bioconda::busco=6.1.0'
    publishDir "${params.outdir}/busco", mode: 'copy'

    input:
    tuple val(meta), path(assembly)

    output:
    path("${meta.id}_busco"), emit: results

    script:
    """
    busco -i "${assembly}" -m genome -l ${params.busco_lineage} -o "${meta.id}_busco" -c ${task.cpus}
    """
}

process quast {
    tag "$meta.id"
    conda 'bioconda::quast=5.3.0'
    publishDir "${params.outdir}/quast", mode: 'copy'

    input:
    tuple val(meta), path(assembly)

    output:
    path("quast_report"), emit: results

    script:
    """
    quast ${assembly} --output-dir quast_report --threads ${task.cpus}
    """
}