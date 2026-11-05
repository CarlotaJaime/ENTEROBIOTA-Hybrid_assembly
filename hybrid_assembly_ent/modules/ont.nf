process nanoplot {
    tag "$meta.id"
    conda 'bioconda::nanoplot=1.48.0'
    publishDir "${params.outdir}/nanopore/nanoplot", mode: 'copy'

    input:
    tuple val(meta), path(ontfile), path(filtered_ont)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.png"), emit: png
    tuple val(meta), path("*.txt"), emit: txt
    path("nanoplot_raw"), emit: raw_qc 
    path("nanoplot_filtered"), emit: filtered_qc

    script:
    """
    # Quality control of unprocessed ONT data
    NanoPlot --fastq "${ontfile}" --outdir nanoplot_raw --tsv_stats

    # Quality control of filtered ONT data
    NanoPlot --fastq "${filtered_ont}" --outdir nanoplot_filtered --tsv_stats   
    """
}

process filtlong {
    tag "$meta.id"
    conda 'bioconda::filtlong=0.3.1'
    publishDir "${params.outdir}/nanopore/filtered", mode: 'copy'

    input:
    tuple val(meta), path(ontfile)

    output:
    tuple val(meta), path("${meta.id}_filtered.fastq.gz"), emit: filtered_ont

    script:
    """
     # Filtering ONT data
    filtlong --min_length 1000 --keep_percent 90 --mean_q_weigth 9 "${ontfile}" | gzip > "${meta.id}_filtered.fastq.gz"
    """
}

process ont_quality {
   tag "$meta.id"
    conda 'bioconda::nanoplot=1.41.0'

    input:
    tuple val(meta), path(ontfile)

    output:
    tuple val(meta), path("${meta.id}_ont_quality.tsv"), emit: quality_ont_values

    script:
    def genome_size_py = 0
    if (params.species == 'Enterococcus_faecalis') {
        genome_size_py = 2990000
    } else if (params.species == 'Enterococcus_faecium'){
        genome_size_py = 3050000
    }

    """
    NanoPlot \
        --fastq ${ontfile} \
        --outdir nanoplot \
        --threads ${task.cpus} \
        --N50 \
        --no_static
        
    # IMPORTANT: Before continuing check if the names of the columns in NanoStats.txt are:
    # - Number of bases
    # - Mean read quality
    # - Read length N50

    # Run python script to determine ONT quality for assembly decision
    # Define specie for later depth calculation

    python3 ont_quality.py "${meta.id}" "${genome_size_py}"
    """ 
}