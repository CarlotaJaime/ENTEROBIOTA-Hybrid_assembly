process raw_fastqc {
    tag "$meta.id"
    conda 'bioconda::fastqc=0.13.0' // Conda or a container? Docker or Conda?
    publishDir "${params.outdir}/fastqc/raw" , mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("*_fastqc.zip"), emit: zip_reports
    tuple val(meta), path("*_fastqc.html"), emit: html_reports

    script:
    """
    fastqc "${r1}" "${r2}" --threads ${task.cpus}
    """

}

process raw_multiqc {
    tag "multiqc_raw"
    conda 'bioconda::multiqc=1.35.0'
    publishDir "${params.outdir}/multiqc/raw" , mode: 'copy'

    input:
    path(reports) 

    output:
    path("multiqc_report.html"), emit: final_report

    script:
    """
	multiqc ${reports} --filename multiqc_report.html
    """

}

process trimmomatic {
    tag "$meta.id"
    conda 'bioconda::trimmomatic=0.41'
    publishDir "${params.outdir}/trimmomatic" , mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}_1_paired.fastq.gz"), path("${meta.id}_2_paired.fastq.gz"), path("${meta.id}_1_unpaired.fastq.gz"), path("${meta.id}_2_unpaired.fastq.gz"), emit: reads
    
    script:
    def adapter = params.adapters ?
        "ILLUMINACLIP:${params.adapters}:2:30:10" : ""
    """
    trimmomatic PE \
        -threads ${task.cpus} \
		"${r1}" "${r2}" \
		"${meta.id}_1_paired.fastq.gz" \
        "${meta.id}_1_unpaired.fastq.gz" \
        "${meta.id}_2_paired.fastq.gz" \
        "${meta.id}_2_unpaired.fastq.gz"\
        ${adapter}\
        LEADING:3 \
		TRAILING:3 \
		SLIDINGWINDOW:5:25 \
		MINLEN:60 \
    """

}

process trimmed_fastqc {
    tag "$meta.id"
    conda 'bioconda::fastqc=0.13.0' // Conda or a container? Docker or Conda?
    publishDir "${params.outdir}/fastqc/trimmed" , mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("*_fastqc.zip"), emit: zip_reports
    tuple val(meta), path("*_fastqc.html"), emit: html_reports

    script:
    """
    fastqc "${r1}" "${r2}" --threads ${task.cpus}
    """

}

process trimmed_multiqc {
    tag "multiqc_trimmed"
    conda 'bioconda::multiqc=1.35.0'
    publishDir "${params.outdir}/multiqc/trimmed" , mode: 'copy'

    input:
    path(reports) 

    output:
    path("multiqc_report.html"), emit: final_report

    script:
    """
	multiqc ${reports} --filename multiqc_report.html
    """

}