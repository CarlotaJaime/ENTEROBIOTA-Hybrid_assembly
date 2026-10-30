process map_ont {
    tag "$meta.id"
    conda 'bioconda::minimap2=2.31 bioconda::samtools=1.24'
    publishDir "${params.outdir}/mapping", mode: 'copy'

    input:
    tuple val(meta), path(assembly), path(ont)

    output:
    path("${meta.id}_ONT.bam"), emit: bam
    path("${meta.id}_ONT.bam.bai"), emit: index
    path("${meta.id}_ONT_flagstat.txt"), emit: flagstat
    path("${meta.id}_ONT_depth.txt"), emit: depth

    script:
    """
    minimap2 -t ${task.cpus} -ax map-ont "${assembly}" "${ont}" | samtools sort -@ ${task.cpus} -o "${meta.id}_ONT.bam"
	samtools index "${meta.id}_ONT.bam"
	samtools flagstat "${meta.id}_ONT.bam" > "${meta.id}_ONT_flagstat.txt"
	samtools depth -a "${meta.id}_ONT.bam" > "${meta.id}_ONT_depths.txt"
    """
}

process map_illumina {
    tag "$meta.id"
    conda 'bioconda::bwa=0.7.19 bioconda::samtools=1.24'
    publishDir "${params.outdir}/mapping", mode: 'copy'

    input:
    tuple val(meta), path(assembly), path(r1), path(r2)

    output:
    path("${meta.id}_Illumina.bam"), emit: bam
    path("${meta.id}_Illumina.bam.bai"), emit: index
    path("${meta.id}_Illumina_flagstat.txt"), emit: flagstat
    path("${meta.id}_Illumina_depth.txt"), emit: depth

    script:
    """
    bwa index "${assembly}"
	bwa mem -t ${task.cpus} "${assembly}" "${r1}" "${r2}" | samtools sort -@ ${task.cpus} -o "${meta.id}_Illumina.bam"
	samtools index "${meta.id}_Illumina.bam"
	samtools flagstat "${meta.id}_Illumina.bam" > "${meta.id}_Illumina_flagstat.txt"
	samtools depth -a "${meta.id}_Illumina.bam" > "${meta.id}_illumina_depths.txt"
    """
}