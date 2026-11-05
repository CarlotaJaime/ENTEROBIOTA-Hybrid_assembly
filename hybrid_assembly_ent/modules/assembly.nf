process unicycler {
    tag "$meta.id"
    conda 'bioconda::unicycler=0.5.1'
    publishDir "${params.outdir}/unicycler", mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2), path(ontfile)

    output:
    tuple val(meta), path("${meta.id}_unicycler_final.fasta"), emit: assembly

    script:
    """
    unicycler -1 "${r1}" -2 "${r2}" -l "${ontfile}" -o "unicycler_output" --mode normal --threads ${task.cpus}
    cp unicycler_output/assembly.fasta "${meta.id}_final.fasta" 
    """
}

process autocycler {
    tag "$meta.id"
    conda 'bioconda::autocycler=0.8.0'
    publishDir "${params.outdir}/autocycler", mode: 'copy'

    input:
    tuple val(meta),  path(ontfile), val(genome_size)

    output:
    tuple val(meta), path("${meta.id}_autocycler_final.fasta"), emit: assembly

    script:
    """
    # Subsample ONT reads
    autocycler subsample \
		--count 4 \
		--reads "${ontfile}" \
		--out_dir subsampled_reads \
		--genome_size "${genome_size}" \
		--min_read_depth 50

    # Generate assemblies with different assemblers
    mkdir -p assemblies
    for assembler in raven myloasm miniasm flye metamdbg necat nextdenovo plassembler canu; do
        for readfile in subsampled_reads/*.fastq; do
            if [[ -f "\${readfile}" ]]; then
                id=\$(basename "\${readfile}" .fastq)
                autocycler helper "\${assembler}" \
					--reads "\${readfile}"\
					--out_prefix "assemblies/\${assembler}_\${id}" \
					--genome_size "${genome_size}" \
                    --threads ${task.cpus}
            fi
        done
    done

    autocycler compress -i assemblies -a autocycler_out -t ${task.cpus}
	autocycler cluster -a autocycler_out

	for cluster in autocyler_out/clustering/qc_pass/cluster_*; do
		autocycler trim -c "\${cluster}"
		autocycler resolve -c "\${cluster}"
	done

	autocycler combine \
		-a autocycler_out \
		-i autocyler_out/clustering/qc_pass/cluster_*/5_final.gfa \
        -r "{ontfile}" \
        -t ${task.cpus}
    cp autocycler_out/concensus_assembly.fasta "${meta.id}_autocycler.fasta"
    """
}

process medaka {
    tag "$meta.id"
    conda 'bioconda::medaka=2.2.2'
    publishDir "${params.outdir}/autocycler", mode: 'copy'

    input:
    tuple val(meta),  path(ontfile), val(genome_size)

    output:
    tuple val(meta), path("${meta.id}_autocycler_final.fasta"), emit: assembly

    script:
    """
    medaka_consensus \
		-i "${ontfile}" \
		-d "${assembly}" \
		-o medaka_output \
        -t ${task.cpus}
    cp medaka_output/concensus.fasta "${meta.id}_medaka.fasta"        
    """
}

process polypolish {
    tag "$meta.id"
    conda 'bioconda::polypolish=0.7.1'
    publishDir "${params.outdir}/polypolish", mode: 'copy'

    input:
    tuple val(meta), path(assembly), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}_final.fasta"), emit: assembly

    script:
    """
    ## Index assembly with BWA
	bwa index "${assembly}"
	# Map Illumina reads
	bwa mem -a "${assembly}" "${r1}" > alignments_R1.sam
	bwa mem -a "${assembly}" "${r2}" > alignments_R2.sam

	# Filter alignments
	polypolish filter \
		--in1 alignments_R1.sam \
		--in2 alignments_R2.sam \
		--out1 filtered_R1.sam \
		--out2 filtered_R2.sam

	# Final polish step
	polypolish polish \
		"${assembly}" \
		filtered_R1.sam \
		filtered_R2.sam > "${meta.id}_final.fasta"
    """
}