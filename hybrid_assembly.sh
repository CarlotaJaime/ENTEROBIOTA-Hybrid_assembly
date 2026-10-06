#!/usr/bin/env/bash

#########################################################################
# HYBRID ASSEMBLY DRAFT
#########################################################################

# Project ENTEROBIOTA
# Hybrid assembly with Illumina short-reads and Oxford Nanopore Technology long-reads sequencing data
# This script is a draft for hybrid assembly and characterization of isolates selected for animal model analyses

# Last modification: 05/10/2026 16:01h

# Set up configurations
set -euo pipefail
shopt -s nullglob

# Help message

Help()
{
echo "ENTEROBIOTA - Hybrid assembly draft"
echo ""
echo "Required flags:"
echo " -i Input directory containing Illumina and ONT files"
echo " -s Species (Enterococcus_faecalis or Enterococcus_faecium)"
echo ""
echo "Optional flags:"
echo " -h Print this help message and exit"
echo ""
}

# Define arguments and parameters

while getopts i:s:h option; do
	case "${option}" in
		i) INPUT_DIR=${OPTARG} ;;
		s) SPECIES=${OPTARG} ;;
		h) Help; exit 0 ;;
		*) Help; exit 1 ;;
	esac
done

# Logfile generation
logstamp=$(date '+%Y%m%d%H%M%S') # The timestamp is the name of the logfile
logfile="hybrid_assembly_${logstamp}.log"
exec > >(tee -a "${logfile}") 2>&1

echo "-------------------------------------------------------------"
echo "ENTEROBIOTA project -Hybrid assembly and characterization"
echo "Log file: ${logfile}"
echo "Date: $(date)"
echo "Input directory: ${INPUT_DIR}"
echo "Species: ${SPECIES}"
echo "--------------------------------------------------------------"

# Set up environment
eval "$(conda shell.bash hook)"
WD=${INPUT_DIR}; cd ${WD}


# Check if the arguments are given

if [[ -z "${INPUT_DIR}" && -z "${SPECIES}" ]]; then # in case both are missing
	echo -e "ERROR: -i and -s are required"; Help; exit 1;
	elif [[ -z "${INPUT_DIR}" ]]; then # only -i is missing
		echo -e "ERROR: input directory is missing"; Help; exit 1;
	elif [[ -z "${SPECIES}" ]]; then # only -s is missing
		echo -e "ERROR: species is missing"; Help; exit 1;
fi

# Define specie for later depth calculation
case "${SPECIES}" in 
	"Enterococcus_faecalis")
		GENOME_SIZE=2990000
		;;
	"Enterococcus_faecium")
		GENOME_SIZE=3050000
		;;
	*)
		echo "Error: Unknown species"; Help ; exit 1
esac

#########################################################################

# Quality of ONT data as decision parameter for pipeline selection

## Low/intermediate ONT quality and coverage:
	### Unicycler short-read-first hybrid assembly

## High ONT quality, coverage and suficient read N50:
	### Autocycler long-read-first concensus assembly
	### Polypolish Illumina assembly

	# Conditions of high quality:
		# ONT depth >= 20x
		# ONT mean Q >=X
		# ONT N50 >=X

#########################################################################

#########################################################################
# ILLUMINA AND NANOPORE QUALITY CONTROL OF RAW READS
#########################################################################
echo ""
echo "----------------------------------------------"
echo "1. Quality control of Illumina raw reads"
echo "----------------------------------------------"
echo "------ 1.1. FastQC analysis ------"

mkdir -p "${WD}/raw_fastqc"
samples_analysed=0

R1_files=("${WD}"/*_R1_001.fastq.gz)
if [[ ${#R1_files[@]} -eq 0 ]]; then
	echo "ERROR: No *_R1_001.fastq.gz files found in ${WD}"
	exit 1
fi

conda activate fastqc

for R1 in "${WD}"/*_R1_001.fastq.gz; do
	R2=${R1/_R1_/_R2_}
	if [[ -f "${R2}" ]]; then
		name_sample=$(basename "${R1}" | sed 's/_R1_001.fastq.gz//')
		echo ""
		echo "Processing sample: ${name_sample}"
		fastqc -o "${WD}/raw_fastqc" "${R1}" "${R2}"
		((samples_analysed++))
	else
		echo "Missing R2 file"
	fi
done
conda deactivate

echo ""
echo "FASTQC completed. Number of samples analysed: ${samples_analysed}"

# Quality control with MultiQC
echo ""
echo "------ 1.2. MultiQC analysis ------"

if [[ -d "${WD}/raw_fastqc" && $(ls -A "${WD}/raw_fastqc") ]]; then
	conda activate multiqc
	multiqc ${WD}/raw_fastqc -o ${WD}/raw_fastqc/multiqc
	conda deactivate
	echo ""
	echo "MultiQC report completed"
else
	echo "MultiQC is not performed: raw_fastqc directory is empty"
fi

# Quality control for ONT raw reads with Nanoplot
echo ""
echo "----------------------------------------------------------------"
echo "2. Quality control of Oxford Nanopore raw reads with Nanoplot"
echo "----------------------------------------------------------------"

conda activate NanoPlot

mkdir -p "${WD}/Nanopore/NanoPlot/raw"

ONT_file=()

# Collecting only ONT files 
for file in "${WD}"/*.fastq.gz; do
	file_name=$(basename "${file}" )
	if [[ "${file_name}" == *_R1_001.fastq.gz || "${file_name}" == *_R2_001.fastq.gz ]]; then
		continue
	fi
	ONT_files+=("${file}")
done



for file in "${ONT_files[@]}"; do # Check name file Nanopore
	if [[ -f "${file}" ]]; then
		id=$(basename "${file}" .fastq.gz)
		echo "Running NanoPlot on ${id}"
        NanoPlot --fastq "${file}" --outdir "${WD}/Nanopore/NanoPlot/raw/${id}" --tsv_stats
	else
		echo "ERROR: No ONT input sequences detected"
		conda deactivate
		exit 1
	fi
done

conda deactivate

# Recovering raw ONT metrics from NanoStats (implemented in NanoPlot)

summary_stats_raw="${WD}/Nanopore/NanoPlot/summary_metrics_raw.tsv"
echo -e "Sample\tTotal_bases\tMean_quality\tN50\Depth" > "${summary_stats_raw}"

for stats_file in "${WD}/Nanopore/NanoPlot/raw"/*/NanoStats.txt; do
	if [[ -f "${stats_file}" ]]; then
		sample=$(basename "$(dirname "${stats_file}")")
		awk -v gsize="${GENOME_SIZE}" -v sample="${sample}" '
		BEGIN { FS=":" }
		/Total bases/ { gsub(/[ \t]/, "", $2); tb=$2 }
		/Mean read quality/ { gsub(/[ \t]/, "", $2); qual=$2 }
		/Read length N50/ { gsub(/[ \t]/, "", $2); n50=$2 }
		END {
			if (gsize > 0 && tb > 0) {
				depth = sprintf("%.2fx", tb / gsize)
			} else {
				depth = "N/A"
			}
			printf "Mostra: %s | Total Bases: %s | Quality: %s | N50: %s | Depth: %s\n", sample, tb, qual, n50, depth
			printf "%s\t%s\t%s\t%s\t%s\n", sample, tb, qual, n50, depth
			}' "${stats_file}" >> "${summary_stats_raw}"
	fi
done

#########################################################################
# FILTERING ILLUMINA AND ONT DATA
#########################################################################

# Define variables
trimmed_R1="${WD}/trimmomatic/${sample}_1_paired.fastq.gz"
trimmed_R2="${WD}/trimmomatic/${sample}_2_paired.fastq.gz"
filtered_ONT="${WD}/Nanopore/Filtlong/${sample}_filtered.fastq.gz"

echo ""
echo "---------------------------------------------------------"
echo "3. Filtering Illumina and ONT reads"
echo "---------------------------------------------------------"

# Trimming for Illumina raw reads with Trimmmomatic
echo ""
echo "------ 3.1.Trimmomatic for trimming raw reads ------"

conda activate trimmomatic

mkdir -p "${WD}/trimmomatic"
for R1 in "${WD}"/*_R1_001.fastq.gz; do
	R2=${R1/_R1_/_R2_}
	if [[ -f "${R2}" ]]; then
		id=$(basename "${R1}" _R1_001.fastq.gz)
		mod_id=$(echo "${id}" | sed -E 's/_S[0-9]+_L001//g')
		echo "Running trimmomatic on ${mod_id}"
		trimmomatic PE \
			"${R1}" "${R2}" \
			"${WD}/trimmomatic/${id}_1_paired.fastq.gz" \
                "${WD}/trimmomatic/${id}_1_unpaired.fastq.gz" \
                "${WD}/trimmomatic/${id}_2_paired.fastq.gz" \
                "${WD}/trimmomatic/${id}_2_unpaired.fastq.gz" \
			ILLUMINACLIP:/home/seqadmin/miniconda3/envs/trimmomatic/share/trimmomatic-0.40-0/adapters/TruSeq3-PE.fa:2:30:10 \
        	LEADING:3 \
			TRAILING:3 \
			SLIDINGWINDOW:5:25 \
			MINLEN:60
		echo "Trimming on ${mod_id} finished"
	fi
done

conda deactivate

# Quality control for Illumina trimmed reads

echo ""
echo "------ 3.2 Quality control after trimming with FastQC ------"

conda activate fastqc

mkdir -p "${WD}/trimmomatic/fastqc"
fastqc -o "${WD}/trimmomatic/fastqc" "${WD}/trimmomatic/"*_paired.fastq.gz

conda deactivate

# Quality control with MultiQC
echo ""
echo "------ 3.3 MultiQC analysis ------"

if [[ -d "${WD}/trimmomatic/fastqc" && $(ls -A "${WD}/trimmomatic/fastqc") ]]; then
	conda activate multiqc
	multiqc ${WD}/trimmomatic/fastqc -o ${WD}/trimmomatic/fastqc/multiqc
	conda deactivate
	echo ""
	echo "MultiQC report completed"
else
	echo "MultiQC is not performed: directory is empty"
fi

# Filtering for ONT raw reads
echo ""
echo "------ 3.4 Filtering Oxford Nanopore raw reads with Filtlong ------"

mkdir -p "${WD}/Nanopore/Filtlong"

conda activate Filtlong

for file in "${ONT_files[@]}"; do # Check name file Nanopore
	if [[ -f "${file}" ]]; then
		id=$(basename "${file}" .fastq.gz)
		echo "Running Filtlong on ${id}"
       	filtlong --min_length 1000 --keep_percent 90 --mean_q_weigth 9 "${file}" | gzip > "${WD}/Nanopore/Filtlong/${id}_filtered.fastq.gz"
	fi
done

conda deactivate

# Quality control of Nanopore filtered reads

echo ""
echo "------ 3.5 Quality control of Oxford NanoporeDespres es torna filtered reads with NanoPlot ------"

mkdir -p "${WD}/Nanopore/NanoPlot/filtered"

conda activate NanoPlot

for file in "${WD}/Nanopore/Filtlong"/*_filtered.fastq.gz; do
	if [[ -f "${file}" ]]; then
		id=$(basename "${file}" _filtered.fastq.gz)
       	NanoPlot --fastq "${file}" --outdir "${WD}/Nanopore/NanoPlot/filtered/${id}" --tsv_stats
	fi
done

conda deactivate

# Recovering raw ONT metrics from NanoStats (implemented in NanoPlot)
summary_stats_filtered="${WD}/Nanopore/NanoPlot/summary_metrics_filtered.tsv"
echo -e "Sample\tTotal_bases\tMean_quality\tN50\Depth" > "${summary_stats_filtered}"

for stats_file in "${WD}/Nanopore/NanoPlot/filtered"/*/NanoStats.txt; do
	if [[ -f "${stats_file}" ]]; then
	sample=$(basename "$(dirname "${stats_file}")")
	awk -v gsize="${GENOME_SIZE}" -v sample="${sample}" '
	BEGIN { FS=":" }
	/Total bases/ { gsub(/[ \t]/, "", $2); tb=$2 }
   	/Mean read quality/ { gsub(/[ \t]/, "", $2); qual=$2 }
   	/Read length N50/ { gsub(/[ \t]/, "", $2); n50=$2 }
   	END {
       	if (gsize > 0 && tb > 0) {
           	depth = sprintf("%.2fx", tb / gsize)
       	} else {
           	depth = "N/A"
       	}
       	printf "Mostra: %s | Total Bases: %s | Quality: %s | N50: %s | Depth: %s\n", sample, tb, qual, n50, depth
       	printf "%s\t%s\t%s\t%s\t%s\n", sample, tb, qual, n50, depth
       	}' "${stats_file}" >> "${summary_stats_filtered}"
	fi
done


#########################################################################
# ASSEMBLY STRATEGY DECISION
#########################################################################
echo ""
echo "----------------------------------------------------------------"
echo "4. Aseembly stratgey decision"
echo "----------------------------------------------------------------"

# Stablishing variables of ONT high quality
ONT_MIN_DEPTH=20
ONT_MIN_MEAN_Q=
ONT_MIN_N50=

mkdir -p "${WD}/final_assemblies"

if [[ ! "${summary_stats_raw}" ]]; then
	echo "Error: No summary for raw ONT data has been found"
	exit 1
fi

# Read the summary file to extract quality values
tail -n +2 "${summary_stats_raw}" | while IFS=$'\t' read -r sample total_bases mean_quality n50 depth; do
	num_depth=$(echo "${depth}" | sed 's/x//g')
	long_read_first="NO" # Predefined value
	if [[ -n "${num_depth}" ]]; then
		if (( $(echo "${num_depth} >= ${ONT_MIN_DEPTH}" | bc -l) )) \
			&& (( $(echo "${mean_quality} >= ${ONT_MIN_MEAN_Q}" | bc -l) )) \
			&& (( n50 >= ${ONT_MIN_N50} )); then
			long_read_first="YES"
		fi
	fi

#########################################################################
# HYBRID ASSEMBLY
#########################################################################

echo ""
echo "----------------------------------------------------------------"
echo "5. Hybrid assembly"
echo "----------------------------------------------------------------"

if [[ "${long_read_first}" == "NO" ]]; then

	# Unicycler assembles the Illumina reads first and then completes it with Nanopore reads
	echo ""
	echo "---------------------------------------------------------"
	echo "Unicycler for hybrid assembly (short-reads-first)"
	echo "---------------------------------------------------------"
	echo ""
	echo "------ Running Unicycler hybrid assembler ------"
	mkdir -p "${WD}/unicycler/${sample}"

	conda activate unicycler

	unicycler -1 "${trimmed_R1}" \
    	-2 "${trimmed_R2}" \
    	-l "${filtered_ONT}"\
    	-o "${WD}/unicycler/${sample}" \
    	--mode normal

	conda deactivate

	if [[ -f "${WD}/unicycler/${sample}/assembly.fasta" ]]; then
		echo "Hybrid assembly generated"
		cp "${WD}/unicycler/${sample}/assembly.fasta" "${WD}/final_assembly/${sample}_final.fasta"
	else
		echo "Error: Unicycler could not generate the hybrid assembly"
	fi
else

	# Autocycler long-read-first assembly
	echo ""
	echo " Hybrid assembly generated with AUtocycler + Polypolish"
	mkdir -p "${WD}/autocycler/${sample}"
	mkdir -p "${WD}/autocycler/polypolish/${sample}"

	conda activate autocycler

	mkdir -p "${WD}/autocycler/${sample}/subsampled_reads"

	autocycler subsample \
		--count 4 \
		--reads "${filtered_ONT}" \
		--out_dir "${WD}/autocycler/${sample}/subsampled_reads" \
		--genome_size "${GENOME_SIZE}" \
		--min_read_depth 50

	# Combining the assemblies into a consensus assembly
	mkdir -p "${WD}/autocycler/${sample}/assemblies"
	for assembler in flye miniasm metamdbg; do
		for read_file in "${WD}/autocycler/${sample}/subsampled_reads"/*.fastq.gz; do
			if [[ -f "${read_file}" ]]; then
				read_id=$(basename "${read_file}" .fastq)
				autocycler helper "${assembler}" \
					--reads "${read_file}"\
					--out_prefix "${WD}/autocycler/${sample}/assemblies/${assembler}_${read_id}" \
					--genome_size "${GENOME_SIZE}"
			fi
		done
	done	

	mkdir -p "${WD}/autocycler/${sample}/autocycler_compress_cluster_output"
	autocycler compress -i "${WD}/autocycler/${sample}/assembllies" -a "${WD}/autocycler/${sample}/autocycler_compress_cluster_output"
	autocycler cluster -a "${WD}/autocycler/${sample}/autocycler_compress_cluster_output"

	for cluster in "${WD}/autocycler/${sample}/autocycler_compress_cluster_output/clustering/qc_pass"/cluster_*; do
		autocycler trim -c "${cluster}"
		autocycler resolve -c "${cluster}"
	done

	autocycler combine \
		-a "${WD}/autocycler/${sample}/autocycler_compress_cluster_output" \
		-i "${WD}/autocycler/${sample}/autocycler_compress_cluster_output/clustering/qc_pass"/cluster_*/5_final.gfa
	
	# Convert GFA file to FASTA file
	gfa_output="${WD}/autocycler/${sample}/autocycler_compress_cluster_output/clustering/qc_pass"/cluster_*/5_final.gfa
	fasta_output="${WD}/autocycler/${sample}/autocycler_assembly_${sample}.fasta"
	autocycler gfa2fasta \
		-i "${gfa_output}" \
		-o "${fasta_output}"
	
	conda deactivate

	# Medaka ONT assembly polishing with ONT reads
	mkdir -p "${WD}/medaka/${sample}"

	conda activate medaka

	medaka_consensus \
		-i "${filtered_ONT}" \
		-d "${fasta_output}" \
		-o "${WD}/medaka/${sample}"

	conda deactivate

	medaka_assembly="{WD}/medaka/${sample}/consensus.fasta"

	# Illumina polish with PolyPolish
	
	mkdir -p "{WD}/polypolsih/${sample}"

	conda activate polypolish

	## Index assembly with BWA
	bwa index "${medaka_assembly}"
	# Map Illumina reads
	bwa mem -a "${medaka_assembly}" "${trimmed_R1}" > "${WD}/polypolish/${sample}/alignments_R1.sam"
	bwa mem -a "${medaka_assembly}" "${trimmed_R2}" > "${WD}/polypolish/${sample}/alignments_R2.sam"

	# Filter alignments
	polypolish filter \
		--in1 "${WD}/polypolish/${sample}/alignments_R1.sam" \
		--in2 "${WD}/polypolish/${sample}/alignments_R2.sam" \
		--out1 "${WD}/polypolish/${sample}/filtered_R1.sam" \
		--out2 "${WD}/polypolish/${sample}/filtered_R2.sam"
	# Final polish step
	mkdir -p "${WD}/final_assembly"
	polypolish polish \
		"${autocycler_output}" \
		"${WD}/polypolish/${sample}/filtered_R1.sam" \
		"${WD}/polypolish/${sample}/filtered_R2.sam" > "${WD}/final_assembly/${sample}_final.fasta"
	
	conda deactivate

fi

#########################################################################
# ASSESING ASSEMBLY QUALITY
#########################################################################
# Quality control assembly with BUSCO
echo ""
echo "------ 5.1 Quality control of hybrid assembly with BUSCO ------"

mkdir -p "${WD}/busco"

conda activate busco

for fa in "${WD}/final_assemblies/"*.fasta; do
	if [[ -f "${fa}" ]]; then
    	id=$(basename "${fa}" _final.fasta)
        busco -i "${fa}" -m genome -l enterococcus_odb12 -o "${id}" -c 4
    fi
done
mkdir -p "${WD}/busco/json"
cp "${WD}"/busco//*/short_summary.*.json "${WD}"/busco/json

# Generating BUSCO plots
echo "Generating BUSCO plot results"
busco --plot "${WD}/busco/json"

conda deactivate

echo ""
echo "BUSCO analysis for evaluation of assembly completed"

# Quality control assembly with QUAST

echo ""
echo "------ 5.2 Quality control of hybrid assembly with QUAST ------"

mkdir -p "${WD}/quast"

conda activate quast

quast_input=("${WD}/final_assemblies/"*.fasta)
if [[ ${#quast_input[@]} -gt 0 ]]; then
	quast "${WD}/final_assemblies/"*.fasta -o "${WD}/quast"
	if [[ -f "${WD}/quast/transposed_report.tsv" ]]; then
		cp "${WD}/quast/transposed_report.tsv" "${WD}/quast/${logstamp}_summary.tsv"
	fi
else
	echo "Skipping QUAST, no hybrid assemblies found"
fi

conda deactivate

#########################################################################
# READ MAPPING TO FINAL ASSEMBLY
#########################################################################
echo ""
echo "----------------------------------------------------------------"
echo "6. Read mapping to final assembly"
echo "----------------------------------------------------------------"

mkdir -p "${WD}/mapping"
mkdir -p "${WD}/mapping/${sample}"
for fa in "${WD}/final_assemblies/"*_final.fasta; do
	if [[ -f "${fa}" ]]; then
		sample=(basename "${fa}" _final.fasta)
		final_size=$(awk '/^>/ {next} {sum += length($0)} END {print sum+0}' "${fa}")
        final_contigs=$(grep -c '^>' "${fa}")
		# Nanopore
		if [[ -f "${filtered_ONT}" ]]; then
			ONT_BAM="${WD}/mapping/${sample}/${sample}_ONT.bam"

			conda activate minimap2

			minimap2 -ax map-ont "${fa}" "${filtered_ONT}" | samtools sort -@ 4 -o "${ont_bam}"
			samtools index "${ont_bam}"
			samtools flagstat "${ont_bam}" > "${WD}/mapping/${sample}/${sample}_ONT_flagstat.txt"

			conda deactivate
		fi
		
		# Illumina
		if [[ -f "${trimmed_R1}" && "${trimmed_R2}" ]]; then
			illumina_bam="${WD}/mapping/${sample}/${sample}_Illumina.bam"

			conda activate bwa

			conda activate samtools

			bwa index "${fa}"
			bwa mem "${fa}" "${trimmed_R1}" "${trimmed_R2}" | samtools sort -@ 4 -o "${illumina_bam}"
			samtools index "${illumina_bam}"
			samtools flagstat "${illumina_bam}" > "${WD}/mapping/${sample}/${sample}_Illumina_flagstat.txt"

			conda deactivate

			conda deactivate
		fi
	fi
done

echo ""
echo "Hybrid assembly pipeline completed successfully"

############################################
# Plasmid detection
############################################

# Plassembler (recover small plasmids)
echo ""
echo "Plasmid detection: Plassembler"

mkdir -p "${WD}/autocycler/plassembler/${sample}"

conda activate plassembler

# Crear un database de plassembler
if [[ -d "/miniconda/envs/plassembler/database" ]]; then
	echo "Database directory is created"
else
	mkdir -p "/miniconda/envs/plassembler/database"
	conda activate plassembler
	plassembler download -d "/miniconda/envs/plassembler/database"
	conda deactivate

plassembler run \
	-d "${plassembler_database}" \
	-l "${filtered_ONT}" \
	-o "${WD}/autocycler/plassembler/${sample}" \
	-1 "${trimmed_R1}" \
	-2 "${trimmed_R2}" \
	-c "${GENOME_SIZE}"

conda deactivate

plassembler_result="${WD}/autocycler/plassembler/${sample}/plassembler.fasta"
if [[ -f "${plassembler_result}" ]]; then
	cp "${plassembler_result}" "${WD}/autocycler/${sample}/assemblies/plassembler.fasta"
else
	echo "Error: Plassembler results not found"	

#################################
# End of draft hybrid assembly
#################################
