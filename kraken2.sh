#!/bin/bash
#SBATCH --job-name=kraken2_full


set -euo pipefail

############################################
# 1. LOAD SOFTWARE
############################################

conda activate kraken2_env

############################################
# 2. USER SETTINGS
############################################

# Directory containing your paired-end FASTQ files
READDIR="data"

# Main output directory
OUTDIR="kraken2_results"

# Existing built Kraken2 database directory
DBDIR="kraken2_db"

# Number of CPU threads
THREADS="${SLURM_CPUS_PER_TASK:-16}"

############################################
# 3. MAKE DIRECTORIES
############################################

mkdir -p "${OUTDIR}"

REPORT_DIR="${OUTDIR}/reports"
KRAKEN_DIR="${OUTDIR}/kraken_output"
CLASSIFIED_DIR="${OUTDIR}/classified_reads"
UNCLASSIFIED_DIR="${OUTDIR}/unclassified_reads"
LOG_DIR="${OUTDIR}/sample_logs"

mkdir -p "${REPORT_DIR}" "${KRAKEN_DIR}" "${CLASSIFIED_DIR}" "${UNCLASSIFIED_DIR}" "${LOG_DIR}"

echo "========================================"
echo "Kraken2 pipeline started: $(date)"
echo "READDIR = ${READDIR}"
echo "OUTDIR  = ${OUTDIR}"
echo "DBDIR   = ${DBDIR}"
echo "THREADS = ${THREADS}"
echo "========================================"

############################################
# 4. CHECK DATABASE
############################################

if [[ ! -f "${DBDIR}/hash.k2d" || ! -f "${DBDIR}/opts.k2d" || ! -f "${DBDIR}/taxo.k2d" ]]; then
    echo "ERROR: ${DBDIR} does not appear to be a valid built Kraken2 database."
    echo "Missing one or more of: hash.k2d, opts.k2d, taxo.k2d"
    exit 1
fi

############################################
# 5. CLASSIFY ALL PAIRED-END SAMPLES
############################################

echo "Starting sample classification: $(date)"

FOUND_ANY="no"

for R1 in "${READDIR}"/*_R1*.fastq.gz; do
    if [[ ! -e "${R1}" ]]; then
        echo "No R1 FASTQ files found in ${READDIR}"
        break
    fi

    FOUND_ANY="yes"

    R2="${R1/_R1/_R2}"

    if [[ ! -f "${R2}" ]]; then
        echo "WARNING: Missing R2 pair for:"
        echo "  ${R1}"
        echo "Expected:"
        echo "  ${R2}"
        echo "Skipping this sample."
        continue
    fi

    BASENAME=$(basename "${R1}")
    SAMPLE=${BASENAME%%_R1*}

    echo "----------------------------------------"
    echo "Processing sample: ${SAMPLE}"
    echo "R1: ${R1}"
    echo "R2: ${R2}"
    echo "Start: $(date)"

    kraken2 \
        --db "${DBDIR}" \
        --threads "${THREADS}" \
        --paired \
        --gzip-compressed \
        --use-names \
        --report "${REPORT_DIR}/${SAMPLE}.report" \
        --output "${KRAKEN_DIR}/${SAMPLE}.kraken" \
        --classified-out "${CLASSIFIED_DIR}/${SAMPLE}_classified_#.fastq" \
        --unclassified-out "${UNCLASSIFIED_DIR}/${SAMPLE}_unclassified_#.fastq" \
        "${R1}" "${R2}" \
        > "${LOG_DIR}/${SAMPLE}.stdout.log" \
        2> "${LOG_DIR}/${SAMPLE}.stderr.log"

    echo "Finished sample: ${SAMPLE}"
    echo "End: $(date)"
done

if [[ "${FOUND_ANY}" == "no" ]]; then
    echo "ERROR: No input files matching *_R1*.fastq.gz were found in ${READDIR}"
    exit 1
fi

############################################
# 6. SUMMARIZE REPORT FILES
############################################

SUMMARY_FILE="${OUTDIR}/kraken2_summary_top_hits.tsv"

echo -e "Sample\tPercent\tReads_clade\tReads_direct\tRank\tTaxID\tName" > "${SUMMARY_FILE}"

for REPORT in "${REPORT_DIR}"/*.report; do
    [[ -e "${REPORT}" ]] || continue

    SAMPLE=$(basename "${REPORT}" .report)

    awk -v sample="${SAMPLE}" '
        BEGIN{OFS="\t"}
        $4 != "0" && $6 !~ /unclassified/ {
            print sample, $1, $2, $3, $4, $5, $6
        }
    ' "${REPORT}" | head -n 10 >> "${SUMMARY_FILE}"
done

############################################
# 7. FINISH
############################################

echo "========================================"
echo "Kraken2 pipeline finished: $(date)"
echo "Main outputs:"
echo "  Reports: ${REPORT_DIR}"
echo "  Per-read output: ${KRAKEN_DIR}"
echo "  Classified reads: ${CLASSIFIED_DIR}"
echo "  Unclassified reads: ${UNCLASSIFIED_DIR}"
echo "  Logs: ${LOG_DIR}"
echo "  Summary: ${SUMMARY_FILE}"
echo "========================================"
