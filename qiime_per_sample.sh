#!/bin/bash
#SBATCH --job-name=qiime_classifier_iqtree

set -euo pipefail

###############################################################################
# USAGE
#   sbatch qiime2_per_sample.sh EF1
#   sbatch qiime2_per_sample.sh ITS1
#   sbatch qiime2_per_sample.sh ITS2
###############################################################################
MARKER="${1:-ITS2}"

###############################################################################
# INPUT FILES
###############################################################################

MANIFEST="manifest.tsv"
METADATA="metadata.tsv"

###############################################################################
# REFERENCE DATABASES
###############################################################################

FUSARIUM_RAW_FASTA="fusarioidIDversion10022026.fas"
UNITE_FASTA="UNITE_QIIME/sh_refs_qiime_ver10_dynamic_19.02.2025.fasta"
UNITE_TAXONOMY="UNITE_QIIME/sh_taxonomy_qiime_ver10_dynamic_19.02.2025.txt"

###############################################################################
# CLASSIFIER OUTPUTS
###############################################################################

CLASSIFIER_DIR="classifiers"
mkdir -p "${CLASSIFIER_DIR}"

EF1_CLASSIFIER="${CLASSIFIER_DIR}/fusarioidid_ef1_classifier.qza"
ITS1_CLASSIFIER="${CLASSIFIER_DIR}/unite_its1_classifier.qza"
ITS2_CLASSIFIER="${CLASSIFIER_DIR}/unite_its2_classifier.qza"

THREADS="${SLURM_CPUS_PER_TASK:-12}"

###############################################################################
# SOFTWARE
###############################################################################

source /home/tntech.edu/ssalimi42/miniconda3/etc/profile.d/conda.sh
conda activate /home/software/qiime2-2023.5

spack load clustal-omega
spack load trimal
spack load iq-tree@2.3.6

###############################################################################
# MARKER-SPECIFIC SETTINGS
###############################################################################

case "${MARKER}" in
  EF1)
    F_PRIMER="CCGGTCACTTGATCTACCAG"
    R_PRIMER="ATGACGGTGACATAGTAGCG"
    TRIM_LEFT_F=0
    TRIM_LEFT_R=0
    TRUNC_LEN_F=220
    TRUNC_LEN_R=220
    CLASSIFIER="${EF1_CLASSIFIER}"
    ;;

  ITS1)
    F_PRIMER="ACCTGCGGARGGATCA"
    R_PRIMER="GAGATCCRTTGYTRAAAGTT"
    TRIM_LEFT_F=0
    TRIM_LEFT_R=0
    TRUNC_LEN_F=200
    TRUNC_LEN_R=200
    CLASSIFIER="${ITS1_CLASSIFIER}"
    ;;

  ITS2)
    F_PRIMER="AACTTTYRRCAAYGGATCWCT"
    R_PRIMER="AGCCTCCGCTTATTGATATGCTTAART"
    TRIM_LEFT_F=0
    TRIM_LEFT_R=0
    TRUNC_LEN_F=200
    TRUNC_LEN_R=200
    CLASSIFIER="${ITS2_CLASSIFIER}"
    ;;

  *)
    echo "ERROR: MARKER must be one of: EF1, ITS1, ITS2"
    exit 1
    ;;
esac

###############################################################################
# GENERAL FILTERING SETTINGS
###############################################################################

MIN_FEATURE_FREQ=5
MIN_SEQ_LEN=100

###############################################################################
# CHECK INPUTS
###############################################################################

for f in "${MANIFEST}" "${METADATA}"; do
  if [[ ! -f "${f}" ]]; then
    echo "ERROR: Required file not found: ${f}"
    exit 1
  fi
done

###############################################################################
# BUILD EF1 CLASSIFIER FROM FusarioidID FILE
###############################################################################

build_ef1_classifier() {
  if [[ -f "${EF1_CLASSIFIER}" ]]; then
    echo "[$(date)] EF1 classifier already exists: ${EF1_CLASSIFIER}"
    return
  fi

  echo "[$(date)] Building EF1 FusarioidID classifier from ${FUSARIUM_RAW_FASTA} ..."

  if [[ ! -f "${FUSARIUM_RAW_FASTA}" ]]; then
    echo "ERROR: FusarioidID FASTA not found: ${FUSARIUM_RAW_FASTA}"
    exit 1
  fi

  local workdir="${CLASSIFIER_DIR}/ef1_build"
  mkdir -p "${workdir}"

  local ef1_clean_fasta="${workdir}/fusarioid_tef1_only.fasta"
  local ef1_ref_qza="${workdir}/ef1_ref_seqs.qza"
  local ef1_tax_tsv="${workdir}/ef1_taxonomy.tsv"
  local ef1_tax_qza="${workdir}/ef1_ref_taxonomy.qza"

  python <<'PY'
import re

infile = "fusarioidIDversion10022026.fas"
fasta_out = "classifiers/ef1_build/fusarioid_tef1_only.fasta"
tax_out = "classifiers/ef1_build/ef1_taxonomy.tsv"

def strip_html(text):
    return re.sub(r"<[^>]+>", "", text)

def clean_seq(seq_lines):
    seq = "".join(seq_lines)
    seq = seq.replace("\r", "").replace("\n", "")
    seq = seq.upper()
    seq = re.sub(r"[^ACGTRYKMSWBDHVN]", "", seq)
    return seq

count = 0
header = None
seq_lines = []

with open(infile, "r", encoding="utf-8", errors="ignore") as fin, \
     open(fasta_out, "w", encoding="utf-8") as fout_fa, \
     open(tax_out, "w", encoding="utf-8") as fout_tax:

    def flush_record(h, s_lines):
        global count
        if h is None:
            return

        h2 = strip_html(h.strip())
        if "|tef1" not in h2.lower():
            return

        seq = clean_seq(s_lines)
        if not seq:
            return

        parts = h2[1:].split("|") if h2.startswith(">") else h2.split("|")
        if len(parts) < 5:
            return

        taxon_field = parts[3].strip()
        primary_name = taxon_field.split(",")[0].strip()
        if not primary_name:
            return

        primary_name_clean = re.sub(r"\s+", "_", primary_name)
        genus = primary_name.split()[0]
        genus_clean = re.sub(r"\s+", "_", genus)

        count += 1
        feature_id = f"EF1_{count}"

        taxonomy = (
            "k__Fungi;"
            "p__Ascomycota;"
            "c__Sordariomycetes;"
            "o__Hypocreales;"
            "f__Nectriaceae;"
            f"g__{genus_clean};"
            f"s__{primary_name_clean}"
        )

        fout_fa.write(f">{feature_id}\n{seq}\n")
        fout_tax.write(f"{feature_id}\t{taxonomy}\n")

    for line in fin:
        if line.startswith(">"):
            flush_record(header, seq_lines)
            header = line.rstrip("\r\n")
            seq_lines = []
        else:
            seq_lines.append(line)

    flush_record(header, seq_lines)
PY

  qiime tools import \
    --type 'FeatureData[Sequence]' \
    --input-path "${ef1_clean_fasta}" \
    --output-path "${ef1_ref_qza}"

  qiime tools import \
    --type 'FeatureData[Taxonomy]' \
    --input-format HeaderlessTSVTaxonomyFormat \
    --input-path "${ef1_tax_tsv}" \
    --output-path "${ef1_tax_qza}"

  qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads "${ef1_ref_qza}" \
    --i-reference-taxonomy "${ef1_tax_qza}" \
    --o-classifier "${EF1_CLASSIFIER}"
}

###############################################################################
# BUILD ITS1 / ITS2 CLASSIFIERS FROM UNITE
###############################################################################

build_its_classifier() {
  local marker_name="$1"
  local f_primer="$2"
  local r_primer="$3"
  local out_classifier="$4"

  if [[ -f "${out_classifier}" ]]; then
    echo "[$(date)] ${marker_name} classifier already exists: ${out_classifier}"
    return
  fi

  local workdir="${CLASSIFIER_DIR}/${marker_name}_build"
  mkdir -p "${workdir}"

  local unite_ref_qza="${workdir}/unite_ref_seqs.qza"
  local unite_tax_qza="${workdir}/unite_ref_taxonomy.qza"
  local unite_trimmed_qza="${workdir}/unite_ref_trimmed.qza"

  qiime tools import \
    --type 'FeatureData[Sequence]' \
    --input-path "${UNITE_FASTA}" \
    --output-path "${unite_ref_qza}"

  qiime tools import \
    --type 'FeatureData[Taxonomy]' \
    --input-format HeaderlessTSVTaxonomyFormat \
    --input-path "${UNITE_TAXONOMY}" \
    --output-path "${unite_tax_qza}"

  qiime feature-classifier extract-reads \
    --i-sequences "${unite_ref_qza}" \
    --p-f-primer "${f_primer}" \
    --p-r-primer "${r_primer}" \
    --p-trunc-len 250 \
    --o-reads "${unite_trimmed_qza}"

  qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads "${unite_trimmed_qza}" \
    --i-reference-taxonomy "${unite_tax_qza}" \
    --o-classifier "${out_classifier}"
}

###############################################################################
# BUILD NEEDED CLASSIFIER
###############################################################################

case "${MARKER}" in
  EF1)
    build_ef1_classifier
    ;;
  ITS1)
    build_its_classifier "ITS1" "ACCTGCGGARGGATCA" "GAGATCCRTTGYTRAAAGTT" "${ITS1_CLASSIFIER}"
    ;;
  ITS2)
    build_its_classifier "ITS2" "AACTTTYRRCAAYGGATCWCT" "AGCCTCCGCTTATTGATATGCTTAART" "${ITS2_CLASSIFIER}"
    ;;
esac

###############################################################################
# PER-SAMPLE ROOT
###############################################################################

PER_SAMPLE_DIR="per_sample_${MARKER}"
mkdir -p "${PER_SAMPLE_DIR}"

###############################################################################
# LOOP OVER SAMPLES
###############################################################################

awk 'BEGIN{FS="\t"} NR>1 && NF>0 {print $1}' "${METADATA}" | while read -r SAMPLE_ID; do
  echo "[$(date)] ========================================================"
  echo "[$(date)] Processing sample: ${SAMPLE_ID}"
  echo "[$(date)] ========================================================"

  SAMPLE_DIR="${PER_SAMPLE_DIR}/${SAMPLE_ID}"
  mkdir -p "${SAMPLE_DIR}"

  SAMPLE_MANIFEST="${SAMPLE_DIR}/${SAMPLE_ID}_manifest.tsv"
  SAMPLE_METADATA="${SAMPLE_DIR}/${SAMPLE_ID}_metadata.tsv"

  # Build one-sample manifest
  awk -v sid="${SAMPLE_ID}" 'BEGIN{FS=OFS="\t"} NR==1 || $1==sid {print}' "${MANIFEST}" > "${SAMPLE_MANIFEST}"

  # Build one-sample metadata
  grep -P "^#|^${SAMPLE_ID}\t" "${METADATA}" > "${SAMPLE_METADATA}"

  # Skip if sample not found in manifest
  if [[ $(wc -l < "${SAMPLE_MANIFEST}") -le 1 ]]; then
    echo "[$(date)] WARNING: Sample ${SAMPLE_ID} not found in manifest. Skipping."
    continue
  fi

  DEMUX_QZA="${SAMPLE_DIR}/${SAMPLE_ID}-${MARKER}-paired-end-demux.qza"
  DEMUX_QZV="${SAMPLE_DIR}/${SAMPLE_ID}-${MARKER}-paired-end-demux.qzv"

  TRIMMED_QZA="${SAMPLE_DIR}/${SAMPLE_ID}-${MARKER}-trimmed-paired-end-demux.qza"
  TRIMMED_QZV="${SAMPLE_DIR}/${SAMPLE_ID}-${MARKER}-trimmed-paired-end-demux.qzv"

  TABLE_QZA="${SAMPLE_DIR}/${SAMPLE_ID}-${MARKER}-trimmed-table.qza"
  REP_QZA="${SAMPLE_DIR}/${SAMPLE_ID}-${MARKER}-trimmed-rep-seqs.qza"
  STATS_QZA="${SAMPLE_DIR}/${SAMPLE_ID}-${MARKER}-trimmed-denoising-stats.qza"

  TABLE_QZV="${SAMPLE_DIR}/${SAMPLE_ID}-${MARKER}-trimmed-table.qzv"
  REP_QZV="${SAMPLE_DIR}/${SAMPLE_ID}-${MARKER}-trimmed-rep-seqs.qzv"
  STATS_QZV="${SAMPLE_DIR}/${SAMPLE_ID}-${MARKER}-trimmed-denoising-stats.qzv"

  RETAINED_TABLE_QZA="${SAMPLE_DIR}/${SAMPLE_ID}-${MARKER}-retained-hits-table.qza"
  RETAINED_TABLE_QZV="${SAMPLE_DIR}/${SAMPLE_ID}-${MARKER}-retained-hits-table.qzv"
  RETAINED_REP_QZA="${SAMPLE_DIR}/${SAMPLE_ID}-${MARKER}-retained-hits-rep-seqs.qza"
  RETAINED_REP_QZV="${SAMPLE_DIR}/${SAMPLE_ID}-${MARKER}-retained-hits-rep-seqs.qzv"

  LENFILTER_REP_QZA="${SAMPLE_DIR}/${SAMPLE_ID}-${MARKER}-retained-more-${MIN_SEQ_LEN}-rep-seqs.qza"
  LENFILTER_REP_QZV="${SAMPLE_DIR}/${SAMPLE_ID}-${MARKER}-retained-more-${MIN_SEQ_LEN}-rep-seqs.qzv"
  LENFILTER_TABLE_QZA="${SAMPLE_DIR}/${SAMPLE_ID}-${MARKER}-retained-more-${MIN_SEQ_LEN}-table.qza"
  LENFILTER_TABLE_QZV="${SAMPLE_DIR}/${SAMPLE_ID}-${MARKER}-retained-more-${MIN_SEQ_LEN}-table.qzv"

  TAXONOMY_QZA="${SAMPLE_DIR}/${SAMPLE_ID}-${MARKER}-taxonomy.qza"
  TAXONOMY_QZV="${SAMPLE_DIR}/${SAMPLE_ID}-${MARKER}-taxonomy.qzv"
  TAXA_BARPLOT_QZV="${SAMPLE_DIR}/${SAMPLE_ID}-${MARKER}-taxa-barplot.qzv"

  TREE_DIR="${SAMPLE_DIR}/iqtree"
  mkdir -p "${TREE_DIR}"
  TREE_FASTA="${TREE_DIR}/${SAMPLE_ID}-${MARKER}-sequences.fasta"
  ALIGNED_FASTA="${TREE_DIR}/${SAMPLE_ID}-${MARKER}_aligned.fa"
  TRIMMED_FASTA="${TREE_DIR}/${SAMPLE_ID}-${MARKER}_trimmed.fa"

  TMP_EXPORT_DIR="${SAMPLE_DIR}/tmp_retained_export"
  TMP_LEN_FASTA="${SAMPLE_DIR}/${SAMPLE_ID}-${MARKER}_tmp_lenfiltered.fasta"
  TMP_TREE_EXPORT="${TREE_DIR}/exported_seqs"

  #############################################################################
  # STEP 1. IMPORT READS
  #############################################################################

  echo "[$(date)] STEP 1: Importing reads for ${SAMPLE_ID}..."
  qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path "${SAMPLE_MANIFEST}" \
    --output-path "${DEMUX_QZA}" \
    --input-format PairedEndFastqManifestPhred33V2

  qiime demux summarize \
    --i-data "${DEMUX_QZA}" \
    --o-visualization "${DEMUX_QZV}"

  #############################################################################
  # STEP 2. TRIM PRIMERS
  #############################################################################

  echo "[$(date)] STEP 2: Trimming primers for ${SAMPLE_ID}..."
  qiime cutadapt trim-paired \
    --i-demultiplexed-sequences "${DEMUX_QZA}" \
    --p-front-f "${F_PRIMER}" \
    --p-front-r "${R_PRIMER}" \
    --p-error-rate 0.1 \
    --p-overlap 3 \
    --p-discard-untrimmed \
    --p-match-adapter-wildcards True \
    --p-match-read-wildcards False \
    --p-minimum-length "${MIN_SEQ_LEN}" \
    --o-trimmed-sequences "${TRIMMED_QZA}" \
    --verbose > "${SAMPLE_DIR}/${SAMPLE_ID}-${MARKER}-cutadapt-output.txt" 2>&1

  qiime demux summarize \
    --i-data "${TRIMMED_QZA}" \
    --o-visualization "${TRIMMED_QZV}"

  #############################################################################
  # STEP 3. DADA2
  #############################################################################

  echo "[$(date)] STEP 3: Running DADA2 for ${SAMPLE_ID}..."
  if ! qiime dada2 denoise-paired \
      --i-demultiplexed-seqs "${TRIMMED_QZA}" \
      --p-trim-left-f "${TRIM_LEFT_F}" \
      --p-trim-left-r "${TRIM_LEFT_R}" \
      --p-trunc-len-f "${TRUNC_LEN_F}" \
      --p-trunc-len-r "${TRUNC_LEN_R}" \
      --p-n-threads "${THREADS}" \
      --o-table "${TABLE_QZA}" \
      --o-representative-sequences "${REP_QZA}" \
      --o-denoising-stats "${STATS_QZA}" \
      --verbose > "${SAMPLE_DIR}/${SAMPLE_ID}-${MARKER}-dada2-output.txt" 2>&1; then
    echo "[$(date)] WARNING: DADA2 failed for ${SAMPLE_ID}. Skipping sample."
    continue
  fi

  qiime feature-table summarize \
    --i-table "${TABLE_QZA}" \
    --o-visualization "${TABLE_QZV}" \
    --m-sample-metadata-file "${SAMPLE_METADATA}"

  qiime feature-table tabulate-seqs \
    --i-data "${REP_QZA}" \
    --o-visualization "${REP_QZV}"

  qiime metadata tabulate \
    --m-input-file "${STATS_QZA}" \
    --o-visualization "${STATS_QZV}"

  #############################################################################
  # STEP 4. FILTER LOW-FREQUENCY FEATURES
  #############################################################################

  echo "[$(date)] STEP 4: Filtering low-frequency features for ${SAMPLE_ID}..."
  qiime feature-table filter-features \
    --i-table "${TABLE_QZA}" \
    --p-min-frequency "${MIN_FEATURE_FREQ}" \
    --o-filtered-table "${RETAINED_TABLE_QZA}"

  if ! qiime feature-table summarize \
      --i-table "${RETAINED_TABLE_QZA}" \
      --o-visualization "${RETAINED_TABLE_QZV}" \
      --m-sample-metadata-file "${SAMPLE_METADATA}" >/dev/null 2>&1; then
    echo "[$(date)] WARNING: No retained features for ${SAMPLE_ID} after frequency filtering. Skipping."
    continue
  fi

  qiime feature-table filter-seqs \
    --i-data "${REP_QZA}" \
    --i-table "${RETAINED_TABLE_QZA}" \
    --o-filtered-data "${RETAINED_REP_QZA}"

  qiime feature-table tabulate-seqs \
    --i-data "${RETAINED_REP_QZA}" \
    --o-visualization "${RETAINED_REP_QZV}"

  #############################################################################
  # STEP 5. FILTER SEQUENCES BY LENGTH
  #############################################################################

  echo "[$(date)] STEP 5: Filtering sequences by length for ${SAMPLE_ID}..."
  rm -rf "${TMP_EXPORT_DIR}"
  mkdir -p "${TMP_EXPORT_DIR}"

  qiime tools export \
    --input-path "${RETAINED_REP_QZA}" \
    --output-path "${TMP_EXPORT_DIR}"

  awk -v minlen="${MIN_SEQ_LEN}" '
  BEGIN { RS=">"; FS="\n" }
  NR > 1 {
    header=$1
    seq=""
    for (i=2; i<=NF; i++) seq=seq $i
    if (length(seq) >= minlen) {
      print ">" header
      print seq
    }
  }
  ' "${TMP_EXPORT_DIR}/dna-sequences.fasta" > "${TMP_LEN_FASTA}"

  if [[ ! -s "${TMP_LEN_FASTA}" ]]; then
    echo "[$(date)] WARNING: No sequences passed length filter for ${SAMPLE_ID}. Skipping."
    continue
  fi

  qiime tools import \
    --type 'FeatureData[Sequence]' \
    --input-path "${TMP_LEN_FASTA}" \
    --output-path "${LENFILTER_REP_QZA}"

  qiime feature-table tabulate-seqs \
    --i-data "${LENFILTER_REP_QZA}" \
    --o-visualization "${LENFILTER_REP_QZV}"

  qiime feature-table filter-features \
    --i-table "${RETAINED_TABLE_QZA}" \
    --m-metadata-file "${LENFILTER_REP_QZA}" \
    --o-filtered-table "${LENFILTER_TABLE_QZA}"

  if ! qiime feature-table summarize \
      --i-table "${LENFILTER_TABLE_QZA}" \
      --o-visualization "${LENFILTER_TABLE_QZV}" \
      --m-sample-metadata-file "${SAMPLE_METADATA}" >/dev/null 2>&1; then
    echo "[$(date)] WARNING: No features left after length filtering for ${SAMPLE_ID}. Skipping."
    continue
  fi

  rm -rf "${TMP_EXPORT_DIR}"
  rm -f "${TMP_LEN_FASTA}"

  #############################################################################
  # STEP 6. CLASSIFICATION
  #############################################################################

  echo "[$(date)] STEP 6: Classifying ${SAMPLE_ID}..."
  qiime feature-classifier classify-sklearn \
    --i-classifier "${CLASSIFIER}" \
    --i-reads "${LENFILTER_REP_QZA}" \
    --o-classification "${TAXONOMY_QZA}"

  qiime metadata tabulate \
    --m-input-file "${TAXONOMY_QZA}" \
    --o-visualization "${TAXONOMY_QZV}"

  qiime taxa barplot \
    --i-table "${LENFILTER_TABLE_QZA}" \
    --i-taxonomy "${TAXONOMY_QZA}" \
    --m-metadata-file "${SAMPLE_METADATA}" \
    --o-visualization "${TAXA_BARPLOT_QZV}"

  #############################################################################
  # EXPORT TABLE / SEQS
  #############################################################################

  qiime tools export \
    --input-path "${LENFILTER_TABLE_QZA}" \
    --output-path "${SAMPLE_DIR}/exported_table"

  if command -v biom >/dev/null 2>&1; then
    biom convert \
      -i "${SAMPLE_DIR}/exported_table/feature-table.biom" \
      -o "${SAMPLE_DIR}/${SAMPLE_ID}-${MARKER}-table.tsv" \
      --to-tsv
  fi

  qiime tools export \
    --input-path "${LENFILTER_REP_QZA}" \
    --output-path "${SAMPLE_DIR}/exported_rep_seqs"

  mv "${SAMPLE_DIR}/exported_rep_seqs/dna-sequences.fasta" \
     "${SAMPLE_DIR}/${SAMPLE_ID}-${MARKER}-rep-seqs.fasta"

  #############################################################################
  # STEP 7. TREE PER SAMPLE
  #############################################################################

  echo "[$(date)] STEP 7: Building tree for ${SAMPLE_ID}..."
  rm -rf "${TMP_TREE_EXPORT}"
  mkdir -p "${TMP_TREE_EXPORT}"

  qiime tools export \
    --input-path "${LENFILTER_REP_QZA}" \
    --output-path "${TMP_TREE_EXPORT}"

  mv "${TMP_TREE_EXPORT}/dna-sequences.fasta" "${TREE_FASTA}"
  rm -rf "${TMP_TREE_EXPORT}"

  # Count sequences; skip tree if fewer than 2
  NSEQ=$(grep -c '^>' "${TREE_FASTA}" || true)
  if [[ "${NSEQ}" -lt 2 ]]; then
    echo "[$(date)] WARNING: Sample ${SAMPLE_ID} has fewer than 2 sequences after filtering. Skipping tree."
    continue
  fi

  clustalo \
    -i "${TREE_FASTA}" \
    -o "${ALIGNED_FASTA}" \
    -t DNA \
    --threads="${THREADS}" \
    --force

  trimal \
    -in "${ALIGNED_FASTA}" \
    -out "${TRIMMED_FASTA}" \
    -automated1

  iqtree2 \
    -s "${TRIMMED_FASTA}" \
    --prefix "${TREE_DIR}/${SAMPLE_ID}-${MARKER}_tree" \
    -m MFP \
    -bb 1000 \
    -alrt 1000 \
    -nt AUTO \
    -redo

  echo "[$(date)] Finished sample: ${SAMPLE_ID}"
done

###############################################################################
# DONE
###############################################################################

echo
echo "[$(date)] DONE."
echo "Per-sample results are in: ${PER_SAMPLE_DIR}"
