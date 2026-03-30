#!/bin/bash
#SBATCH --job-name=qiime_classifier_iqtree


set -euo pipefail

###############################################################################
# USAGE
#   sbatch qiime2.sh EF1
#   sbatch qiime2.sh ITS1
#   sbatch qiime2.sh ITS2
###############################################################################

MARKER="${1:-ITS2}"

###############################################################################
# INPUT FILES
###############################################################################

MANIFEST="manifest_fixed.tsv"
METADATA="metadata.tsv"

###############################################################################
# REFERENCE DATABASES
###############################################################################

# New multi-locus FusarioidID file
FUSARIUM_RAW_FASTA="fusarioidIDversion10022026.fas"

# UNITE v10 dynamic fungal dataset
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
# BUILD EF1 CLASSIFIER FROM NEW FusarioidID FILE
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

  echo "[$(date)] Extracting TEF1 records only and creating safe IDs..."
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
        nonlocal_count = None
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

print(f"Wrote {count} TEF1 reference sequences.")
PY

  if [[ ! -s "${ef1_clean_fasta}" || ! -s "${ef1_tax_tsv}" ]]; then
    echo "ERROR: EF1 cleaned FASTA or taxonomy file is empty."
    exit 1
  fi

  echo "[$(date)] Preview of cleaned EF1 FASTA:"
  head -4 "${ef1_clean_fasta}"

  echo "[$(date)] Preview of EF1 taxonomy:"
  head -4 "${ef1_tax_tsv}"

  echo "[$(date)] Importing EF1 reference sequences..."
  qiime tools import \
    --type 'FeatureData[Sequence]' \
    --input-path "${ef1_clean_fasta}" \
    --output-path "${ef1_ref_qza}"

  echo "[$(date)] Importing EF1 taxonomy..."
  qiime tools import \
    --type 'FeatureData[Taxonomy]' \
    --input-format HeaderlessTSVTaxonomyFormat \
    --input-path "${ef1_tax_tsv}" \
    --output-path "${ef1_tax_qza}"

  echo "[$(date)] Training EF1 classifier directly from cleaned TEF1 references..."
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

  echo "[$(date)] Building ${marker_name} UNITE classifier..."

  if [[ ! -f "${UNITE_FASTA}" || ! -f "${UNITE_TAXONOMY}" ]]; then
    echo "ERROR: UNITE files not found."
    echo "  FASTA: ${UNITE_FASTA}"
    echo "  TAXONOMY: ${UNITE_TAXONOMY}"
    exit 1
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
# OUTPUT FILE NAMES
###############################################################################

DEMUX_QZA="${MARKER}-paired-end-demux.qza"
DEMUX_QZV="${MARKER}-paired-end-demux.qzv"

TRIMMED_QZA="${MARKER}-trimmed-paired-end-demux.qza"
TRIMMED_QZV="${MARKER}-trimmed-paired-end-demux.qzv"

TABLE_QZA="${MARKER}-trimmed-table.qza"
REP_QZA="${MARKER}-trimmed-rep-seqs.qza"
STATS_QZA="${MARKER}-trimmed-denoising-stats.qza"

TABLE_QZV="${MARKER}-trimmed-table.qzv"
REP_QZV="${MARKER}-trimmed-rep-seqs.qzv"
STATS_QZV="${MARKER}-trimmed-denoising-stats.qzv"

RETAINED_TABLE_QZA="${MARKER}-retained-hits-table.qza"
RETAINED_TABLE_QZV="${MARKER}-retained-hits-table.qzv"
RETAINED_REP_QZA="${MARKER}-retained-hits-rep-seqs.qza"
RETAINED_REP_QZV="${MARKER}-retained-hits-rep-seqs.qzv"

LENFILTER_REP_QZA="${MARKER}-retained-more-${MIN_SEQ_LEN}-rep-seqs.qza"
LENFILTER_REP_QZV="${MARKER}-retained-more-${MIN_SEQ_LEN}-rep-seqs.qzv"
LENFILTER_TABLE_QZA="${MARKER}-retained-more-${MIN_SEQ_LEN}-table.qza"
LENFILTER_TABLE_QZV="${MARKER}-retained-more-${MIN_SEQ_LEN}-table.qzv"

TAXONOMY_QZA="${MARKER}-taxonomy.qza"
TAXONOMY_QZV="${MARKER}-taxonomy.qzv"
TAXA_BARPLOT_QZV="${MARKER}-taxa-barplot.qzv"

TREE_DIR="iqtree_${MARKER}"
TREE_FASTA="${TREE_DIR}/${MARKER}-sequences.fasta"
ALIGNED_FASTA="${TREE_DIR}/${MARKER}_aligned.fa"
TRIMMED_FASTA="${TREE_DIR}/${MARKER}_trimmed.fa"

mkdir -p "${TREE_DIR}"

###############################################################################
# STEP 1. IMPORT READS
###############################################################################

echo "[$(date)] STEP 1: Importing paired-end reads..."

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path "${MANIFEST}" \
  --output-path "${DEMUX_QZA}" \
  --input-format PairedEndFastqManifestPhred33V2

qiime demux summarize \
  --i-data "${DEMUX_QZA}" \
  --o-visualization "${DEMUX_QZV}"

###############################################################################
# STEP 2. TRIM PRIMERS / ADAPTERS
###############################################################################

echo "[$(date)] STEP 2: Trimming primers with cutadapt..."

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
  --verbose > "${MARKER}-cutadapt-output.txt" 2>&1

qiime demux summarize \
  --i-data "${TRIMMED_QZA}" \
  --o-visualization "${TRIMMED_QZV}"

###############################################################################
# STEP 3. DADA2 DENOISING + MERGING
###############################################################################

echo "[$(date)] STEP 3: Running DADA2..."

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs "${TRIMMED_QZA}" \
  --p-trim-left-f "${TRIM_LEFT_F}" \
  --p-trim-left-r "${TRIM_LEFT_R}" \
  --p-trunc-len-f "${TRUNC_LEN_F}" \
  --p-trunc-len-r "${TRUNC_LEN_R}" \
  --p-n-threads "${THREADS}" \
  --o-table "${TABLE_QZA}" \
  --o-representative-sequences "${REP_QZA}" \
  --o-denoising-stats "${STATS_QZA}" \
  --verbose > "${MARKER}-dada2-output.txt" 2>&1

qiime feature-table summarize \
  --i-table "${TABLE_QZA}" \
  --o-visualization "${TABLE_QZV}" \
  --m-sample-metadata-file "${METADATA}"

qiime feature-table tabulate-seqs \
  --i-data "${REP_QZA}" \
  --o-visualization "${REP_QZV}"

qiime metadata tabulate \
  --m-input-file "${STATS_QZA}" \
  --o-visualization "${STATS_QZV}"

###############################################################################
# STEP 4. FILTER LOW-FREQUENCY FEATURES
###############################################################################

echo "[$(date)] STEP 4: Filtering ASVs with frequency < ${MIN_FEATURE_FREQ}..."

qiime feature-table filter-features \
  --i-table "${TABLE_QZA}" \
  --p-min-frequency "${MIN_FEATURE_FREQ}" \
  --o-filtered-table "${RETAINED_TABLE_QZA}"

qiime feature-table summarize \
  --i-table "${RETAINED_TABLE_QZA}" \
  --o-visualization "${RETAINED_TABLE_QZV}" \
  --m-sample-metadata-file "${METADATA}"

qiime feature-table filter-seqs \
  --i-data "${REP_QZA}" \
  --i-table "${RETAINED_TABLE_QZA}" \
  --o-filtered-data "${RETAINED_REP_QZA}"

qiime feature-table tabulate-seqs \
  --i-data "${RETAINED_REP_QZA}" \
  --o-visualization "${RETAINED_REP_QZV}"

###############################################################################
# STEP 5. FILTER REPRESENTATIVE SEQUENCES BY LENGTH
###############################################################################

echo "[$(date)] STEP 5: Keeping ASVs with length >= ${MIN_SEQ_LEN} bp..."

TMP_EXPORT_DIR="${TREE_DIR}/tmp_retained_export"
TMP_LEN_FASTA="${TREE_DIR}/tmp_lenfiltered.fasta"

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

qiime feature-table summarize \
  --i-table "${LENFILTER_TABLE_QZA}" \
  --o-visualization "${LENFILTER_TABLE_QZV}" \
  --m-sample-metadata-file "${METADATA}"

rm -rf "${TMP_EXPORT_DIR}"
rm -f "${TMP_LEN_FASTA}"

###############################################################################
# STEP 6. TAXONOMIC CLASSIFICATION
###############################################################################

echo "[$(date)] STEP 6: Classifying ${MARKER} ASVs..."

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
  --m-metadata-file "${METADATA}" \
  --o-visualization "${TAXA_BARPLOT_QZV}"

###############################################################################
# STEP 7. EXPORT SEQUENCES FOR TREE
###############################################################################

echo "[$(date)] STEP 7: Exporting representative sequences for tree..."

TMP_TREE_EXPORT="${TREE_DIR}/exported_seqs"
rm -rf "${TMP_TREE_EXPORT}"
mkdir -p "${TMP_TREE_EXPORT}"

qiime tools export \
  --input-path "${LENFILTER_REP_QZA}" \
  --output-path "${TMP_TREE_EXPORT}"

mv "${TMP_TREE_EXPORT}/dna-sequences.fasta" "${TREE_FASTA}"
rm -rf "${TMP_TREE_EXPORT}"

###############################################################################
# STEP 8. ALIGNMENT
###############################################################################

echo "[$(date)] STEP 8: Aligning sequences with Clustal Omega..."

clustalo \
  -i "${TREE_FASTA}" \
  -o "${ALIGNED_FASTA}" \
  -t DNA \
  --threads="${THREADS}" \
  --force

###############################################################################
# STEP 9. TRIM ALIGNMENT
###############################################################################

echo "[$(date)] STEP 9: Trimming alignment with trimAl..."

trimal \
  -in "${ALIGNED_FASTA}" \
  -out "${TRIMMED_FASTA}" \
  -automated1

###############################################################################
# STEP 10. BUILD TREE
###############################################################################

echo "[$(date)] STEP 10: Building phylogenetic tree with IQ-TREE2..."

iqtree2 \
  -s "${TRIMMED_FASTA}" \
  --prefix "${TREE_DIR}/${MARKER}_tree" \
  -m MFP \
  -bb 1000 \
  -alrt 1000 \
  -nt AUTO \
  -redo

###############################################################################
# DONE
###############################################################################

echo
echo "[$(date)] DONE."
echo
echo "===================== MAIN OUTPUTS ====================="
echo "Classifier used               : ${CLASSIFIER}"
echo "QIIME artifacts / visualizations:"
echo "  Imported demux              : ${DEMUX_QZA}"
echo "  Demux summary               : ${DEMUX_QZV}"
echo "  Trimmed demux               : ${TRIMMED_QZA}"
echo "  Trimmed demux summary       : ${TRIMMED_QZV}"
echo "  DADA2 table                 : ${TABLE_QZA}"
echo "  DADA2 rep seqs              : ${REP_QZA}"
echo "  DADA2 stats                 : ${STATS_QZA}"
echo "  DADA2 table summary         : ${TABLE_QZV}"
echo "  DADA2 rep seq summary       : ${REP_QZV}"
echo "  DADA2 stats summary         : ${STATS_QZV}"
echo "  Frequency-filtered table    : ${RETAINED_TABLE_QZA}"
echo "  Frequency-filtered rep seqs : ${RETAINED_REP_QZA}"
echo "  Length-filtered table       : ${LENFILTER_TABLE_QZA}"
echo "  Length-filtered rep seqs    : ${LENFILTER_REP_QZA}"
echo "  Taxonomy                    : ${TAXONOMY_QZA}"
echo "  Taxonomy table              : ${TAXONOMY_QZV}"
echo "  Taxa barplot                : ${TAXA_BARPLOT_QZV}"
echo
echo "Tree outputs:"
echo "  Tree FASTA                  : ${TREE_FASTA}"
echo "  Alignment                   : ${ALIGNED_FASTA}"
echo "  Trimmed alignment           : ${TRIMMED_FASTA}"
echo "  IQ-TREE prefix              : ${TREE_DIR}/${MARKER}_tree"
echo "========================================================"