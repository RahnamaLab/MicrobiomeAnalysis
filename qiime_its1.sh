#!/bin/bash
#SBATCH --job-name=qiime_ITS1

set -euo pipefail


###############################################################################
# INPUT FILES
###############################################################################

MANIFEST="manifest_fixed.tsv"
METADATA="metadata.tsv"

###############################################################################
# ITS1 SETTINGS
###############################################################################

MARKER="ITS1"

F_PRIMER="ACCTGCGGARGGATCA"
R_PRIMER="GAGATCCRTTGYTRAAAGTT"
F_PRIMER_RC="TGATCCYTCCGCAGGT"      # revcomp of ACCTGCGGARGGATCA
R_PRIMER_RC="AACTTTYARCAAYGGATCTC"  # revcomp of GAGATCCRTTGYTRAAAGTT

TRIM_LEFT_F=0
TRIM_LEFT_R=0
TRUNC_LEN_F=200
TRUNC_LEN_R=200

MIN_FEATURE_FREQ=5
MIN_SEQ_LEN=100

THREADS="${SLURM_CPUS_PER_TASK:-12}"

###############################################################################
# REFERENCE DATABASE
###############################################################################

UNITE_FASTA="UNITE_QIIME/sh_refs_qiime_ver10_dynamic_19.02.2025.fasta"
UNITE_TAXONOMY="UNITE_QIIME/sh_taxonomy_qiime_ver10_dynamic_19.02.2025.txt"

CLASSIFIER_DIR="classifiers"
mkdir -p "${CLASSIFIER_DIR}"

CLASSIFIER="${CLASSIFIER_DIR}/unite_its1_classifier.qza"

###############################################################################
# SOFTWARE
###############################################################################

conda activate /home/software/qiime2-2023.5

spack load clustal-omega
spack load trimal
spack load iq-tree@2.3.6

###############################################################################
# BUILD ITS1 CLASSIFIER
###############################################################################

if [[ ! -f "${CLASSIFIER}" ]]; then

    WORKDIR="${CLASSIFIER_DIR}/ITS1_build"
    mkdir -p "${WORKDIR}"

    qiime tools import \
      --type 'FeatureData[Sequence]' \
      --input-path "${UNITE_FASTA}" \
      --output-path "${WORKDIR}/ref_seqs.qza"

    qiime tools import \
      --type 'FeatureData[Taxonomy]' \
      --input-format HeaderlessTSVTaxonomyFormat \
      --input-path "${UNITE_TAXONOMY}" \
      --output-path "${WORKDIR}/taxonomy.qza"

    qiime feature-classifier extract-reads \
      --i-sequences "${WORKDIR}/ref_seqs.qza" \
      --p-f-primer "${F_PRIMER}" \
      --p-r-primer "${R_PRIMER}" \
      --p-trunc-len 250 \
      --o-reads "${WORKDIR}/trimmed_refs.qza"

    qiime feature-classifier fit-classifier-naive-bayes \
      --i-reference-reads "${WORKDIR}/trimmed_refs.qza" \
      --i-reference-taxonomy "${WORKDIR}/taxonomy.qza" \
      --o-classifier "${CLASSIFIER}"
fi

###############################################################################
# OUTPUT FILES
###############################################################################

DEMUX_QZA="ITS1-paired-end-demux.qza"
DEMUX_QZV="ITS1-paired-end-demux.qzv"

TRIMMED_QZA="ITS1-trimmed-paired-end-demux.qza"
TRIMMED_QZV="ITS1-trimmed-paired-end-demux.qzv"

TABLE_QZA="ITS1-trimmed-table.qza"
REP_QZA="ITS1-trimmed-rep-seqs.qza"
STATS_QZA="ITS1-trimmed-denoising-stats.qza"

TABLE_QZV="ITS1-trimmed-table.qzv"
REP_QZV="ITS1-trimmed-rep-seqs.qzv"
STATS_QZV="ITS1-trimmed-denoising-stats.qzv"

FILTER_TABLE="ITS1-retained-hits-table.qza"
FILTER_REP="ITS1-retained-hits-rep-seqs.qza"

LEN_TABLE="ITS1-retained-more-${MIN_SEQ_LEN}-table.qza"
LEN_REP="ITS1-retained-more-${MIN_SEQ_LEN}-rep-seqs.qza"

TAXONOMY_QZA="ITS1-taxonomy.qza"
TAXONOMY_QZV="ITS1-taxonomy.qzv"
BARPLOT_QZV="ITS1-taxa-barplot.qzv"

TREE_DIR="iqtree_ITS1"
mkdir -p "${TREE_DIR}"

TREE_FASTA="${TREE_DIR}/ITS1-sequences.fasta"
ALIGN_FASTA="${TREE_DIR}/ITS1_aligned.fa"
TRIM_FASTA="${TREE_DIR}/ITS1_trimmed.fa"

###############################################################################
# STEP 1 IMPORT
###############################################################################

qiime tools import \
 --type 'SampleData[PairedEndSequencesWithQuality]' \
 --input-path "${MANIFEST}" \
 --input-format PairedEndFastqManifestPhred33V2 \
 --output-path "${DEMUX_QZA}"

qiime demux summarize \
 --i-data "${DEMUX_QZA}" \
 --o-visualization "${DEMUX_QZV}"

###############################################################################
# STEP 2 CUTADAPT
###############################################################################

qiime cutadapt trim-paired \
 --i-demultiplexed-sequences "${DEMUX_QZA}" \
 --p-front-f "${F_PRIMER}" \
 --p-front-r "${R_PRIMER}" \
 --p-adapter-f "${R_PRIMER_RC}" \
 --p-adapter-r "${F_PRIMER_RC}" \
 --p-times 2 \
 --p-error-rate 0.1 \
 --p-overlap 3 \
 --p-discard-untrimmed \
 --p-match-adapter-wildcards True \
 --p-minimum-length "${MIN_SEQ_LEN}" \
 --o-trimmed-sequences "${TRIMMED_QZA}"

qiime demux summarize \
 --i-data "${TRIMMED_QZA}" \
 --o-visualization "${TRIMMED_QZV}"

###############################################################################
# STEP 3 DADA2
###############################################################################

qiime dada2 denoise-paired \
 --i-demultiplexed-seqs "${TRIMMED_QZA}" \
 --p-trim-left-f "${TRIM_LEFT_F}" \
 --p-trim-left-r "${TRIM_LEFT_R}" \
 --p-trunc-len-f "${TRUNC_LEN_F}" \
 --p-trunc-len-r "${TRUNC_LEN_R}" \
 --p-n-threads "${THREADS}" \
 --o-table "${TABLE_QZA}" \
 --o-representative-sequences "${REP_QZA}" \
 --o-denoising-stats "${STATS_QZA}"

###############################################################################
# STEP 4 FILTER FEATURES
###############################################################################

qiime feature-table filter-features \
 --i-table "${TABLE_QZA}" \
 --p-min-frequency "${MIN_FEATURE_FREQ}" \
 --o-filtered-table "${FILTER_TABLE}"

qiime feature-table filter-seqs \
 --i-data "${REP_QZA}" \
 --i-table "${FILTER_TABLE}" \
 --o-filtered-data "${FILTER_REP}"

###############################################################################
# STEP 5 LENGTH FILTER
###############################################################################

TMP_EXPORT="${TREE_DIR}/export"

mkdir -p "${TMP_EXPORT}"

qiime tools export \
 --input-path "${FILTER_REP}" \
 --output-path "${TMP_EXPORT}"

awk -v minlen="${MIN_SEQ_LEN}" '
BEGIN{RS=">";FS="\n"}
NR>1{
seq=""
for(i=2;i<=NF;i++)seq=seq $i
if(length(seq)>=minlen){
print ">"$1
print seq
}
}' "${TMP_EXPORT}/dna-sequences.fasta" > "${TREE_DIR}/filtered.fasta"

qiime tools import \
 --type 'FeatureData[Sequence]' \
 --input-path "${TREE_DIR}/filtered.fasta" \
 --output-path "${LEN_REP}"

qiime feature-table filter-features \
 --i-table "${FILTER_TABLE}" \
 --m-metadata-file "${LEN_REP}" \
 --o-filtered-table "${LEN_TABLE}"

###############################################################################
# STEP 6 TAXONOMY
###############################################################################

qiime feature-classifier classify-sklearn \
 --i-classifier "${CLASSIFIER}" \
 --i-reads "${LEN_REP}" \
 --o-classification "${TAXONOMY_QZA}"

qiime metadata tabulate \
 --m-input-file "${TAXONOMY_QZA}" \
 --o-visualization "${TAXONOMY_QZV}"

qiime taxa barplot \
 --i-table "${LEN_TABLE}" \
 --i-taxonomy "${TAXONOMY_QZA}" \
 --m-metadata-file "${METADATA}" \
 --o-visualization "${BARPLOT_QZV}"

###############################################################################
# STEP 7 EXPORT FASTA
###############################################################################

qiime tools export \
 --input-path "${LEN_REP}" \
 --output-path "${TREE_DIR}"

mv "${TREE_DIR}/dna-sequences.fasta" "${TREE_FASTA}"


echo "ITS1 analysis completed successfully."

