# MicrobiomeAnalysis
The goal of this project is to classify Fusarium species at the highest possible taxonomic resolution from crude cannabis DNA extracts. These results are integrated with mycotoxin presence/absence and additional metadata to better understand the role of Fusarium in cannabis production.

# Data
You can download the data from this link.*****

# Kraken 2
Within the Kraken2 workflow, raw reads are taxonomically classified using a prebuilt database that includes the standard database along with RefSeq, protozoa, fungi, and plant sequences.

```bash

## 1. Download Kraken2 Database
wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20260226.tar.gz
tar -xvzf k2_pluspf_20260226.tar.gz

# Check database content
cut -f1 library_report.tsv | sort | uniq -c

## 2. Run Kraken2
sbatch kraken2.sh

```
# QIIME 2

Within the QIIME 2 workflow, reads are first trimmed to remove adapter and primer sequences, then separated into independent datasets based on target amplicons (EF1α, ITS1, and ITS2). Reads are quality-filtered, and paired-end reads are merged during denoising.

ITS1 and ITS2 sequences are classified using a Naive Bayes classifier trained on fungal sequences from the UNITE database clustered at 99% identity.

EF1α amplicons are classified using a Naive Bayes classifier trained on the Fusarioid ID database (latest version), filtered to include only EF1α sequences:
https://www.fusarium.org/page/Sequencesindatabase

Evaluation of the EF1α classifier on the reference dataset achieved a precision of 0.98 and recall of 0.92 for species-level classification.

## Phylogenetic Tree Construction

Phylogenetic trees were generated from the filtered representative ASV sequences obtained after QIIME 2 processing. For the combined analysis, sequences from all samples were aligned using Clustal Omega, and poorly aligned regions were removed with trimAl. A maximum likelihood phylogenetic tree was then constructed using IQ-TREE2 with automatic model selection (ModelFinder) and branch support assessed using ultrafast bootstrap and SH-aLRT tests (1000 replicates each). In addition to the combined dataset, trees were also constructed separately for each sample by first filtering sample-specific ASVs and repeating the same alignment and tree-building procedure. This approach enables both global comparison of diversity across all samples and detailed evaluation of within-sample variation.

```bash

## 1. Download EF1α Database
wget http://www.fusarium.org/images/alignments/sequencesindatabase10022026.zip
unzip sequencesindatabase10022026.zip

## 2. Download ITS Database (UNITE)
mkdir UNITE_QIIME && cd UNITE_QIIME
# Download: sh_qiime_release_19.02.2025.tgz from:
# https://doi.plutof.ut.ee/doi/10.15156/BIO/3301241
tar -xvzf sh_qiime_release_19.02.2025.tgz
cd ..

## 3. Run QIIME2 (all samples)
sbatch qiime2.sh EF1
sbatch qiime2.sh ITS1
sbatch qiime2.sh ITS2

## 4. Run QIIME2 (per sample)
sbatch qiime2_per_sample.sh EF1
sbatch qiime2_per_sample.sh ITS1
sbatch qiime2_per_sample.sh ITS2

```
