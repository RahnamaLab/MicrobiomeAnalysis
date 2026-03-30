# MicrobiomeAnalysis
The goal of the project was to classify, at the highest taxonomic resolution possible, strains of
Fusarium spp. from crude cannabis DNA extracts. These data would be paired with mycotoxin
presence/absence and additional datapoints to further understand the role of Fusarium on
cannabis production.
# Kraken 2
Within the Kraken 2 workflow, raw reads were classified using a prebuilt database that includes the standard database plus RefSeq, protozoa, fungi, and plants.
## 1. Download Kraken2 Database

```bash
wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20260226.tar.gz
tar -xvzf k2_pluspf_20260226.tar.gz
#Check Database Content
cut -f1 library_report.tsv | sort | uniq -c

# QIIME 2
Within the QIIME 2 workflow, reads were trimmed of adapter sequences, split into independent data sets based on target amplicon (i.e., EF1α, ITS1, ITS2), quality trimmed, and processed for quality control. During quality control, paired-end reads were merged. ITS-1 and ITS-2 sequences were classified using a Naive Bayes classifier trained on fungal sequences from the UNITE database clustered at 99% identity. EF1α amplicons were classified using a Native Bayes classifier trained on the Fuseroid ID database* (version 12 March 2025) obtained from https://www.fusarium.org/page/Sequencesindatabase and filtered to only contain EF1α sequences. When the classifier was applied to the reference database to assess the best case classification of EF1α amplicons, it achieved a precision and recall of 0.98 and 0.92 (scale of 0.0-1.0), respectively, for species-level classification.
