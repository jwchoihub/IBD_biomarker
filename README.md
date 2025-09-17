# IBD_biomarker

### Introduction
Welcome to the repository for the analysis codes used in our study on gut microbiota-based biomarkers for inflammatory bowel disease (IBD), including Crohn’s disease (CD) and ulcerative colitis (UC).

### Overview
This repository contains scripts and resources for analyzing and visualizing gut microbiota data obtained from a large-scale, multi-center cohort.
Our research aimed to identify reproducible and biologically meaningful microbial features capable of classifying IBD and its subtypes using fecal 16S rRNA gene sequencing data.

### Data and Methods
- **Sample Collection**: Fecal samples were collected from participants across multiple centers. DNA was extracted and subjected to 16S rRNA gene sequencing.
- **Data Processing**: Raw sequence data were processed using QIIME2 and DADA2 to obtain high-quality ASVs (amplicon sequence variants). Raw FASTQ files are processed in QIIME2 to generate QIIME2 artifacts (e.g., `feature-table.qza`, `taxonomy.qza`). 
- **Analysis**: The downstream analysis was conducted in R. Import the resulting `.qza` files into R and convert to a `phyloseq` object.

### Note
Some parts of the code (e.g., detailed configurations or parameters) are not publicly available at this time, as the manuscript is currently under submission. These components will be released after publication to ensure proper citation and responsible reuse. In addition, we plan to provide a small test dataset (processed genus-level abundance table and metadata) together with example scripts to facilitate reproducibility and allow users to replicate the key analyses.
For academic or review-related inquiries, please feel free to contact us.
