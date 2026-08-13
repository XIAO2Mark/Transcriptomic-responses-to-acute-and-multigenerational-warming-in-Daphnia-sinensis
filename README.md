# Transcriptomic responses to acute and multigenerational warming in Daphnia sinensis

This repository contains the R scripts used for RNA-seq processing, differential-expression analysis, expression-pattern classification, concordant/divergent response analysis, variance-sensitivity analysis, transcription-factor regulatory-complexity analysis, and generation of the study figures.

## Scripts

- `00_helpers.R`: shared functions and input/output settings
- `01_primary_deseq2.R`: primary differential-expression analysis
- `02_expression_patterns.R`: classification of expression patterns
- `03_concordant_divergent.R`: concordant and divergent response analysis
- `04_variance_sensitivity.R`: sensitivity analysis excluding highly variable genes
- `05_tf_regulatory_complexity.R`: predicted TF regulatory-complexity analysis
- `99_session_info.R`: R session information and package versions

Raw RNA-seq data are available from the NCBI Sequence Read Archive under BioProject accession no. PRJNAxxxxxx.
