# MetaXcan Paper Analysis

This repository contains software and instructions to reproduce all material and analysis in the MetaXcan paper ([Nature Communications](https://www.nature.com/articles/s41467-018-03621-1)).

You can find MetaXcan implementation in [Github](https://github.com/hakyimlab/MetaXcan).

## Abstract

<p> <em> Scalable, integrative methods to understand mechanisms that link genetic variants with phenotypes are needed. Here we derive a mathematical expression to compute PrediXcan (a gene mapping approach) results using summary data (S-PrediXcan) and show its accuracy and general robustness to misspecified reference sets. We apply this framework to 44 GTEx tissues and 100+ phenotypes from GWAS and meta-analysis studies, creating a growing public catalog of associations that seeks to capture the effects of gene expression variation on human phenotypes. Replication in an independent cohort is shown. Most of the associations are tissue specific, suggesting context specificity of the trait etiology. Colocalized significant associations in unexpected tissues underscore the need for an agnostic scanning of multiple contexts to improve our ability to detect causal regulatory mechanisms. Monogenic disease genes are enriched among significant associations for related traits, suggesting that smaller alterations of these genes may cause a spectrum of milder phenotypes. </em> </p>

## Instructions

You should build the R package `ratools`, included in a subfolder in this repository. It's readme list isntalation instructions.

Download the data release from [https://doi.org/10.5281/zenodo.1406323](https://zenodo.org/record/1406324#.W4hCCBgna90).
Unpack it into `data` subfolder so that it looks like:
```
$ tree -d -L 1 data
data
├── alternative_hypothesis
├── clinvar
├── coloc_clinvar
├── ecaviar
├── GEUVADIS
├── GIANT_HEIGHT
├── GWAS
├── images
├── metaxcan_results
├── predixcan_results
├── simulation
└── v6p_unfiltered_dbs
```

Then, knit the R markdown files in succession (You can use RStudio for that, or R command line):
- 00_figures.Rmd
- 01_predixcan_metaxcan_comparison.Rmd
- 02_Miscellaneous.Rmd
- 03_Coloc_Tern.Rmd
- 04_GEUVADIS.Rmd
- 05_twas_smr_coloc.Rmd

Finally, run:
```
#compiles figures and supplementary data
./paper_figures.py
```

## Data Availability

### Reproducible Analysis

The data strictly necessary to reproduce this paper is publicly available in zenodo:
- [Data](https://zenodo.org/record/1406324#.W4hCCBgna90), doi `10.5281/zenodo.1406323`

This release contains all method results, and supporting display information; no other data is needed for the analysis scripts in this repository.

### Analyzed Data

The methods were run on data obtained from publicly available resources:

- GTEx release version 6p, both eQTL and individual level data, obtained with dbGaP accession phs000424.v6.p1.
- De-identified genotype and phenotype information from [WTCCC](https://www.wtccc.org.uk/)
- De-identified genotype and phenotype information from [GERA](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000674.v1.p1)
- Publicly available GWAS summary statistics
The summary statistics download links  and relevant publications are included in the file `data/pheno_info.txt` from the zenodo data release (doi `10.5281/zenodo.1406323`)

