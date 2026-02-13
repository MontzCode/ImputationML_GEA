# ImputationML_GEA

Evaluating machine learning and statistical imputation methods for handling missing data in RNA-seq gene expression datasets, implemented in R as part of a group project at Queen Mary University of London.

The study uses transcriptomic data from Mack et al. (2018) on environmental adaptation in house mice, testing eight imputation algorithms across three levels of simulated missingness (5%, 10%, 30%) and assessing their downstream impact on differential gene expression analysis.

---

## Table of Contents

- [Overview](#overview)
- [Project Structure](#project-structure)
- [Dataset](#dataset)
- [Installation](#installation)
- [Usage](#usage)
- [Imputation Methods](#imputation-methods)
- [Evaluation Pipeline](#evaluation-pipeline)
- [Key Findings](#key-findings)
- [Limitations](#limitations)

---

## Overview

Missing data in RNA-seq count matrices can introduce bias and reduce statistical power in downstream genomic analyses. This project systematically evaluates eight imputation strategies — ranging from simple statistical methods to advanced machine learning approaches — to determine which best preserves biological signal when applied to gene expression data.

The pipeline covers:

- Missingness pattern identification (Little's MCAR test, visualisation)
- Parallel imputation across all eight methods
- Kolmogorov-Smirnov (KS) test to assess distribution preservation
- PCA-based variance analysis comparing imputed data to complete cases
- Differential expression analysis (DESeq2) to evaluate downstream biological accuracy
- Gene recovery benchmarking against a predefined set of expected marker genes

---

## Project Structure

```
ImputationML_GEA/
│
├── Functions.R                        # All custom functions: preprocessing, imputation, DEA, visualisation
├── ImputationScript.R                 # Main pipeline: runs all imputations and DEA
├── Visualization.R                    # Loads saved results and generates all evaluation plots
├── LICENSE
│
├── data/                              # Input data (not tracked)
│   ├── countdata5.txt                 # Gene expression counts with 5% missing data
│   ├── countdata10.txt                # Gene expression counts with 10% missing data
│   ├── countdata30.txt                # Gene expression counts with 30% missing data
│   └── metadata.txt                  # Sample metadata including latitude and physiological traits
│
├── ImputedDatasets/                   # Imputed output CSVs (auto-generated, not tracked)
│   ├── Imputed_Median_30%.csv
│   ├── Imputed_Mean_30%.csv
│   ├── Imputed_KNN_30%.csv
│   ├── Imputed_SVD_30%.csv
│   ├── Imputed_EM_30%.csv
│   ├── Imputed_MissForest_30%.csv
│   ├── Imputed_BPCA_30%.csv
│   └── Imputed_pca_30%.csv
│
└── Results/
    └── DEA_Results/                   # Differential expression results per method (auto-generated, not tracked)
        ├── results_Median_30%.csv
        ├── results_KNN_30%.csv
        └── ...
```

---

## Dataset

The project uses RNA-seq count data from:

> Mack KL, Ballinger MA, Phifer-Rixey M, Nachman MW. *Gene regulation underlies environmental adaptation in house mice.* Genome Research. 2018;28(11):1636–45.

Samples are drawn from Northern (latitude > 40°) and Southern (latitude < 35°) mouse populations along the east coast of North America. The goal of the differential expression analysis is to identify genes up- or down-regulated between these two groups.

Three versions of the count matrix were used, each with a different level of simulated missingness: 5%, 10% and 30%. Missingness was confirmed as Missing Completely At Random (MCAR) via Little's test and visual inspection.

A set of 11 biologically expected marker genes was defined based on the original study findings and used as a benchmark for evaluating gene recovery:

```r
GenesExpected = c("Lcn12", "Nlrp12", "ENSMUSG00000084349", "ENSMUSG00000089698",
                  "Mup8", "ENSMUSG00000081684", "ENSMUSG000000110453", "Sdk2",
                  "Cyp2d11", "ENSMUSG00000078965", "ENSMUSG00000106576")
```

---

## Installation

This project requires R (≥ 4.0) and Bioconductor.

```r
# Install CRAN packages
install.packages(c(
  "data.table", "dplyr", "tidyverse", "tibble", "magrittr",
  "parallel", "doParallel", "ggplot2", "gridExtra",
  "metan", "VIM", "missForest", "mice", "imputeR",
  "naniar", "na.tools", "missMDA", "caret"
))

# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "impute", "pcaMethods", "edgeR",
                       "EnhancedVolcano", "gprofiler2"))
```

---

## Usage

### Step 1 — Run imputations and DEA

Open `ImputationScript.R` and set the input dataset to the missingness level you want to analyse:

```r
countdata <- read.table("data/countdata30.txt", header = T, row.names = 1)
```

Run the full script. All eight imputations will execute in sequence (computationally intensive methods use parallel processing). Uncomment the `write.csv` lines to save outputs to `ImputedDatasets/` and `Results/DEA_Results/`.

### Step 2 — Visualise results

Open `Visualization.R`, update the import paths to match the missingness level you analysed, and run the script to generate:

- KS test distribution plots per imputation method
- PCA variance comparison across all methods vs complete cases
- Volcano plots from DESeq2 DEA
- Bar plot of expected gene recovery per method

---

## Imputation Methods

All imputation functions are defined in `Functions.R`. Data is transposed before imputation so that genes are treated as features (columns) and samples as observations (rows).

| Method | Function | Package | Notes |
|---|---|---|---|
| Median | `impute_median()` | `caret` | `preProcess` with `medianImpute` |
| Mean | `impute_mean()` | base R | Column-wise mean replacement |
| KNN | `impute_knn()` | `impute` (Bioconductor) | k=2 nearest neighbours |
| MissForest | `impute_MissForest()` | `missForest` | Iterative random forest, parallelised by variable |
| BPCA | `impute_BPCA()` | `pcaMethods` | Bayesian PCA with 6 principal components |
| SVD (EM-SVD) | `impute_svd()` | `metan` | EM-SVD algorithm, 20,000 iterations |
| EM-AMMI | `impute_em()` | `metan` | EM-AMMI statistical model, 20,000 iterations |
| PCA (missMDA) | `impute_pca()` | `missMDA` | `imputePCA` with 5 components |

Parallel processing is configured automatically based on available physical cores:

```r
physical_cores <- detectCores(logical = FALSE)
cl <- makeCluster(physical_cores - 1)
registerDoParallel(cl)
```

Simple methods (mean, median, KNN) do not use parallel computation given their speed.

---

## Evaluation Pipeline

**1. Missingness pattern analysis**

Little's MCAR test confirms that missing values are distributed independently of observed data. Visual inspection via `vis_miss()` supports this across all three missingness levels.

**2. Kolmogorov-Smirnov test**

The `ks_test()` function compares the distribution of each imputed gene column against the original observed values. Results are expressed as the percentage of genes with non-significant differences (p > 0.05), where higher is better.

**3. PCA variance analysis**

`VariancePCAGraph()` runs PCA on each imputed dataset after TMM normalisation (via edgeR) and plots the first two principal components for all methods alongside complete cases, allowing visual assessment of how well each method preserves the global data structure.

**4. Differential gene expression analysis**

`DEA()` runs DESeq2 on each imputed dataset using `~ Sex + group` as the design formula. Significant genes (padj < 0.05) are converted from Ensembl IDs to gene names via `gprofiler2::gconvert()` and compared against the expected gene list.

**5. Gene recovery benchmarking**

`compile_gene_list()` and `visualize_DEA()` compile all significant genes per method and generate a stacked bar plot showing how many expected marker genes each imputation method successfully recovered.

---

## Key Findings

Median and KNN imputation consistently outperformed all other methods across all three missingness levels.

At 5% and 10% missingness, most methods performed reasonably well in both KS score and gene recovery. The separation between methods became most apparent at 30% missingness, where only median and KNN continued to identify meaningful genes in downstream DEA. Methods like MissForest, BPCA and SVD, despite their complexity, showed significant degradation at higher sparsity levels.

Median imputation identified one additional expected gene compared to KNN at 30% missingness, recovering Cyp2d11 where KNN did not. This was the clearest demonstration that, for high-dimensional sparse genomic data under substantial data loss, simpler imputation strategies can outperform iterative and model-based approaches.

---

## Limitations

- **No ground truth**: Without the original complete dataset, imputation accuracy relies on statistical proxies (KS test, PCA) rather than direct comparison to true values
- **Computational cost**: BPCA and SVD required 20-45 minutes per run even with parallel processing, making them impractical for large-scale genomic pipelines
- **Keyword-based labelling**: Not applicable here, but missingness was simulated rather than naturally occurring, which may not fully represent real experimental data loss patterns
- **Small sample size**: The relatively small number of samples in the Mack et al. dataset limits how much nearest-neighbour and forest-based methods can leverage local structure

---

## Reference

Mack KL, Ballinger MA, Phifer-Rixey M, Nachman MW. Gene regulation underlies environmental adaptation in house mice. *Genome Research*. 2018 Nov 1;28(11):1636–45.
