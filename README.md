# RSV-Microbiome-Analysis

Analysis of nasal microbiome in relation to RSV infection severity and clinical outcomes.

## Project Overview

This repository contains analysis code and workflows for studying the nasal microbiome in the context of RSV infections. The project supports multiple microbiome profiling tools including Kraken2/Bracken and Sylph to investigate relationships between nasal microbiome composition and RSV disease severity, symptoms, and clinical outcomes across different time points (Prior, Acute, Post).

## Directory Structure

```
RSV-Microbiome-Analysis/
├── data/                # Data files
│   ├── raw/             # Raw profiling output files (Kraken, Sylph)
│   └── processed/       # Processed abundance tables
├── metadata.csv         # Sample and subject metadata
├── scripts/             # Analysis scripts
│   ├── kraken/          # Kraken/Bracken analysis scripts
│   │   ├── cooccurence_analysis.py             # Co-occurrence patterns analysis
│   │   ├── feature_selection.py                # Variable importance analysis
│   │   ├── fixed/                              # Fixed versions of analysis scripts
│   │   ├── kraken_differential_abundance.py    # Differential abundance testing
│   │   ├── kraken_permanova.py                 # PERMANOVA analysis
│   │   ├── kraken_rf_shap.py                   # Random Forest with SHAP analysis
│   │   ├── kraken_rsv_microbiome.sh            # Kraken2/Bracken workflow script
│   │   ├── kraken_tsne.py                      # t-SNE visualization
│   │   ├── process_kraken_data.py              # Process raw Kraken/Bracken data
│   │   └── rf_shap_fixed.py                    # Fixed RF-SHAP implementation
│   ├── sylph/           # Sylph analysis scripts
│   │   ├── 01_sylph_parse_output.py            # Parse Sylph output files
│   │   ├── 02_analyze_sylph_data.py            # General Sylph data analysis
│   │   ├── 03_sylph_differential_abundance.py  # Differential abundance testing
│   │   ├── 04_cooccurence_analysis.py          # Co-occurrence analysis
│   │   ├── 05_generate_sylph_report.py         # Report generation
│   │   └── sylph_rsv_microbiome.sh             # Sylph profiling workflow
│   └── utils/           # Utility functions and helpers
├── results/             # Analysis outputs (figures, tables)
│   ├── figures/         # Generated figures
│   │   ├── species_boxplots/      # Boxplots for significant species
│   │   └── longitudinal_plots/    # Longitudinal analysis visualizations
│   ├── kraken_analysis/           # Kraken analysis results
│   │   ├── filtered_kraken_s_abundance.tsv     # Filtered species abundance
│   │   ├── tables/                # Statistical results tables
│   │   └── figures/               # Visualization figures
│   ├── sylph_analysis/            # Sylph analysis results 
│   └── tables/          # Generated data tables
├── config/              # Configuration files
└── docs/                # Documentation
```

## Getting Started

### Prerequisites

- Python 3.12 or higher
- Git
- Required Python packages (pandas, numpy, matplotlib, seaborn, scikit-bio, scikit-learn, shap, statsmodels)
- For raw data analysis:
  - Kraken2 and Bracken
  - Sylph

### Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/RSV-Microbiome-Analysis.git
   cd RSV-Microbiome-Analysis
   ```

2. Create and activate a conda environment:
   ```bash
   conda create -n rsv-microbiome python=3.9 pandas numpy matplotlib seaborn scikit-learn statsmodels
   conda activate rsv-microbiome
   pip install scikit-bio shap biom-format
   ```

3. Create necessary directories (if they don't exist):
   ```bash
   mkdir -p data/raw data/processed results/kraken_analysis/{tables,figures} results/sylph_analysis
   ```

## Kraken2/Bracken Analysis Workflow

The Kraken2/Bracken analysis workflow consists of several steps, from processing raw sequencing data to advanced statistical analysis.

### 1. Running Kraken2 and Bracken

The `kraken_rsv_microbiome.sh` script runs Kraken2 and Bracken on raw sequencing data:

```bash
# Run Kraken2 and Bracken on samples
./scripts/kraken/kraken_rsv_microbiome.sh
```

This script:
- Processes each sample with Kraken2 to assign taxonomic labels
- Runs Bracken to improve abundance estimation at species, genus, and family levels
- Outputs classification results (.kraken files) and reports (.kreport files)
- Generates abundance estimates at species, genus, and family levels

### 2. Processing Kraken/Bracken Data

The `process_kraken_data.py` script combines individual sample files into a unified abundance table:

```bash
python scripts/kraken/process_kraken_data.py \
  --kreport-dir data/raw/KrakenReports \
  --output-dir results/kraken_analysis \
  --metadata metadata.csv \
  --taxonomic-level S \
  --min-abundance 0.01 \
  --min-prevalence 0.1 \
  --normalization clr
```

Key features:
- Combines abundance data from multiple samples
- Filters out human reads
- Filters low-abundance and low-prevalence taxa
- Applies normalization (CLR, rarefaction, or other methods)
- Joins with metadata for downstream analysis
- Creates abundance tables at different taxonomic levels

### 3. Co-occurrence Analysis

The `cooccurence_analysis.py` script analyzes relationships between specific bacterial species:

```bash
python scripts/kraken/cooccurence_analysis.py \
  --input-file results/kraken_analysis/filtered_kraken_s_abundance.tsv \
  --metadata metadata.csv \
  --output-dir results/kraken_analysis/cooccurrence \
  --target-species "Streptococcus pneumoniae,Haemophilus influenzae,Moraxella catarrhalis"
```

This analysis:
- Examines co-occurrence patterns between specified species (especially S. pneumoniae and H. influenzae)
- Analyzes abundance patterns across time points (Prior, Acute, Post)
- Compares abundance between severity and symptom groups
- Generates visualizations for species correlations
- Creates faceted boxplots for each species across time points and clinical groups

### 4. Differential Abundance Analysis

The `kraken_differential_abundance.py` script identifies taxa that differ significantly between groups:

```bash
python scripts/kraken/kraken_differential_abundance.py \
  --abundance-file results/kraken_analysis/filtered_kraken_s_abundance.tsv \
  --output-dir results/kraken_analysis/differential_abundance \
  --group-col Timing \
  --method kruskal \
  --p-threshold 0.05
```

Features:
- Performs statistical testing (Kruskal-Wallis, ANOVA, t-test)
- Applies multiple testing correction
- Generates boxplots for significant species
- Supports different grouping variables (Timing, Severity, Symptoms)

### 5. PERMANOVA Analysis

The `kraken_permanova.py` script tests for significant associations between metadata variables and microbiome composition:

```bash
python scripts/kraken/kraken_permanova.py \
  --abundance-file results/kraken_analysis/filtered_kraken_s_abundance.tsv \
  --metadata metadata.csv \
  --output-dir results/kraken_analysis/permanova \
  --categorical-vars Timing,Severity,Symptoms \
  --distance-metric bray \
  --transform clr
```

This script:
- Calculates beta diversity distances between samples
- Performs PERMANOVA tests for each categorical variable
- Generates ordination plots (PCoA and NMDS)
- Creates summary reports with p-values and R² values

### 6. Random Forest with SHAP Analysis

The `kraken_rf_shap.py` script identifies important variables using Random Forest and SHAP values:

```bash
python scripts/kraken/kraken_rf_shap.py \
  --abundance-file results/kraken_analysis/filtered_kraken_s_abundance.tsv \
  --metadata metadata.csv \
  --output-dir results/kraken_analysis/rf_shap \
  --target-taxa "top" \
  --predictors Timing,Severity,Symptoms \
  --random-effects SubjectID
```

Features:
- Trains Random Forest models to predict microbe abundance
- Calculates SHAP values to explain model predictions
- Identifies the most important clinical variables for each taxon
- Fits mixed-effects models to account for repeated measurements
- Generates SHAP summary plots and feature importance visualizations

### 7. Feature Selection

The `feature_selection.py` script identifies clinical variables associated with microbiome differences:

```bash
python scripts/kraken/feature_selection.py \
  --abundance-file results/kraken_analysis/filtered_kraken_s_abundance.tsv \
  --metadata metadata.csv \
  --output-dir results/kraken_analysis/feature_selection \
  --predictors "Age,Gender,Severity,Symptoms" \
  --distance-metric bray \
  --transform clr
```

This analysis:
- Creates pairwise differences between samples for both microbiome and metadata
- Trains a Random Forest model to predict microbiome distances from metadata differences
- Identifies which variables best explain differences in microbiome composition
- Generates visualizations of feature importance

### 8. t-SNE Visualization

The `kraken_tsne.py` script creates dimensionality-reduced visualizations of microbiome data:

```bash
python scripts/kraken/kraken_tsne.py \
  --abundance-file results/kraken_analysis/filtered_kraken_s_abundance.tsv \
  --metadata metadata.csv \
  --output-dir results/kraken_analysis/tsne \
  --categorical-vars Timing,Severity,Symptoms \
  --transform clr \
  --perplexity 30
```

Features:
- Creates t-SNE visualizations to reveal patterns in high-dimensional data
- Colors points by different metadata variables
- Helps visualize clustering of samples by clinical factors

## Key Analyses and Findings

### Diversity Analysis

The Kraken workflow includes comprehensive diversity analysis:
- **Alpha Diversity**: Shannon, Simpson, observed species
- **Beta Diversity**: Bray-Curtis, Jaccard, Euclidean distances
- **Ordination**: PCoA and NMDS for visualization

### Differential Abundance

Identify significant differences in microbial abundance:
- Between severity groups (mild vs. severe)
- Between symptom categories
- Across time points (Prior, Acute, Post)

### Co-occurrence Analysis

Special focus on relationships between specific bacterial species:
- S. pneumoniae and H. influenzae interactions
- M. catarrhalis relationships with other species
- Cross-timepoint analysis of species co-occurrence
- Correlation networks of interacting species

### Longitudinal Analysis

Track changes in the microbiome over the course of RSV infection:
- Visualize species abundance across time points
- Create subject trajectory plots to track individual changes
- Analyze changes from Prior to Acute and from Acute to Post timepoints
- Identify species with significant temporal patterns

## Advanced Features

### Data Normalization

The Kraken/Bracken workflow includes several normalization methods to transform abundance data for statistical analysis:

```bash
# Apply CLR normalization (default)
python scripts/kraken/process_kraken_data.py --normalization clr

# Apply rarefaction normalization (for count data)
python scripts/kraken/process_kraken_data.py --normalization rarefaction --rarefaction-depth 10000
```

Available normalization methods:
- **CLR (Centered Log-Ratio)**: Handles compositional data by log-transforming after dividing by the geometric mean
- **Rarefaction**: Subsamples to equal depth across samples
- **Relative Abundance**: Converts counts to percentages

### Mixed-Effects Modeling

The RF-SHAP analysis includes mixed-effects modeling to account for repeated measurements:

```bash
python scripts/kraken/kraken_rf_shap.py --random-effects SubjectID
```

This helps account for subject-specific effects when analyzing longitudinal data.

### PERMANOVA Analysis

PERMANOVA helps identify which factors (e.g., disease severity, timing) significantly influence microbiome structure:

```bash
python scripts/kraken/kraken_permanova.py --categorical-vars Timing,Severity,Symptoms
```

The analysis reports R² values (effect sizes) and p-values for each variable.

### Random Forest with SHAP

The RF-SHAP analysis helps identify which bacterial species are most predictive of clinical outcomes:

```bash
python scripts/kraken/kraken_rf_shap.py --predictors Timing,Severity,Symptoms
```

SHAP values provide both global importance and the direction of effect for each variable.

## Sylph Analysis

The repository also includes similar analysis capabilities for Sylph, a modern sketching-based profiler with support for bacteria, viruses, and fungi:

```bash
# Run Sylph profiling
./scripts/sylph/sylph_rsv_microbiome.sh

# Process Sylph outputs
python scripts/sylph/01_sylph_parse_output.py

# Run analyses
python scripts/sylph/02_analyze_sylph_data.py
python scripts/sylph/03_sylph_differential_abundance.py
python scripts/sylph/04_cooccurence_analysis.py
```

The Sylph workflow provides similar capabilities to the Kraken workflow, with the added benefit of viral and fungal profiling.

## License

This project is licensed under the MIT License - see the LICENSE file for details.
