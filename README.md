# RSV-microbiome-2025

Analysis of nasal microbiome in relation to RSV infection severity and clinical outcomes.

## Project Overview

This repository contains analysis code and workflows for studying the nasal microbiome in the context of RSV infections. The project supports multiple microbiome profiling tools including MetaPhlAn, Kraken2/Bracken, and Sylph to investigate relationships between nasal microbiome composition and RSV disease severity, symptoms, and clinical outcomes.

## Directory Structure

```
RSV-microbiome-2025/
├── data/               # Data files
│   ├── raw/            # Raw profiling output files (MetaPhlAn, Kraken, Sylph)
│   └── processed/      # Processed abundance tables
├── metadata.csv        # Sample and subject metadata
├── notebooks/          # Jupyter notebooks for analysis
├── scripts/            # Analysis scripts
│   ├── bash_workflow.sh           # Complete MetaPhlAn workflow script
│   ├── metaphlan/                 # MetaPhlAn analysis scripts
│   │   ├── 01_process_metaphlan_files.py
│   │   ├── 02_calculate_diversity.py
│   │   ├── 03_differential_abundance.py
│   │   ├── 04_longitudinal_analysis.py
│   │   └── 05_generate_report.py
│   ├── kraken/                    # Kraken/Bracken analysis scripts
│   │   ├── 01_process_kraken_data.py
│   │   ├── 02_kraken_differential_abundance.py
│   │   ├── 03_generate_kraken_report.py
│   │   ├── 04_cooccurence_analysis.py
│   │   └── kraken_rsv_microbiome.sh
│   ├── sylph/                     # Sylph analysis scripts
│   │   ├── 01_sylph_parse_output.py
│   │   ├── 02_analyze_sylph_data.py
│   │   ├── 03_sylph_differential_abundance.py
│   │   ├── 04_cooccurence_analysis.py
│   │   ├── 05_generate_sylph_report.py
│   │   └── sylph_rsv_microbiome.sh
│   └── utils/          # Utility functions and parser patches
│       └── parser_patch.py
├── results/            # Analysis outputs (figures, tables)
│   ├── figures/        # Generated figures
│   │   ├── species_boxplots/      # Boxplots for significant species
│   │   └── longitudinal_plots/    # Longitudinal analysis visualizations
│   ├── kraken_analysis/           # Kraken analysis results
│   ├── sylph_analysis/            # Sylph analysis results 
│   └── tables/         # Generated data tables
├── tools/              # Analysis tools
│   ├── metaphlan_tools/  # Tools for MetaPhlAn analysis
│   ├── kraken_tools/     # Tools for Kraken/Bracken analysis
│   └── sylph_tools/      # Tools for Sylph microbiome analysis
├── config/             # Configuration files
│   └── analysis_parameters.yml    # Parameters for all analysis workflows
└── docs/               # Documentation
```

## Getting Started

### Prerequisites

- Python 3.12 or higher
- Git (with Git LFS for tracking large files)
- Conda or Mamba (recommended for environment management)
- For raw data analysis:
  - MetaPhlAn v4.0+
  - Kraken2 and Bracken
  - Sylph

### Installation

1. Clone the repository with submodules:
   ```bash
   git clone --recurse-submodules https://github.com/yourusername/RSV-microbiome-2025.git
   cd RSV-microbiome-2025
   ```

2. Create and activate the conda environment:
   ```bash
   conda env create -f environment.yml
   conda activate rsv-microbiome
   ```

3. Install the analysis tools packages:
   ```bash
   # For all analysis tools
   conda install conda-build
   conda develop tools/metaphlan_tools
   conda develop tools/kraken_tools
   conda develop tools/sylph_tools
   ```
   
   Alternatively, if you encounter issues, you can use pip:
   ```bash
   pip install -e tools/metaphlan_tools
   pip install -e tools/kraken_tools
   pip install -e tools/sylph_tools
   ```

4. Create necessary directories (if they don't exist):
   ```bash
   mkdir -p data/raw data/processed results/figures/species_boxplots results/figures/longitudinal_plots results/tables config
   ```

5. Configure the analysis parameters:
   ```bash
   # If not already present
   cp config/analysis_parameters.yml.example config/analysis_parameters.yml
   # Edit the parameters as needed
   ```

## Analysis Workflows

### MetaPhlAn Analysis

```bash
# Full MetaPhlAn workflow
./scripts/bash_workflow.sh
```

The MetaPhlAn workflow includes:
1. **Data Processing**: Combine MetaPhlAn outputs into abundance tables
2. **Diversity Analysis**: Calculate and compare diversity metrics
3. **Differential Abundance**: Identify taxa that differ between clinical groups
4. **Longitudinal Analysis**: Track microbiome changes over the course of infection
5. **Report Generation**: Create summary visualizations and findings

Individual steps can be run separately:
```bash
python scripts/metaphlan/01_process_metaphlan_files.py
python scripts/metaphlan/02_calculate_diversity.py
python scripts/metaphlan/03_differential_abundance.py
python scripts/metaphlan/04_longitudinal_analysis.py
python scripts/metaphlan/05_generate_report.py
```

### Kraken/Bracken Analysis

```bash
# Run Kraken2 and Bracken on raw sequencing data
./scripts/kraken/kraken_rsv_microbiome.sh

# Process Kraken data and run all analyses with a single command
./scripts/kraken/run_all_analyses.sh --kreport-dir path/to/kreports --bracken-dir path/to/bracken

# Or run individual analysis steps
python scripts/kraken/process_kraken_data.py --kreport-dir path/to/kreports --bracken-dir path/to/bracken
python scripts/kraken/kraken_differential_abundance.py --abundance-file results/kraken_analysis/normalized_abundance.tsv
python scripts/kraken/kraken_permanova.py --abundance-file results/kraken_analysis/normalized_abundance.tsv
python scripts/kraken/kraken_rf_shap.py --abundance-file results/kraken_analysis/normalized_abundance.tsv
python scripts/kraken/kraken_tsne.py --abundance-file results/kraken_analysis/normalized_abundance.tsv
```

The Kraken/Bracken workflow includes:
1. **Data Processing**: Run Kraken2/Bracken on raw reads and process outputs
2. **Data Normalization**: Apply CLR (centered log-ratio) or other normalization methods to the abundance data
3. **Diversity Analysis**: Calculate alpha and beta diversity metrics
4. **Differential Abundance**: Identify differentially abundant taxa between groups
5. **PERMANOVA Analysis**: Test for significant associations between metadata variables and microbiome composition
6. **Random Forest with SHAP**: Identify important taxa for predicting clinical outcomes
7. **t-SNE Visualization**: Create dimensionality-reduced visualizations of microbiome data
8. **Co-occurrence Analysis**: Study relationships between specific species

### Sylph Analysis

```bash
# Run Sylph on raw sequencing data
./scripts/sylph/sylph_rsv_microbiome.sh

# Process Sylph data and perform analyses
python scripts/sylph/01_sylph_parse_output.py
python scripts/sylph/02_analyze_sylph_data.py
python scripts/sylph/03_sylph_differential_abundance.py
python scripts/sylph/04_cooccurence_analysis.py
```

The Sylph workflow includes:
1. **Data Processing**: Parse Sylph profiling outputs for bacteria, viruses, and fungi
2. **Diversity Analysis**: Calculate diversity metrics across samples
3. **Differential Abundance**: Identify differentially abundant taxa between clinical groups
4. **Co-occurrence Analysis**: Examine specific species relationships (similar to Kraken workflow)
5. **Temporal Analysis**: Track microbiome changes over time (Prior, Acute, Post)

## Tool Comparison

The repository supports three different profiling tools, each with its strengths:

- **MetaPhlAn**: Marker gene-based approach with detailed taxonomic resolution
- **Kraken2/Bracken**: Fast k-mer based classification with good accuracy
- **Sylph**: Modern sketching-based profiler with support for bacteria, viruses, and fungi

## Metadata

The metadata.csv file contains the following key variables:

- **SampleID**: Unique identifier for each sample
- **SubjectID**: Identifier for each subject (patient)
- **CollectionDate**: Date of sample collection
- **Timing**: Timing relative to infection (Prior, Acute, Post)
- **Severity**: Severity score of RSV infection
- **Symptoms**: Symptom classification (Asymptomatic, Mild, Severe)

## Key Analyses

### Diversity Analysis

All three workflows include comprehensive diversity analysis:
- **Alpha Diversity**: Shannon, Simpson, observed species
- **Beta Diversity**: Bray-Curtis, weighted/unweighted UniFrac
- **Ordination**: PCoA and NMDS for visualization

### Differential Abundance

Identify significant differences in microbial abundance:
- Between severity groups (mild vs. severe)
- Between symptom categories
- Across time points (Prior, Acute, Post)

### Co-occurrence Analysis

The repository includes specialized scripts for co-occurrence analysis:
- Built-in focus on S. pneumoniae and H. influenzae interactions
- Cross-timepoint analysis of species co-occurrence

### Longitudinal Analysis

Track changes in the microbiome over the course of RSV infection:
- Visualize species abundance across time points
- Identify species with significant temporal patterns

## Advanced Features

### Data Normalization

The Kraken/Bracken workflow includes several normalization methods to transform abundance data for statistical analysis:

```bash
# Apply CLR normalization (default)
./scripts/kraken/run_all_analyses.sh --kreport-dir path/to/kreports --norm-method clr

# Apply relative abundance normalization (0-1 scale)
./scripts/kraken/run_all_analyses.sh --norm-method relabundance

# Apply counts per million normalization
./scripts/kraken/run_all_analyses.sh --norm-method cpm

# Apply log10 transformation
./scripts/kraken/run_all_analyses.sh --norm-method log10

# Disable normalization
./scripts/kraken/run_all_analyses.sh --no-normalize
```

Available normalization methods:
- **CLR (Centered Log-Ratio)**: Handles compositional data by log-transforming after dividing by the geometric mean (default method)
- **Relative Abundance**: Converts counts to proportions (0-1 scale)
- **CPM (Counts Per Million)**: Scales data to counts per million
- **Log10**: Log10 transformation with a small pseudocount to avoid log(0)

The normalization is automatically applied after filtering and before all downstream analyses.

### Correlation Networks

```python
# For any of the analysis tools
from metaphlan_tools.viz import plot_correlation_network
from kraken_tools.analysis.network import plot_species_correlation_network
from sylph_tools import plot_correlation_network
```

### Robust Parsing

The repository includes enhanced parsers for handling various file formats:

```python
from scripts.utils.parser_patch import patched_parse_metaphlan_file, patched_combine_samples
```

### Taxonomic Analysis

Extract and analyze taxonomy information:

```python
from kraken_tools.analysis.taxonomy import read_and_process_taxonomy
from sylph_tools import add_taxonomy_metadata
```

### PERMANOVA Analysis

The Kraken workflow includes PERMANOVA (Permutational Multivariate Analysis of Variance) to test how metadata variables explain microbiome composition:

```bash
# Run PERMANOVA analysis with specific categorical variables
python scripts/kraken/kraken_permanova.py --abundance-file results/kraken_analysis/normalized_abundance.tsv \
  --categorical-vars Timing,Severity,Symptoms --distance-metric bray --transform clr
```

PERMANOVA helps identify which factors (e.g., disease severity, timing) significantly influence microbiome structure.

### Random Forest with SHAP

The RF-SHAP analysis combines Random Forest machine learning with SHAP (SHapley Additive exPlanations) values to identify the most important taxa for prediction:

```bash
# Run RF-SHAP analysis
python scripts/kraken/kraken_rf_shap.py --abundance-file results/kraken_analysis/normalized_abundance.tsv \
  --predictors Timing,Severity,Symptoms --random-effects SubjectID
```

This analysis helps identify which bacterial species are most predictive of clinical outcomes.

### t-SNE Visualization

The t-SNE (t-Distributed Stochastic Neighbor Embedding) visualization creates low-dimensional representations of the microbiome data:

```bash
# Generate t-SNE plots
python scripts/kraken/kraken_tsne.py --abundance-file results/kraken_analysis/normalized_abundance.tsv \
  --categorical-vars Timing,Severity,Symptoms
```

t-SNE helps visualize sample clustering and relationships between samples in different clinical groups.

## Contributing

Please follow these steps to contribute to the project:

1. Create a new branch for your feature or bugfix
2. Make your changes
3. Run tests to ensure functionality
4. Submit a pull request with a clear description of the changes

## License

This project is licensed under the MIT License - see the LICENSE file for details.
