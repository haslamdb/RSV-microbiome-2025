# RSV-microbiome-2025

Analysis of nasal microbiome in relation to RSV infection severity and clinical outcomes.

## Project Overview

This repository contains analysis code and workflows for studying the nasal microbiome in the context of RSV infections. The project utilizes both MetaPhlAn and Sylph microbial profiling data to investigate relationships between nasal microbiome composition and RSV disease severity, symptoms, and clinical outcomes.

## Directory Structure

```
RSV-microbiome-2025/
├── data/               # Data files
│   ├── raw/            # Raw MetaPhlAn and Sylph output files
│   └── processed/      # Processed abundance tables
├── metadata.csv        # Sample and subject metadata
├── notebooks/          # Jupyter notebooks for analysis
├── scripts/            # Analysis scripts
│   └── utils/          # Utility functions and parser patches
├── results/            # Analysis outputs (figures, tables)
│   ├── figures/        # Generated figures
│   │   ├── species_boxplots/   # Boxplots for significant species
│   │   └── longitudinal_plots/ # Longitudinal analysis visualizations
│   └── tables/         # Generated data tables
├── tools/              # Analysis tools
│   ├── metaphlan_tools/  # Git submodule for MetaPhlAn analysis tools
│   └── sylph_tools/      # Sylph microbiome analysis tools
├── config/             # Configuration files
└── docs/               # Documentation
```

## Getting Started

### Prerequisites

- Python 3.12 or higher
- Git (with Git LFS for tracking large files)
- Conda or Mamba (recommended for environment management)

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
   # For MetaPhlAn tools
   conda install conda-build
   conda develop tools/metaphlan_tools
   
   # For Sylph tools
   conda develop tools/sylph_tools
   ```
   
   Alternatively, if you encounter issues, you can use pip with the necessary flag:
   ```bash
   pip install -e tools/metaphlan_tools --break-system-packages
   pip install -e tools/sylph_tools --break-system-packages
   ```

4. Create necessary directories (if they don't exist):
   ```bash
   mkdir -p data/raw data/processed results/figures/species_boxplots results/figures/longitudinal_plots results/tables config
   ```

5. Copy the analysis configuration file:
   ```bash
   # If not already present
   cp config/analysis_parameters.yml.example config/analysis_parameters.yml
   # Edit the parameters as needed
   ```

## Analysis Workflow

The repository supports analysis of both MetaPhlAn and Sylph microbial profiling data through specialized tool modules.

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
5. **Visualization**: Generate boxplots, ordination plots, and other visualizations
6. **Report Generation**: Create summary visualizations and findings

### Sylph Analysis

For Sylph data analysis, the repository provides a suite of tools in the `sylph_tools` module with equivalent functionality:

1. **Diversity Analysis**:
   ```python
   from sylph_tools import calculate_alpha_diversity, calculate_beta_diversity, compare_alpha_diversity
   ```

2. **Statistical Analysis**:
   ```python
   from sylph_tools import perform_permanova, differential_abundance_analysis
   ```

3. **Visualization**:
   ```python
   from sylph_tools import plot_alpha_diversity_boxplot, plot_ordination, plot_stacked_bar, plot_correlation_network
   ```

4. **Data Processing**:
   ```python
   from sylph_tools import load_metadata, preprocess_abundance_data, filter_low_abundance
   ```

### Interactive Analysis

Explore the analysis interactively through the Jupyter notebooks in the `notebooks/` directory.

## Metadata

The metadata.csv file contains the following key variables:

- **SampleID**: Unique identifier for each sample
- **SubjectID**: Identifier for each subject (patient)
- **CollectionDate**: Date of sample collection
- **Timing**: Timing relative to infection (Prior, Acute, Post)
- **Severity**: Severity score of RSV infection
- **Symptoms**: Symptom classification (Asymptomatic, Mild, Severe)

## Advanced Features

### Robust Parsing

The repository includes enhanced parsers for handling various MetaPhlAn file formats and edge cases:

```python
from scripts.utils.parser_patch import patched_parse_metaphlan_file, patched_combine_samples
```

### Network Analysis

Correlation network analysis to identify interactions between microbial species:

```python
from sylph_tools import plot_correlation_network
```

### Heatmap Visualization

Create relative abundance heatmaps with sample clustering and group annotations:

```python
from sylph_tools import plot_relative_abundance_heatmap
```

## Additional Utilities

### Taxonomy Metadata

Extract and add taxonomy information from taxon names:

```python
from sylph_tools import add_taxonomy_metadata
```

### Diversity Comparisons

Statistical comparison of alpha diversity between groups:

```python
from sylph_tools import compare_alpha_diversity
```

### Export Formats

Export abundance data to various formats, including BIOM:

```python
from sylph_tools import export_biom_format
```

## Contributing

Please follow these steps to contribute to the project:

1. Create a new branch for your feature or bugfix
2. Make your changes
3. Run tests to ensure functionality
4. Submit a pull request with a clear description of the changes

## License

This project is licensed under the MIT License - see the LICENSE file for details.
