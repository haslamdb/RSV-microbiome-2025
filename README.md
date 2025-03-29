# RSV-microbiome-2025

Analysis of nasal microbiome in relation to RSV infection severity and clinical outcomes.

## Project Overview

This repository contains analysis code and workflows for studying the nasal microbiome in the context of RSV infections. The project utilizes MetaPhlAn-processed microbial profiling data to investigate relationships between nasal microbiome composition and RSV disease severity, symptoms, and clinical outcomes.

## Key Research Questions

1. How does the nasal microbiome composition differ between patients with varying RSV severity?
2. Are there specific bacterial taxa associated with more severe RSV symptoms?
3. How does the nasal microbiome change over the course of RSV infection (prior, acute, post)?
4. Can microbiome features predict clinical outcomes in RSV infections?

## Directory Structure

```
RSV-microbiome-2025/
├── data/               # Data files
│   ├── raw/            # Raw MetaPhlAn output files
│   └── processed/      # Processed abundance tables
├── metadata.csv        # Sample and subject metadata
├── notebooks/          # Jupyter notebooks for analysis
├── scripts/            # Analysis scripts
├── results/            # Analysis outputs (figures, tables)
├── tools/              # Analysis tools
│   └── metaphlan_tools/  # Git submodule for microbiome analysis tools
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

3. Install the metaphlan_tools package using conda-build:
   ```bash
   conda install conda-build
   conda develop tools/metaphlan_tools
   ```
   
   Alternatively, if you encounter issues, you can use pip with the necessary flag:
   ```bash
   pip install -e tools/metaphlan_tools --break-system-packages
   ```

4. Create necessary directories (if they don't exist):
   ```bash
   mkdir -p data/raw data/processed results/figures results/tables config
   ```

5. Copy the analysis configuration file:
   ```bash
   # If not already present
   cp config/analysis_parameters.yml.example config/analysis_parameters.yml
   # Edit the parameters as needed
   ```

### Analysis Workflow

The analysis follows these main steps:

1. **Data Processing**: Combine MetaPhlAn outputs into abundance tables
   ```bash
   python scripts/01_process_metaphlan_files.py
   ```

2. **Diversity Analysis**: Calculate and compare diversity metrics
   ```bash
   python scripts/02_calculate_diversity.py
   ```

3. **Differential Abundance**: Identify taxa that differ between clinical groups
   ```bash
   python scripts/03_identify_differential_taxa.py
   ```

4. **Longitudinal Analysis**: Track microbiome changes over the course of infection
   ```bash
   python scripts/04_analyze_longitudinal_changes.py
   ```

5. **Generate Report**: Create summary visualizations and findings
   ```bash
   python scripts/05_generate_report.py
   ```

Alternatively, explore the analysis interactively through the Jupyter notebooks in the `notebooks/` directory.

## Metadata

The metadata.csv file contains the following key variables:

- **SampleID**: Unique identifier for each sample
- **SubjectID**: Identifier for each subject (patient)
- **CollectionDate**: Date of sample collection
- **Timing**: Timing relative to infection (Prior, Acute, Post)
- **Severity**: Severity score of RSV infection
- **Symptoms**: Symptom classification (Asymptomatic, Mild, Severe)

## Using metaphlan_tools

This project uses the metaphlan_tools package (included as a Git submodule) for microbiome data analysis. The package provides functions for:

- Parsing MetaPhlAn output files
- Calculating diversity metrics
- Performing statistical analyses
- Generating visualizations

See the [metaphlan_tools README](tools/metaphlan_tools/README.md) for more details.

## Contributing

Please follow these steps to contribute to the project:

1. Create a new branch for your feature or bugfix
2. Make your changes
3. Run tests to ensure functionality
4. Submit a pull request with a clear description of the changes

## License

This project is licensed under the MIT License - see the LICENSE file for details.
