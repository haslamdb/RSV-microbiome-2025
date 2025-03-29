# RSV-microbiome-2025 Project Structure

This document outlines the recommended directory structure for the RSV-microbiome-2025 project, designed to facilitate microbiome analysis using the metaphlan_tools package.

## Directory Structure

```
RSV-microbiome-2025/
├── README.md                     # Project overview, goals, and usage instructions
├── metadata.csv                  # Subject and sample metadata
├── requirements.txt              # Project-specific Python dependencies
├── .gitignore                    # Files to exclude from git
├── environment.yml               # Conda environment specification
│
├── tools/                        # Analysis tools
│   └── metaphlan_tools/          # Git submodule for metaphlan_tools package
│
├── data/                         # Data directory (consider gitignoring large files)
│   ├── raw/                      # Raw MetaPhlAn output files
│   │   └── README.md             # Description of data sources and format
│   ├── processed/                # Processed abundance tables
│   └── external/                 # External data resources (e.g., reference databases)
│
├── notebooks/                    # Jupyter notebooks for exploratory analysis
│   ├── 01_data_preprocessing.ipynb
│   ├── 02_diversity_analysis.ipynb
│   ├── 03_differential_abundance.ipynb
│   └── 04_longitudinal_analysis.ipynb
│
├── scripts/                      # Analysis scripts
│   ├── 01_process_metaphlan_files.py
│   ├── 02_calculate_diversity.py
│   ├── 03_identify_differential_taxa.py
│   ├── 04_analyze_longitudinal_changes.py
│   ├── 05_generate_report.py
│   └── utils/                    # Utility scripts
│       ├── __init__.py
│       └── helpers.py            # Project-specific helper functions
│
├── results/                      # Analysis outputs
│   ├── figures/                  # Generated plots and visualizations
│   ├── tables/                   # Statistical results and data tables
│   └── reports/                  # Analysis reports and summaries
│
├── docs/                         # Documentation
│   ├── workflow.md               # Analysis workflow documentation
│   ├── data_dictionary.md        # Description of metadata variables
│   └── submodule_workflow.md     # Copy of the submodule workflow document
│
└── config/                       # Configuration files
    ├── analysis_parameters.yml   # Parameters for analysis scripts
    └── visualization_settings.yml # Settings for plot appearance
```

## Key Components

### Data Management

- **data/raw/**: Store original MetaPhlAn output files here
- **data/processed/**: Store combined abundance tables and transformed data
- **metadata.csv**: Keep centralized sample metadata in the project root

### Analysis Workflow

The project follows a sequential analysis workflow:

1. **Data Processing**: Convert raw MetaPhlAn outputs to abundance tables
2. **Diversity Analysis**: Calculate alpha and beta diversity metrics
3. **Differential Abundance**: Identify taxa that differ between clinical groups
4. **Longitudinal Analysis**: Track changes in microbiome composition over time
5. **Reporting**: Generate summary visualizations and reports

### Using metaphlan_tools

The metaphlan_tools package is integrated as a submodule in the tools/ directory. You can use it in your analysis scripts and notebooks:

```python
# Example script usage
import sys
import os

# Add the tools directory to Python path
tools_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'tools')
sys.path.append(tools_dir)

# Import functions from metaphlan_tools
from metaphlan_tools import (
    parse_metaphlan_file, 
    combine_samples, 
    load_metadata,
    calculate_alpha_diversity,
    plot_relative_abundance_heatmap
)

# Your analysis code here
```

### Configuration Management

Store analysis parameters in configuration files to:
- Make analyses reproducible
- Allow easy adjustment of parameters
- Document analysis settings

### Documentation

Maintain comprehensive documentation:
- README.md for project overview
- Data dictionary explaining metadata variables
- Workflow documentation describing analysis steps
- Submodule workflow for managing the metaphlan_tools integration

## Getting Started

1. Clone the repository and initialize the submodule:
   ```bash
   git clone --recurse-submodules https://github.com/yourusername/RSV-microbiome-2025.git
   ```

2. Create and activate the conda environment:
   ```bash
   conda env create -f environment.yml
   conda activate rsv-microbiome
   ```

3. Install the metaphlan_tools package in development mode:
   ```bash
   pip install -e tools/metaphlan_tools
   ```

4. Begin your analysis with the notebooks or scripts
