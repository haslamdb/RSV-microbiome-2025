# RSV Microbiome Analysis: Command-line Workflow

This document outlines a step-by-step command-line workflow for analyzing nasal microbiome data in the RSV-microbiome-2025 project using the metaphlan_tools package.

## Prerequisites

1. The metaphlan_tools package is installed or added as a submodule
2. Raw MetaPhlAn output files are in the `data/raw` directory
3. The metadata.csv file is in the project root directory

## Setup

First, create the necessary directories if they don't exist:

```bash
# Create directory structure
mkdir -p data/raw data/processed results/figures results/tables
```

## 1. Reading Metadata and Extracting Sample IDs

First, let's extract sample IDs from the metadata file:

```bash
# Extract sample IDs to a file for reference
python -c "import pandas as pd; pd.read_csv('metadata.csv')['SampleID'].to_csv('data/processed/sample_ids.txt', index=False, header=False)"

# View the sample IDs
cat data/processed/sample_ids.txt

# Check how many samples we have
wc -l data/processed/sample_ids.txt
```

## 2. Processing MetaPhlAn Files into a Combined Abundance Table

Use the `process` command from metaphlan_tools CLI to combine all MetaPhlAn output files:

```bash
# Navigate to the project root
cd /path/to/RSV-microbiome-2025

# Process all MetaPhlAn files in the raw data directory
metaphlan_tools process --input-dir data/raw --output-dir data/processed --file-pattern "*.txt"

# Check the combined abundance table
head -n 5 data/processed/combined_abundance.csv
```

This creates a combined abundance table with species as rows and samples as columns.

## 3. Analyzing Alpha and Beta Diversity

Next, analyze the diversity metrics based on clinical variables (Timing and Symptoms):

```bash
# Calculate and analyze diversity metrics by Timing
metaphlan_tools diversity --metadata-file metadata.csv --output-dir results --group-var Timing

# Calculate and analyze diversity metrics by Symptoms
metaphlan_tools diversity --metadata-file metadata.csv --output-dir results --group-var Symptoms
```

This will:
- Calculate alpha diversity metrics (Shannon, Simpson, observed species)
- Compare alpha diversity between groups using statistical tests
- Calculate beta diversity distance matrix
- Perform PERMANOVA to test for community composition differences
- Generate visualization of alpha diversity comparisons

## 4. Differential Abundance Testing

Identify species that differ significantly based on Timing or Symptoms:

```bash
# Differential abundance analysis by Timing
metaphlan_tools differential --metadata-file metadata.csv --output-dir results --group-var Timing

# Differential abundance analysis by Symptoms
metaphlan_tools differential --metadata-file metadata.csv --output-dir results --group-var Symptoms
```

This will:
- Perform statistical testing for each species
- Apply multiple testing correction
- Generate a table of differentially abundant species with statistics
- Create visualizations for top significant species

## 5. Longitudinal Analysis

Analyze changes in microbiome over time points (Prior, Acute, Post):

```bash
# Longitudinal analysis with Symptoms as additional grouping variable
metaphlan_tools longitudinal --metadata-file metadata.csv --output-dir results \
  --time-var Timing --subject-var SubjectID --group-var Symptoms
```

This will:
- Track changes in top species across time points
- Group subjects by symptoms severity
- Generate longitudinal plots showing trends

## 6. Generate Summary Report

Create a summary report with key visualizations:

```bash
# Generate comprehensive report
metaphlan_tools report --metadata-file metadata.csv --output-dir results --group-var Symptoms
```

## 7. Custom Analysis: Box Plots for Significant Species

To create box plots for specific species that differ significantly based on Timing or Symptoms:

```bash
# First, identify top significant species from differential abundance results
TOP_SPECIES=$(head -n 6 results/differential_abundance_Timing.csv | tail -n 5 | cut -d',' -f1)

# Create a script to generate box plots for these species
cat > scripts/plot_significant_species.py << 'EOF'
#!/usr/bin/env python3
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from metaphlan_tools import load_metadata
from metaphlan_tools.stats import plot_abundance_boxplot

# Load data
abundance_df = pd.read_csv('data/processed/combined_abundance.csv', index_col=0)
metadata_df = load_metadata('metadata.csv')

# Get species from command line arguments
species_list = sys.argv[1:]
variables = ['Timing', 'Symptoms']

# Create output directory
os.makedirs('results/figures/species_boxplots', exist_ok=True)

# Generate boxplots for each species and variable
for species in species_list:
    for var in variables:
        print(f"Plotting {species} by {var}")
        try:
            fig = plot_abundance_boxplot(abundance_df, metadata_df, species, var)
            output_file = f'results/figures/species_boxplots/{species.replace(" ", "_")}_{var}.png'
            fig.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close(fig)
            print(f"Saved to {output_file}")
        except Exception as e:
            print(f"Error plotting {species}: {str(e)}")
EOF

# Make the script executable
chmod +x scripts/plot_significant_species.py

# Run the script with the top species
python scripts/plot_significant_species.py $TOP_SPECIES
```

## 8. One-line Complete Analysis Pipeline

For convenience, here's a one-line sequence to run the entire analysis:

```bash
# Complete analysis pipeline
mkdir -p data/raw data/processed results/figures results/tables && \
metaphlan_tools process --input-dir data/raw --output-dir data/processed && \
metaphlan_tools diversity --metadata-file metadata.csv --output-dir results --group-var Timing && \
metaphlan_tools diversity --metadata-file metadata.csv --output-dir results --group-var Symptoms && \
metaphlan_tools differential --metadata-file metadata.csv --output-dir results --group-var Timing && \
metaphlan_tools differential --metadata-file metadata.csv --output-dir results --group-var Symptoms && \
metaphlan_tools longitudinal --metadata-file metadata.csv --output-dir results --time-var Timing --subject-var SubjectID --group-var Symptoms && \
metaphlan_tools report --metadata-file metadata.csv --output-dir results --group-var Symptoms
```

## Output File Structure

After running this workflow, your results directory should contain:

```
results/
├── alpha_diversity.csv
├── alpha_diversity_Timing_boxplot.png
├── alpha_diversity_Timing_stats.csv
├── alpha_diversity_Symptoms_boxplot.png
├── alpha_diversity_Symptoms_stats.csv
├── permanova_Timing_results.txt
├── permanova_Symptoms_results.txt
├── differential_abundance_Timing.csv
├── differential_abundance_Symptoms.csv
├── figures/
│   ├── abundance_heatmap.png
│   ├── abundance_barplot.png
│   ├── species_boxplots/
│   │   ├── Species_Name_Timing.png
│   │   └── Species_Name_Symptoms.png
│   └── longitudinal_plots/
└── tables/
```

## Notes

- This workflow assumes the structure of your metaphlan_tools CLI as shown in the files you shared
- Adjust file paths as needed based on your actual directory structure
- If any command fails, check the error message and make sure the input files exist
- The metaphlan_tools commands may need additional parameters depending on your specific data characteristics
