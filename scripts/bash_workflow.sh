#!/bin/bash
# rsv_microbiome_analysis.sh
#
# Complete workflow for analyzing RSV microbiome data using metaphlan_tools CLI
#
# Usage: ./rsv_microbiome_analysis.sh

# Stop on errors
set -e

# Print commands as they are executed
set -x

# Set variables
PROJECT_DIR=$(pwd)
RAW_DATA_DIR="${PROJECT_DIR}/data/raw"
PROCESSED_DATA_DIR="${PROJECT_DIR}/data/processed"
RESULTS_DIR="${PROJECT_DIR}/results"
FIGURES_DIR="${RESULTS_DIR}/figures"
TABLES_DIR="${RESULTS_DIR}/tables"
METADATA_FILE="${PROJECT_DIR}/metadata.csv"

# Create directory structure
echo "Creating directory structure..."
mkdir -p "${RAW_DATA_DIR}" "${PROCESSED_DATA_DIR}" "${FIGURES_DIR}" "${TABLES_DIR}"
mkdir -p "${FIGURES_DIR}/species_boxplots" "${FIGURES_DIR}/longitudinal_plots"

# Step 1: Extract sample IDs from metadata
echo "Step 1: Extracting sample IDs from metadata..."
python -c "import pandas as pd; pd.read_csv('${METADATA_FILE}')['SampleID'].to_csv('${PROCESSED_DATA_DIR}/sample_ids.txt', index=False, header=False)"
echo "Found $(wc -l < ${PROCESSED_DATA_DIR}/sample_ids.txt) samples in metadata"

# Step 2: Process MetaPhlAn files into a combined abundance table
echo "Step 2: Processing MetaPhlAn files..."
metaphlan_tools process \
  --input-dir "${RAW_DATA_DIR}" \
  --output-dir "${PROCESSED_DATA_DIR}" \
  --file-pattern "*.txt"

# Check if combined abundance table was created
if [ ! -f "${PROCESSED_DATA_DIR}/combined_abundance.csv" ]; then
  echo "Error: Combined abundance table not created!"
  exit 1
fi
echo "Combined abundance table created successfully"

# Step 3: Analyze alpha and beta diversity
echo "Step 3: Analyzing diversity metrics..."

# By Timing
echo "Analyzing diversity by Timing..."
metaphlan_tools diversity \
  --metadata-file "${METADATA_FILE}" \
  --output-dir "${RESULTS_DIR}" \
  --group-var "Timing"

# By Symptoms
echo "Analyzing diversity by Symptoms..."
metaphlan_tools diversity \
  --metadata-file "${METADATA_FILE}" \
  --output-dir "${RESULTS_DIR}" \
  --group-var "Symptoms"

# Step 4: Differential abundance testing
echo "Step 4: Performing differential abundance testing..."

# By Timing
echo "Testing differential abundance by Timing..."
metaphlan_tools differential \
  --metadata-file "${METADATA_FILE}" \
  --output-dir "${RESULTS_DIR}" \
  --group-var "Timing"

# By Symptoms
echo "Testing differential abundance by Symptoms..."
metaphlan_tools differential \
  --metadata-file "${METADATA_FILE}" \
  --output-dir "${RESULTS_DIR}" \
  --group-var "Symptoms"

# Step 5: Longitudinal analysis
echo "Step 5: Performing longitudinal analysis..."
metaphlan_tools longitudinal \
  --metadata-file "${METADATA_FILE}" \
  --output-dir "${RESULTS_DIR}" \
  --time-var "Timing" \
  --subject-var "SubjectID" \
  --group-var "Symptoms"

# Step 6: Generate summary report
echo "Step 6: Generating summary report..."
metaphlan_tools report \
  --metadata-file "${METADATA_FILE}" \
  --output-dir "${RESULTS_DIR}" \
  --group-var "Symptoms"

# Step 7: Create custom box plots for significant species
echo "Step 7: Creating box plots for significant species..."

# Create a temporary Python script for box plots
cat > "${PROJECT_DIR}/scripts/plot_significant_species.py" << 'EOF'
#!/usr/bin/env python3
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from metaphlan_tools import load_metadata
from metaphlan_tools.stats import plot_abundance_boxplot
from scipy import stats

# Set up directories
project_dir = os.getcwd()
processed_dir = os.path.join(project_dir, 'data', 'processed')
results_dir = os.path.join(project_dir, 'results')
figures_dir = os.path.join(results_dir, 'figures', 'species_boxplots')
tables_dir = os.path.join(results_dir, 'tables')

# Create output directories if they don't exist
os.makedirs(figures_dir, exist_ok=True)
os.makedirs(tables_dir, exist_ok=True)

# Load data
abundance_file = os.path.join(processed_dir, 'combined_abundance.csv')
metadata_file = os.path.join(project_dir, 'metadata.csv')

print(f"Loading abundance data from {abundance_file}")
abundance_df = pd.read_csv(abundance_file, index_col=0)
print(f"Loading metadata from {metadata_file}")
metadata_df = load_metadata(metadata_file)

# Get top significant species from differential abundance results
def get_top_species(diff_file, n=5):
    if not os.path.exists(diff_file):
        print(f"Warning: {diff_file} not found")
        return []
    diff_df = pd.read_csv(diff_file)
    if diff_df.empty:
        return []
    # Ensure 'Species' and 'Adjusted P-value' columns exist
    if 'Species' not in diff_df.columns:
        species_col = diff_df.columns[0]  # Assume first column is species
    else:
        species_col = 'Species'
    
    if 'Adjusted P-value' not in diff_df.columns:
        # Try to find a p-value column
        pval_cols = [col for col in diff_df.columns if 'p-value' in col.lower()]
        if pval_cols:
            pval_col = pval_cols[0]
        else:
            # If no p-value column, sort by fold change or just take first few
            print("No p-value column found, using first few species")
            return diff_df[species_col].head(n).tolist()
    else:
        pval_col = 'Adjusted P-value'
    
    # Filter significant species and take top n
    sig_species = diff_df[diff_df[pval_col] < 0.05][species_col].tolist()
    return sig_species[:min(n, len(sig_species))]

# Get top species from each differential analysis
timing_file = os.path.join(results_dir, 'differential_abundance_Timing.csv')
symptoms_file = os.path.join(results_dir, 'differential_abundance_Symptoms.csv')

timing_species = get_top_species(timing_file)
symptoms_species = get_top_species(symptoms_file)

# Combine unique species
species_list = list(set(timing_species + symptoms_species))

print(f"Analyzing {len(species_list)} significant species")

# Variables to analyze
variables = ['Timing', 'Symptoms']

# Create a summary table for statistical results
stats_results = []

# Generate boxplots for each species and variable
for species in species_list:
    print(f"\nAnalyzing {species}:")
    for var in variables:
        print(f"  Testing by {var}")
        
        # Get data for this species and join with metadata
        if species not in abundance_df.index:
            print(f"  Warning: Species '{species}' not found in abundance data")
            continue
            
        species_data = abundance_df.loc[species]
        plot_df = pd.DataFrame({'Abundance': species_data})
        plot_df = plot_df.join(metadata_df[[var]])
        
        # Calculate statistics
        groups = plot_df.groupby(var)['Abundance'].apply(list).values()
        if len(groups) < 2:
            print(f"  Warning: Not enough groups for statistical testing")
            continue
            
        # Choose statistical test based on number of groups
        if len(groups) == 2:
            # Use Mann-Whitney U test for two groups
            try:
                stat, pval = stats.mannwhitneyu(*groups)
                test_name = 'Mann-Whitney U'
            except ValueError:
                # Fall back to Kruskal-Wallis if mann-whitney fails
                stat, pval = stats.kruskal(*groups)
                test_name = 'Kruskal-Wallis'
        else:
            # Use Kruskal-Wallis for more than two groups
            stat, pval = stats.kruskal(*groups)
            test_name = 'Kruskal-Wallis'
            
        # Add to results
        stats_results.append({
            'Species': species,
            'Variable': var,
            'Test': test_name,
            'Statistic': stat,
            'P-value': pval,
            'Significant': pval < 0.05
        })
        
        print(f"  {test_name}: p-value = {pval:.4f}")
        
        # Create and save boxplot
        try:
            fig = plot_abundance_boxplot(abundance_df, metadata_df, species, var)
            
            # Add p-value annotation to plot
            plt.annotate(f"p = {pval:.4f}", xy=(0.5, 0.01), xycoords='axes fraction',
                         ha='center', va='bottom', fontsize=10,
                         bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=0.7))
            
            # Save plot
            output_file = os.path.join(figures_dir, f"{species.replace(' ', '_')}_{var}.png")
            fig.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close(fig)
            print(f"  Saved boxplot to {output_file}")
        except Exception as e:
            print(f"  Error plotting {species}: {str(e)}")

# Save statistical results
stats_df = pd.DataFrame(stats_results)
stats_file = os.path.join(tables_dir, 'species_statistical_tests.csv')
stats_df.to_csv(stats_file, index=False)
print(f"\nSaved statistical test results to {stats_file}")

print("\nAnalysis complete!")
EOF

# Make script executable
chmod +x "${PROJECT_DIR}/scripts/plot_significant_species.py"

# Run the script
python "${PROJECT_DIR}/scripts/plot_significant_species.py"

echo "Analysis complete! Results are in ${RESULTS_DIR}"
echo "Figures are in ${FIGURES_DIR}"
echo "Tables are in ${TABLES_DIR}"

# List key output files
echo "Key output files:"
find "${RESULTS_DIR}" -type f -name "*.csv" -o -name "*.png" | sort
