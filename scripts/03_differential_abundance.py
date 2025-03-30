#!/usr/bin/env python3
"""
Identify differentially abundant microbial taxa between clinical groups.

This script:
1. Loads the combined abundance table and metadata
2. Performs statistical testing for each species across clinical groups
3. Applies multiple testing correction to control false discovery rate
4. Identifies significantly differentially abundant species
5. Creates visualizations of key differential species
6. Ranks species by statistical significance and effect size

Usage:
    python scripts/03_identify_differential_taxa.py [--config CONFIG_FILE]
"""

import os
import sys
import argparse
import yaml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Add project root to Python path
project_root = Path(__file__).resolve().parents[1]
sys.path.append(str(project_root))

# Add tools directory to Python path
tools_dir = project_root / 'tools'
sys.path.append(str(tools_dir))


# Import functions from metaphlan_tools
from metaphlan_tools import (
    load_metadata,
    differential_abundance_analysis,
    plot_stacked_bar
)

# Import additional visualization function
from metaphlan_tools.stats import plot_abundance_boxplot


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Identify differentially abundant taxa')
    parser.add_argument('--config', type=str, default='config/analysis_parameters.yml',
                       help='Path to configuration file')
    return parser.parse_args()


def filter_species(abundance_df, min_prevalence=0.1, min_abundance=0.01):
    """
    Filter species based on prevalence and abundance thresholds.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Species abundance DataFrame with species as index, samples as columns
    min_prevalence : float
        Minimum fraction of samples in which a species must be present
    min_abundance : float
        Minimum mean relative abundance a species must have
        
    Returns:
    --------
    pandas.DataFrame
        Filtered abundance DataFrame
    """
    # Calculate prevalence (fraction of samples where species is present)
    prevalence = (abundance_df > 0).mean(axis=1)
    
    # Calculate mean abundance
    mean_abundance = abundance_df.mean(axis=1)
    
    # Filter based on thresholds
    keep_species = (prevalence >= min_prevalence) & (mean_abundance >= min_abundance)
    
    print(f"Filtering from {len(abundance_df)} to {keep_species.sum()} species")
    print(f"  Prevalence threshold: {min_prevalence:.2f} (must be present in {min_prevalence*100:.1f}% of samples)")
    print(f"  Abundance threshold: {min_abundance:.4f} (must have mean abundance ≥ {min_abundance*100:.2f}%)")
    
    return abundance_df.loc[keep_species]


def main():
    """Main function to identify differential taxa."""
    # Parse arguments
    args = parse_args()
    
    # Load configuration
    config_path = project_root / args.config
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    # Set up paths
    processed_data_dir = project_root / 'data' / 'processed'
    results_dir = project_root / 'results'
    figures_dir = results_dir / 'figures'
    tables_dir = results_dir / 'tables'
    species_boxplots_dir = figures_dir / 'species_boxplots'
    
    # Create directories if they don't exist
    figures_dir.mkdir(exist_ok=True, parents=True)
    tables_dir.mkdir(exist_ok=True, parents=True)
    species_boxplots_dir.mkdir(exist_ok=True, parents=True)
    
    # Load abundance data
    abundance_file = processed_data_dir / 'combined_abundance.csv'
    if not abundance_file.exists():
        print(f"Error: Abundance file not found at {abundance_file}")
        print("Please run 01_process_metaphlan_files.py first.")
        sys.exit(1)
    
    print(f"Loading abundance data from {abundance_file}")
    abundance_df = pd.read_csv(abundance_file, index_col=0)
    
    # Load metadata
    metadata_file = project_root / config['metadata']['filename']
    if not metadata_file.exists():
        print(f"Error: Metadata file not found at {metadata_file}")
        sys.exit(1)
    
    print(f"Loading metadata from {metadata_file}")
    metadata_df = load_metadata(metadata_file, config['metadata']['sample_id_column'])
    
    # Check data dimensions
    print(f"Abundance data: {abundance_df.shape[0]} species, {abundance_df.shape[1]} samples")
    print(f"Metadata: {metadata_df.shape[0]} samples, {metadata_df.shape[1]} variables")
    
    # Check sample overlap
    common_samples = set(abundance_df.columns).intersection(set(metadata_df.index))
    print(f"Samples with both abundance and metadata: {len(common_samples)}")
    
    # Filter to common samples if needed
    if len(common_samples) < len(abundance_df.columns):
        print(f"Filtering to {len(common_samples)} common samples")
        abundance_df = abundance_df[list(common_samples)]
    
    # Filter species by prevalence and abundance
    min_prevalence = config['differential_abundance']['min_prevalence']
    min_abundance = config['differential_abundance']['min_abundance']
    filtered_abundance_df = filter_species(abundance_df, min_prevalence, min_abundance)
    
    # Save filtered abundance table
    filtered_file = processed_data_dir / 'filtered_abundance.csv'
    filtered_abundance_df.to_csv(filtered_file)
    print(f"Filtered abundance table saved to {filtered_file}")
    
    # Perform differential abundance analysis for each group variable
    group_vars = config['metadata']['group_variables']
    
    for var in group_vars:
        if var in metadata_df.columns:
            print(f"\nPerforming differential abundance analysis by {var}")
            
            # Check if variable has enough groups
            n_groups = metadata_df[var].nunique()
            if n_groups < 2:
                print(f"  Skipping {var}: needs at least 2 groups, found {n_groups}")
                continue
            
            # Perform differential abundance analysis
            results = differential_abundance_analysis(
                filtered_abundance_df, 
                metadata_df, 
                var,
                method='wilcoxon' if n_groups == 2 else 'kruskal'
            )
            
            # Save results
            da_file = tables_dir / f"differential_abundance_{var}.csv"
            results.to_csv(da_file)
            print(f"  Results saved to {da_file}")
            
            # Print top significant results
            alpha = config['differential_abundance']['p_value_threshold']
            sig_species = results[results['Adjusted P-value'] < alpha]
            print(f"  Found {len(sig_species)} significantly different species (adj. p < {alpha})")
            
            if not sig_species.empty:
                top_n = min(5, len(sig_species))
                print(f"\n  Top {top_n} significant species:")
                for i, (_, row) in enumerate(sig_species.head(top_n).iterrows()):
                    print(f"    {i+1}. {row['Species']}: adj p-value = {row['Adjusted P-value']:.4e}, " +
                          f"fold change = {row['Fold Change']:.2f}")
                
                # Create visualizations for top significant species
                print(f"\n  Creating visualizations for top {top_n} species")
                for i, (_, row) in enumerate(sig_species.head(top_n).iterrows()):
                    species = row['Species']
                    
                    # Create boxplot
                    try:
                        fig = plot_abundance_boxplot(filtered_abundance_df, metadata_df, species, var)
                        # Add p-value annotation
                        plt.annotate(f"Adj. p = {row['Adjusted P-value']:.2e}", 
                                    xy=(0.5, 0.01), xycoords='axes fraction',
                                    ha='center', va='bottom', fontsize=10,
                                    bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=0.7))
                        
                        fig_file = species_boxplots_dir / f"{species.replace(' ', '_')}_{var}.png"
                        fig.savefig(fig_file, dpi=config['visualization']['figure_dpi'], bbox_inches='tight')
                        plt.close(fig)
                        print(f"    Boxplot saved to {fig_file}")
                    except Exception as e:
                        print(f"    Error creating boxplot for {species}: {str(e)}")
                
                # Create a stacked bar plot for species composition
                try:
                    top_species = sig_species.head(10)['Species'].tolist()
                    species_subset = filtered_abundance_df.loc[filtered_abundance_df.index.isin(top_species)]
                    
                    fig = plot_stacked_bar(species_subset, metadata_df, var)
                    fig_file = figures_dir / f"composition_by_{var}.png"
                    fig.savefig(fig_file, dpi=config['visualization']['figure_dpi'], bbox_inches='tight')
                    plt.close(fig)
                    print(f"  Stacked bar plot saved to {fig_file}")
                except Exception as e:
                    print(f"  Error creating stacked bar plot: {str(e)}")
            else:
                print("  No significantly different species found")
        else:
            print(f"Warning: Variable '{var}' not found in metadata")
    
    # Create a heatmap of top species
    print("\nCreating heatmap of top variable species...")
    try:
        # Get list of all significant species across all variables
        all_sig_species = []
        
        for var in group_vars:
            if var in metadata_df.columns:
                da_file = tables_dir / f"differential_abundance_{var}.csv"
                if da_file.exists():
                    df = pd.read_csv(da_file)
                    if 'Adjusted P-value' in df.columns and 'Species' in df.columns:
                        sig = df[df['Adjusted P-value'] < alpha]['Species'].tolist()
                        all_sig_species.extend(sig)
        
        # Get unique species
        all_sig_species = list(set(all_sig_species))
        
        if all_sig_species:
            # Limit to top 30 species if there are too many
            if len(all_sig_species) > 30:
                print(f"  Limiting heatmap to top 30 species out of {len(all_sig_species)} significant species")
                
                # Get mean abundance for each significant species
                mean_abundance = filtered_abundance_df.loc[
                    filtered_abundance_df.index.isin(all_sig_species)
                ].mean(axis=1)
                
                # Sort by abundance and take top 30
                all_sig_species = mean_abundance.sort_values(ascending=False).head(30).index.tolist()
            
            # Import specialized heatmap function
            from metaphlan_tools import plot_relative_abundance_heatmap
            
            # Create heatmap for significant species
            species_subset = filtered_abundance_df.loc[filtered_abundance_df.index.isin(all_sig_species)]
            
            fig = plot_relative_abundance_heatmap(
                species_subset, 
                metadata_df, 
                group_vars[0] if len(group_vars) > 0 else None,
                cluster_samples=True,
                cluster_taxa=True,
                cmap=config['visualization']['heatmap_colormap'],
                figsize=(12, 10)
            )
            
            heatmap_file = figures_dir / "differential_species_heatmap.png"
            fig.savefig(heatmap_file, dpi=config['visualization']['figure_dpi'], bbox_inches='tight')
            plt.close(fig)
            print(f"Heatmap saved to {heatmap_file}")
        else:
            print("  No significant species found for heatmap")
    except Exception as e:
        print(f"Error creating heatmap: {str(e)}")
    
    print("\nDifferential abundance analysis complete!")


if __name__ == "__main__":
    main()
