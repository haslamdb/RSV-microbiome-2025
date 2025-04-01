#!/usr/bin/env python3
"""
Analyze RSV nasal microbiome data from Kraken2/Bracken across timing points based on severity and symptoms.

This script:
1. Loads processed Kraken2/Bracken abundance data
2. Performs differential abundance testing:
   - Between timing points (Prior, Acute, Post)
   - Between severity groups at each time point
   - Between symptom groups at each time point
3. Generates faceted boxplots for differentially abundant species
4. Uses kraken_tools for enhanced analysis capabilities

Usage:
    python scripts/kraken/02_kraken_differential_abundance.py [--abundance-file ABUNDANCE_FILE] [--output-dir OUTPUT_DIR]
"""

import os
import sys
import glob
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import re
from scipy.stats import mannwhitneyu, kruskal
from statsmodels.stats.multitest import multipletests
import traceback

# Add project root to Python path
project_root = Path(__file__).resolve().parents[2]
sys.path.append(str(project_root))

# Add kraken_tools to Python path
kraken_tools_dir = Path.home() / "Documents" / "Code" / "kraken_tools"
sys.path.append(str(kraken_tools_dir))

# Import functions from kraken_tools
try:
    from kraken_tools.analysis.differential import run_differential_abundance_analysis
    from kraken_tools.analysis.metadata import read_and_process_metadata
    from kraken_tools.analysis.normalization_and_filtering import filter_low_abundance_taxa
    kraken_tools_available = True
except ImportError:
    kraken_tools_available = False
    print("Warning: kraken_tools module not available. Using built-in functions.")

# Set plotting style
plt.style.use('seaborn-v0_8-whitegrid')
sns.set(font_scale=1.1)


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Analyze RSV microbiome data from Kraken2/Bracken')
    parser.add_argument('--abundance-file', type=str, 
                        default='results/kraken_analysis/filtered_kraken_s_abundance.tsv',
                        help='Path to processed Kraken abundance file')
    parser.add_argument('--output-dir', type=str, default='results/kraken_differential',
                       help='Directory to save analysis results')
    parser.add_argument('--metadata', type=str, default='metadata.csv',
                       help='Path to metadata file')
    parser.add_argument('--sample-id-column', type=str, default='SampleID',
                       help='Column name for sample IDs in metadata')
    parser.add_argument('--taxonomic-level', type=str, default='S',
                       choices=['D', 'P', 'C', 'O', 'F', 'G', 'S'],
                       help='Taxonomic level for analysis')
    parser.add_argument('--group-by', type=str, default='Severity,Symptoms,Timing',
                       help='Comma-separated list of grouping variables for differential testing')
    parser.add_argument('--min-prevalence', type=float, default=0.1,
                       help='Minimum prevalence to keep taxon')
    parser.add_argument('--min-abundance', type=float, default=0.01,
                       help='Minimum relative abundance to keep taxon')
    parser.add_argument('--methods', type=str, default='aldex2,ancom,ancom-bc',
                       help='Comma-separated list of differential abundance methods to use')
    parser.add_argument('--p-value-threshold', type=float, default=0.05,
                       help='P-value threshold for significance')
    parser.add_argument('--top-taxa', type=int, default=10,
                       help='Number of top differential taxa to visualize')
    parser.add_argument('--stratify-by-timing', action='store_true',
                       help='Stratify analyses by timing point (Prior, Acute, Post)')
    return parser.parse_args()


# We need a function to process Kraken taxonomic names
def parse_kraken_taxonomy(taxon_name):
    """
    Parse Kraken/Bracken taxonomic name into components.
    
    Parameters:
    -----------
    taxon_name : str
        Taxonomic name from Kraken/Bracken
        
    Returns:
    --------
    dict
        Dictionary with taxonomic information
    """
    # Check if we have a standard format
    if not taxon_name or not isinstance(taxon_name, str):
        return {'name': str(taxon_name), 'genus': 'Unknown', 'species': 'Unknown'}
    
    # Bracken usually outputs names like 'Escherichia coli' or 'Homo sapiens'
    # or with full taxonomic path like 'd__Bacteria|p__Proteobacteria|...|s__Escherichia coli'
    
    # For full path format
    if '|' in taxon_name:
        components = taxon_name.split('|')
        # Get the last component which should be species
        species_part = components[-1]
        if species_part.startswith('s__'):
            species = species_part[3:]  # Remove the s__ prefix
        else:
            species = species_part
            
        # Try to get genus
        genus = "Unknown"
        for comp in components:
            if comp.startswith('g__'):
                genus = comp[3:]  # Remove the g__ prefix
                break
        
        # Return parsed information
        return {
            'name': taxon_name,
            'genus': genus,
            'species': species if species != '' else 'Unknown'
        }
    
    # For standard binomial format
    parts = taxon_name.split()
    if len(parts) >= 2:
        genus = parts[0]
        species = ' '.join(parts[1:])
        return {'name': taxon_name, 'genus': genus, 'species': species}
    
    # For genus-only or other formats
    return {'name': taxon_name, 'genus': taxon_name, 'species': 'Unknown'}


def load_abundance_data(file_path):
    """
    Load Kraken/Bracken abundance data from file.
    
    Parameters:
    -----------
    file_path : str
        Path to abundance file (can be TSV or CSV)
        
    Returns:
    --------
    pandas.DataFrame
        Abundance DataFrame with taxa as index, samples as columns
    """
    try:
        # Determine file format based on extension
        if file_path.lower().endswith('.csv'):
            df = pd.read_csv(file_path, index_col=0)
        else:
            # Default to tab-separated
            df = pd.read_csv(file_path, sep='\t', index_col=0)
        
        # Quick validation
        if df.shape[0] == 0 or df.shape[1] == 0:
            raise ValueError(f"Abundance file has {df.shape[0]} rows and {df.shape[1]} columns")
        
        # Fill NAs with zeros
        df = df.fillna(0)
        
        return df
        
    except Exception as e:
        print(f"Error loading abundance file: {str(e)}")
        traceback.print_exc()
        sys.exit(1)
        
        
def run_differential_analysis(abundance_df, metadata_df, group_var, output_dir, methods=None, min_abundance=0.01, min_prevalence=0.1, p_value_threshold=0.05):
    """
    Run differential abundance analysis using kraken_tools if available, or built-in functions.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Abundance DataFrame with taxa as index, samples as columns
    metadata_df : pandas.DataFrame
        Metadata DataFrame with samples as index
    group_var : str
        Metadata variable for grouping
    output_dir : str
        Directory to save results
    methods : list
        List of methods to use (aldex2, ancom, ancom-bc)
    min_abundance : float
        Minimum mean abundance threshold
    min_prevalence : float
        Minimum prevalence threshold  
    p_value_threshold : float
        Threshold for significance
        
    Returns:
    --------
    dict
        Results from the differential abundance analysis
    """
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # If kraken_tools available, use it
    if kraken_tools_available:
        try:
            print(f"Running differential abundance analysis with kraken_tools for {group_var}")
            
            # Filter abundance data if needed
            if min_abundance > 0 or min_prevalence > 0:
                filtered_df = filter_low_abundance_taxa(
                    abundance_df, 
                    min_abundance=min_abundance, 
                    min_prevalence=min_prevalence
                )
            else:
                filtered_df = abundance_df.copy()
                
            # Run differential abundance analysis
            results = run_differential_abundance_analysis(
                filtered_df,
                metadata_df,
                output_dir,
                group_col=group_var,
                methods=methods
            )
            
            return results
            
        except Exception as e:
            print(f"Error using kraken_tools for differential analysis: {str(e)}")
            print("Falling back to built-in methods")
    
    # Fall back to built-in methods
    print(f"Running built-in differential abundance analysis for {group_var}")
    
    # Get values for grouping
    if group_var not in metadata_df.columns:
        print(f"Error: Group variable '{group_var}' not found in metadata")
        return None
        
    # Find common samples
    common_samples = list(set(abundance_df.columns).intersection(set(metadata_df.index)))
    if len(common_samples) < 5:
        print(f"Error: Too few samples ({len(common_samples)}) for differential analysis")
        return None
        
    # Filter to common samples
    filtered_abundance = abundance_df[common_samples]
    
    # Get group information
    groups = metadata_df.loc[common_samples, group_var]
    unique_groups = groups.unique()
    
    # Filter by abundance and prevalence
    if min_abundance > 0 or min_prevalence > 0:
        # Calculate mean abundance per taxon
        mean_abundance = filtered_abundance.mean(axis=1)
        
        # Calculate prevalence (proportion of samples where taxon is present)
        prevalence = (filtered_abundance > 0).mean(axis=1)
        
        # Apply filters
        before_count = filtered_abundance.shape[0]
        filtered_abundance = filtered_abundance.loc[
            (mean_abundance >= min_abundance) & (prevalence >= min_prevalence)
        ]
        after_count = filtered_abundance.shape[0]
        
        print(f"Filtered from {before_count} to {after_count} taxa based on abundance/prevalence thresholds")
    
    # Decide which test to use based on number of groups
    if len(unique_groups) == 2:
        print(f"Running Mann-Whitney U test for {group_var}")
        results = run_mann_whitney_test(filtered_abundance, groups, unique_groups, p_value_threshold)
    else:
        print(f"Running Kruskal-Wallis test for {group_var} with {len(unique_groups)} groups")
        results = run_kruskal_wallis_test(filtered_abundance, groups, unique_groups, p_value_threshold)
    
    # Save results
    results_file = os.path.join(output_dir, f"diff_abundance_{group_var}.tsv")
    results.to_csv(results_file, sep='\t')
    print(f"Saved differential abundance results to {results_file}")
    
    return results


def run_mann_whitney_test(abundance_df, groups, unique_groups, p_value_threshold=0.05):
    """
    Run Mann-Whitney U test for differential abundance between two groups.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Abundance DataFrame with taxa as index, samples as columns
    groups : pandas.Series
        Group assignment for each sample
    unique_groups : array-like
        Unique group values
    p_value_threshold : float
        Threshold for significance
        
    Returns:
    --------
    pandas.DataFrame
        Results of differential abundance testing
    """
    results = []
    
    # Get sample IDs for each group
    group1_samples = groups[groups == unique_groups[0]].index
    group2_samples = groups[groups == unique_groups[1]].index
    
    # Count number of samples in each group
    n_group1 = len(group1_samples)
    n_group2 = len(group2_samples)
    
    if n_group1 < 3 or n_group2 < 3:
        print(f"Warning: Small sample sizes (group1: {n_group1}, group2: {n_group2})")
    
    # Loop through each taxon
    for taxon in abundance_df.index:
        # Get abundance values for each group
        values1 = abundance_df.loc[taxon, group1_samples]
        values2 = abundance_df.loc[taxon, group2_samples]
        
        # Calculate mean abundance in each group
        mean1 = values1.mean()
        mean2 = values2.mean()
        
        # Calculate fold change (log2)
        # Add small pseudocount to avoid division by zero
        pseudocount = 1e-5
        fold_change = np.log2((mean2 + pseudocount) / (mean1 + pseudocount))
        
        # Perform Mann-Whitney U test
        try:
            stat, p_value = mannwhitneyu(values1, values2, alternative='two-sided')
        except Exception as e:
            # If test fails, set p-value to 1
            p_value = 1.0
        
        # Parse taxonomy (for easier interpretation)
        tax_info = parse_kraken_taxonomy(taxon)
        
        results.append({
            'Taxon': taxon,
            'Genus': tax_info['genus'],
            'Species': tax_info['species'],
            'P-value': p_value,
            'Group1': unique_groups[0],
            'Group2': unique_groups[1],
            'Mean in Group1': mean1,
            'Mean in Group2': mean2,
            'Log2 Fold Change': fold_change,
            'Test': 'Mann-Whitney U'
        })
    
    # Create DataFrame from results
    results_df = pd.DataFrame(results)
    
    # Apply multiple testing correction
    if not results_df.empty and len(results_df) > 1:
        try:
            # Use Benjamini-Hochberg procedure for FDR control
            results_df['Adjusted P-value'] = multipletests(
                results_df['P-value'], method='fdr_bh'
            )[1]
        except Exception as e:
            print(f"Error applying multiple testing correction: {str(e)}")
            results_df['Adjusted P-value'] = results_df['P-value']
    else:
        results_df['Adjusted P-value'] = results_df['P-value']
    
    # Mark significant taxa
    results_df['Significant'] = results_df['Adjusted P-value'] < p_value_threshold
    
    # Sort by adjusted p-value
    results_df = results_df.sort_values('Adjusted P-value')
    
    return results_df


def run_kruskal_wallis_test(abundance_df, groups, unique_groups, p_value_threshold=0.05):
    """
    Run Kruskal-Wallis test for differential abundance between multiple groups.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Abundance DataFrame with taxa as index, samples as columns
    groups : pandas.Series
        Group assignment for each sample
    unique_groups : array-like
        Unique group values
    p_value_threshold : float
        Threshold for significance
        
    Returns:
    --------
    pandas.DataFrame
        Results of differential abundance testing
    """
    results = []
    
    # Loop through each taxon
    for taxon in abundance_df.index:
        # Initialize group values and means
        group_values = {}
        group_means = {}
        
        # Get values and calculate means
        for group in unique_groups:
            group_samples = groups[groups == group].index
            group_values[group] = abundance_df.loc[taxon, group_samples]
            group_means[group] = group_values[group].mean()
        
        # Perform Kruskal-Wallis test
        try:
            # Prepare data for the test
            all_values = [group_values[group] for group in unique_groups]
            stat, p_value = kruskal(*all_values)
        except Exception as e:
            # If test fails, set p-value to 1
            p_value = 1.0
        
        # Calculate maximum fold change between any two groups
        max_fold_change = 0
        max_group_pair = None
        
        for i, group1 in enumerate(unique_groups):
            for group2 in unique_groups[i+1:]:
                # Add small pseudocount to avoid division by zero
                pseudocount = 1e-5
                fold_change = np.abs(np.log2((group_means[group2] + pseudocount) / (group_means[group1] + pseudocount)))
                
                if fold_change > max_fold_change:
                    max_fold_change = fold_change
                    max_group_pair = (group1, group2)
        
        # Parse taxonomy (for easier interpretation)
        tax_info = parse_kraken_taxonomy(taxon)
        
        # Create result entry
        result = {
            'Taxon': taxon,
            'Genus': tax_info['genus'],
            'Species': tax_info['species'],
            'P-value': p_value,
            'Test': 'Kruskal-Wallis'
        }
        
        # Add mean for each group
        for group in unique_groups:
            result[f'Mean in {group}'] = group_means[group]
        
        # Add fold change information
        if max_group_pair:
            result['Max Log2 Fold Change'] = max_fold_change
            result['Max Fold Change Groups'] = f'{max_group_pair[0]} vs {max_group_pair[1]}'
        
        results.append(result)
    
    # Create DataFrame from results
    results_df = pd.DataFrame(results)
    
    # Apply multiple testing correction
    if not results_df.empty and len(results_df) > 1:
        try:
            # Use Benjamini-Hochberg procedure for FDR control
            results_df['Adjusted P-value'] = multipletests(
                results_df['P-value'], method='fdr_bh'
            )[1]
        except Exception as e:
            print(f"Error applying multiple testing correction: {str(e)}")
            results_df['Adjusted P-value'] = results_df['P-value']
    else:
        results_df['Adjusted P-value'] = results_df['P-value']
    
    # Mark significant taxa
    results_df['Significant'] = results_df['Adjusted P-value'] < p_value_threshold
    
    # Sort by adjusted p-value
    results_df = results_df.sort_values('Adjusted P-value')
    
    return results_df


def plot_taxa_boxplots(abundance_df, metadata_df, group_var, results_df, output_dir, top_n=5):
    """
    Create boxplots for top differentially abundant taxa.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Abundance DataFrame with taxa as index, samples as columns
    metadata_df : pandas.DataFrame
        Metadata DataFrame with samples as index
    group_var : str
        Metadata variable for grouping
    results_df : pandas.DataFrame
        Differential abundance results 
    output_dir : str
        Directory to save plots
    top_n : int
        Number of top taxa to visualize
        
    Returns:
    --------
    list
        Paths to created plot files
    """
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Get significant taxa
    sig_taxa = results_df[results_df['Significant'] == True]
    
    if sig_taxa.empty:
        print(f"No significant taxa found for {group_var}")
        return []
    
    # Get top taxa by adjusted p-value
    top_taxa = sig_taxa.head(top_n)['Taxon'].tolist()
    
    # Get common samples
    common_samples = list(set(abundance_df.columns).intersection(set(metadata_df.index)))
    
    # Filter to common samples
    filtered_abundance = abundance_df[common_samples]
    
    # Create a boxplot for each taxon
    plot_files = []
    
    for taxon in top_taxa:
        try:
            # Get taxonomic info for better plot titles
            tax_info = parse_kraken_taxonomy(taxon)
            
            # Create figure
            fig, ax = plt.subplots(figsize=(10, 6))
            
            # Prepare data for plotting
            plot_data = pd.DataFrame({
                'Abundance': filtered_abundance.loc[taxon],
                group_var: metadata_df.loc[filtered_abundance.columns, group_var]
            })
            
            # Create boxplot
            sns.boxplot(x=group_var, y='Abundance', data=plot_data, ax=ax)
            
            # Add individual points
            sns.stripplot(x=group_var, y='Abundance', data=plot_data, 
                         color='black', size=4, alpha=0.5, ax=ax)
            
            # Get p-value info
            p_value = results_df[results_df['Taxon'] == taxon]['Adjusted P-value'].values[0]
            
            # Set titles and labels
            if tax_info['species'] != 'Unknown':
                title = f"{tax_info['genus']} {tax_info['species']}"
            else:
                title = taxon
                
            ax.set_title(f"{title}\n(adj. p-value: {p_value:.2e})")
            ax.set_xlabel(group_var)
            ax.set_ylabel('Relative Abundance')
            
            # Save figure
            safe_taxon = taxon.replace(' ', '_').replace('/', '_').replace(':', '_').replace('|', '_')
            plot_file = os.path.join(output_dir, f"{safe_taxon}_{group_var}.pdf")
            plt.tight_layout()
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            plt.close(fig)
            
            plot_files.append(plot_file)
            
        except Exception as e:
            print(f"Error creating boxplot for {taxon}: {str(e)}")
    
    return plot_files


def main():
    """Main function to run differential abundance analysis."""
    # Parse arguments
    args = parse_args()
    
    # Create output directories
    output_dir = Path(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)
    figures_dir = output_dir / 'figures'
    tables_dir = output_dir / 'tables'
    boxplots_dir = figures_dir / 'boxplots'
    
    for dir_path in [output_dir, figures_dir, tables_dir, boxplots_dir]:
        os.makedirs(dir_path, exist_ok=True)
    
    # Load abundance data
    print(f"Loading abundance data from {args.abundance_file}")
    abundance_df = load_abundance_data(args.abundance_file)
    print(f"Loaded data with {abundance_df.shape[0]} taxa and {abundance_df.shape[1]} samples")
    
    # Load metadata
    print(f"Loading metadata from {args.metadata}")
    if kraken_tools_available:
        metadata_df = read_and_process_metadata(args.metadata)
    else:
        metadata_df = pd.read_csv(args.metadata)
        metadata_df.set_index(args.sample_id_column, inplace=True)
    
    print(f"Loaded metadata with {metadata_df.shape[0]} samples and {metadata_df.shape[1]} columns")
    
    # Get analysis variables
    group_vars = args.group_by.split(',')
    methods = args.methods.split(',')
    
    # Check for timing variable if stratification is requested
    if args.stratify_by_timing:
        if 'Timing' not in metadata_df.columns:
            print("Warning: 'Timing' column not found in metadata. Cannot stratify by timing.")
            args.stratify_by_timing = False
    
    # Run differential abundance analysis for each grouping variable
    for group_var in group_vars:
        if group_var not in metadata_df.columns:
            print(f"Warning: {group_var} not found in metadata. Skipping.")
            continue
            
        print(f"\nAnalyzing differential abundance by {group_var}")
        
        # If stratifying by timing, analyze each timing point separately
        if args.stratify_by_timing and group_var != 'Timing':
            timing_values = metadata_df['Timing'].unique()
            
            for timing in timing_values:
                print(f"\n  Analyzing {group_var} within timing: {timing}")
                
                # Filter metadata to this timing
                timing_metadata = metadata_df[metadata_df['Timing'] == timing]
                
                # Skip if too few samples
                if timing_metadata.shape[0] < 5:
                    print(f"  Too few samples ({timing_metadata.shape[0]}) for timing {timing}. Skipping.")
                    continue
                
                # Create output subdirectory
                timing_dir = output_dir / f"timing_{timing}"
                os.makedirs(timing_dir, exist_ok=True)
                
                # Run differential analysis
                results = run_differential_analysis(
                    abundance_df, 
                    timing_metadata, 
                    group_var, 
                    timing_dir,
                    methods=methods,
                    min_abundance=args.min_abundance,
                    min_prevalence=args.min_prevalence,
                    p_value_threshold=args.p_value_threshold
                )
                
                if results is not None:
                    # Create boxplots
                    boxplots_path = timing_dir / 'boxplots'
                    plot_files = plot_taxa_boxplots(
                        abundance_df,
                        timing_metadata,
                        group_var,
                        results,
                        boxplots_path,
                        top_n=args.top_taxa
                    )
                    
                    # Report significant taxa
                    sig_count = sum(results['Significant'] == True)
                    print(f"  Found {sig_count} significantly different taxa for {group_var} in {timing}")
                
        else:
            # Create output subdirectory
            group_dir = output_dir / f"group_{group_var}"
            os.makedirs(group_dir, exist_ok=True)
            
            # Run differential analysis
            results = run_differential_analysis(
                abundance_df, 
                metadata_df, 
                group_var, 
                group_dir,
                methods=methods,
                min_abundance=args.min_abundance,
                min_prevalence=args.min_prevalence,
                p_value_threshold=args.p_value_threshold
            )
            
            if results is not None:
                # Create boxplots
                boxplots_path = group_dir / 'boxplots'
                plot_files = plot_taxa_boxplots(
                    abundance_df,
                    metadata_df,
                    group_var,
                    results,
                    boxplots_path,
                    top_n=args.top_taxa
                )
                
                # Report significant taxa
                sig_count = sum(results['Significant'] == True)
                print(f"  Found {sig_count} significantly different taxa for {group_var}")
    
    print("\nDifferential abundance analysis complete!")


if __name__ == "__main__":
    main()