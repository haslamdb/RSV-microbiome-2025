#!/usr/bin/env python3
"""
Compare results from Sylph and MetaPhlAn microbiome profiling.

This script:
1. Loads processed abundance tables from both Sylph and MetaPhlAn
2. Standardizes taxonomy names for comparison
3. Analyzes concordance between the two profilers
4. Identifies taxa found by only one method
5. Compares relative abundance distributions
6. Generates visualizations of the comparison

Usage:
    python compare_sylph_metaphlan.py [--config CONFIG_FILE]
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import re
from scipy.stats import spearmanr, pearsonr

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Compare Sylph and MetaPhlAn results')
    parser.add_argument('--sylph-dir', type=str, default='data/processed',
                       help='Directory containing processed Sylph files')
    parser.add_argument('--metaphlan-dir', type=str, default='data/processed',
                       help='Directory containing processed MetaPhlAn files')
    parser.add_argument('--output-dir', type=str, default='results/comparison',
                       help='Directory to save comparison results')
    parser.add_argument('--metadata', type=str, default='metadata.csv',
                       help='Path to metadata file')
    parser.add_argument('--sample-id-column', type=str, default='SampleID',
                       help='Column name for sample IDs in metadata')
    return parser.parse_args()

def load_abundance_table(file_path):
    """
    Load abundance table from CSV file.
    
    Parameters:
    -----------
    file_path : str
        Path to abundance CSV file
        
    Returns:
    --------
    pandas.DataFrame
        Abundance table with taxa as index and samples as columns
    """
    if not os.path.exists(file_path):
        print(f"Warning: File not found: {file_path}")
        return pd.DataFrame()
    
    try:
        df = pd.read_csv(file_path, index_col=0)
        return df
    except Exception as e:
        print(f"Error loading {file_path}: {str(e)}")
        return pd.DataFrame()

def standardize_taxonomy(df, method):
    """
    Standardize taxonomy names for comparison between methods.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        Abundance table with taxa as index and samples as columns
    method : str
        Source of the data ('sylph' or 'metaphlan')
        
    Returns:
    --------
    pandas.DataFrame
        DataFrame with standardized taxonomy index
    """
    # Create a copy to avoid modifying the original
    result_df = df.copy()
    
    # New standardized index
    new_index = []
    
    for taxon in df.index:
        if method == 'sylph':
            # For Sylph, strip any prefixes like "Bacteria: "
            if ':' in taxon:
                taxon = taxon.split(':', 1)[1].strip()
                
            # Convert to lowercase for case-insensitive comparison
            std_taxon = taxon.lower()
        
        elif method == 'metaphlan':
            # For MetaPhlAn, the index might already be the species name
            # Convert to lowercase for case-insensitive comparison
            std_taxon = taxon.lower()
        
        new_index.append(std_taxon)
    
    # Update the index
    result_df.index = new_index
    
    # Handle duplicates if any
    if result_df.index.duplicated().any():
        print(f"Warning: Found duplicate standardized taxa in {method} data, aggregating by sum")
        result_df = result_df.groupby(level=0).sum()
    
    return result_df

def find_common_samples(df1, df2):
    """
    Find common samples between two DataFrames.
    
    Parameters:
    -----------
    df1, df2 : pandas.DataFrame
        DataFrames to compare
        
    Returns:
    --------
    list
        List of common sample IDs
    """
    samples1 = set(df1.columns)
    samples2 = set(df2.columns)
    common = samples1.intersection(samples2)
    return sorted(list(common))

def calculate_concordance(sylph_df, metaphlan_df, common_samples):
    """
    Calculate concordance between Sylph and MetaPhlAn results.
    
    Parameters:
    -----------
    sylph_df : pandas.DataFrame
        Abundance table from Sylph
    metaphlan_df : pandas.DataFrame
        Abundance table from MetaPhlAn
    common_samples : list
        List of common sample IDs
        
    Returns:
    --------
    dict
        Dictionary with concordance statistics
    """
    # Filter to common samples
    sylph_subset = sylph_df[common_samples]
    metaphlan_subset = metaphlan_df[common_samples]
    
    # Get present taxa (non-zero abundance) in each dataset
    sylph_present = set()
    metaphlan_present = set()
    
    for sample in common_samples:
        sylph_taxa = sylph_subset[sylph_subset[sample] > 0].index
        metaphlan_taxa = metaphlan_subset[metaphlan_subset[sample] > 0].index
        
        sylph_present.update(sylph_taxa)
        metaphlan_present.update(metaphlan_taxa)
    
    # Calculate overlap
    all_taxa = sylph_present.union(metaphlan_present)
    common_taxa = sylph_present.intersection(metaphlan_present)
    sylph_only = sylph_present - metaphlan_present
    metaphlan_only = metaphlan_present - sylph_present
    
    # Calculate Jaccard index (measure of taxonomic overlap)
    jaccard = len(common_taxa) / len(all_taxa) if all_taxa else 0
    
    results = {
        'total_taxa': len(all_taxa),
        'common_taxa': len(common_taxa),
        'sylph_only_taxa': len(sylph_only),
        'metaphlan_only_taxa': len(metaphlan_only),
        'jaccard_index': jaccard,
        'sylph_only_taxa_list': sorted(list(sylph_only)),
        'metaphlan_only_taxa_list': sorted(list(metaphlan_only)),
        'common_taxa_list': sorted(list(common_taxa))
    }
    
    return results

def compare_abundance_correlation(sylph_df, metaphlan_df, common_samples, common_taxa):
    """
    Compare abundance correlations between Sylph and MetaPhlAn for common taxa.
    
    Parameters:
    -----------
    sylph_df : pandas.DataFrame
        Abundance table from Sylph
    metaphlan_df : pandas.DataFrame
        Abundance table from MetaPhlAn
    common_samples : list
        List of common sample IDs
    common_taxa : list
        List of common taxa between methods
        
    Returns:
    --------
    dict
        Dictionary with correlation statistics
    """
    # Check if there are any common taxa
    if not common_taxa:
        print("No common taxa found between Sylph and MetaPhlAn. Skipping correlation analysis.")
        return {
            'taxa_correlations': pd.DataFrame(),
            'mean_pearson': float('nan'),
            'mean_spearman': float('nan'),
            'significant_correlations': 0
        }
    
    # Filter to common samples and taxa
    sylph_subset = sylph_df.loc[sylph_df.index.isin(common_taxa), common_samples]
    metaphlan_subset = metaphlan_df.loc[metaphlan_df.index.isin(common_taxa), common_samples]
    
    # Ensure we're comparing the same taxa in the same order
    common_taxa_present = sorted(list(set(sylph_subset.index).intersection(set(metaphlan_subset.index))))
    
    if not common_taxa_present:
        print("No common taxa found in both datasets after filtering. Skipping correlation analysis.")
        return {
            'taxa_correlations': pd.DataFrame(),
            'mean_pearson': float('nan'),
            'mean_spearman': float('nan'),
            'significant_correlations': 0
        }
    
    sylph_subset = sylph_subset.loc[common_taxa_present]
    metaphlan_subset = metaphlan_subset.loc[common_taxa_present]
    
    # Calculate correlations for each taxon across samples
    taxa_correlations = []
    for taxon in common_taxa_present:
        sylph_values = sylph_subset.loc[taxon].values
        metaphlan_values = metaphlan_subset.loc[taxon].values
        
        # Calculate Pearson and Spearman correlations
        try:
            pearson, p_pearson = pearsonr(sylph_values, metaphlan_values)
        except:
            pearson, p_pearson = np.nan, np.nan
            
        try:
            spearman, p_spearman = spearmanr(sylph_values, metaphlan_values)
        except:
            spearman, p_spearman = np.nan, np.nan
        
        taxa_correlations.append({
            'taxon': taxon,
            'pearson_r': pearson,
            'pearson_p': p_pearson,
            'spearman_r': spearman,
            'spearman_p': p_spearman,
            'sylph_mean': np.mean(sylph_values),
            'metaphlan_mean': np.mean(metaphlan_values),
            'sylph_nonzero': np.sum(sylph_values > 0),
            'metaphlan_nonzero': np.sum(metaphlan_values > 0)
        })
    
    # Create a DataFrame of correlation results
    correlation_df = pd.DataFrame(taxa_correlations)
    
    # Calculate summary statistics, handling the case of an empty DataFrame
    if correlation_df.empty:
        mean_pearson = float('nan')
        mean_spearman = float('nan')
        significant_correlations = 0
    else:
        mean_pearson = correlation_df['pearson_r'].mean() if 'pearson_r' in correlation_df else float('nan')
        mean_spearman = correlation_df['spearman_r'].mean() if 'spearman_r' in correlation_df else float('nan')
        significant_correlations = sum((correlation_df['pearson_p'] < 0.05) & (~pd.isna(correlation_df['pearson_p']))) if 'pearson_p' in correlation_df else 0
    
    return {
        'taxa_correlations': correlation_df,
        'mean_pearson': mean_pearson,
        'mean_spearman': mean_spearman,
        'significant_correlations': significant_correlations
    }


def plot_abundance_comparison(sylph_df, metaphlan_df, common_samples, common_taxa, output_dir):
    """
    Create plots comparing abundance distributions between methods.
    
    Parameters:
    -----------
    sylph_df : pandas.DataFrame
        Abundance table from Sylph
    metaphlan_df : pandas.DataFrame
        Abundance table from MetaPhlAn
    common_samples : list
        List of common sample IDs
    common_taxa : list
        List of common taxa between methods
    output_dir : str
        Directory to save plots
        
    Returns:
    --------
    None
    """
    # Create figures directory
    figures_dir = os.path.join(output_dir, 'figures')
    os.makedirs(figures_dir, exist_ok=True)
    
    # Filter to common samples and taxa
    common_taxa_list = sorted(list(set(common_taxa).intersection(set(sylph_df.index)).intersection(set(metaphlan_df.index))))
    
    if not common_taxa_list:
        print("No common taxa found for plotting")
        return
    
    # Filter to top 10 most abundant taxa in either method
    sylph_mean = sylph_df.loc[common_taxa_list, common_samples].mean(axis=1)
    metaphlan_mean = metaphlan_df.loc[common_taxa_list, common_samples].mean(axis=1)
    
    # Combine and sort by maximum mean abundance
    combined_mean = pd.DataFrame({
        'Sylph': sylph_mean,
        'MetaPhlAn': metaphlan_mean
    })
    
    combined_mean['max'] = combined_mean.max(axis=1)
    top_taxa = combined_mean.sort_values('max', ascending=False).head(10).index.tolist()
    
    # 1. Scatterplot of mean abundances
    plt.figure(figsize=(10, 8))
    plt.scatter(combined_mean['Sylph'], combined_mean['MetaPhlAn'], alpha=0.7)
    
    # Highlight top taxa
    for taxon in top_taxa:
        x = combined_mean.loc[taxon, 'Sylph']
        y = combined_mean.loc[taxon, 'MetaPhlAn']
        plt.scatter(x, y, s=100, color='red')
        plt.annotate(taxon, (x, y), textcoords="offset points", xytext=(0,10), ha='center')
    
    plt.title('Mean Abundance Comparison: Sylph vs. MetaPhlAn')
    plt.xlabel('Mean Abundance (%) - Sylph')
    plt.ylabel('Mean Abundance (%) - MetaPhlAn')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(True, alpha=0.3)
    
    # Add identity line
    min_val = min(combined_mean['Sylph'].min(), combined_mean['MetaPhlAn'].min())
    max_val = max(combined_mean['Sylph'].max(), combined_mean['MetaPhlAn'].max())
    plt.plot([min_val, max_val], [min_val, max_val], '--', color='gray')
    
    plt.savefig(os.path.join(figures_dir, 'mean_abundance_comparison.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Bar plot of top taxa
    plt.figure(figsize=(12, 8))
    
    data_to_plot = combined_mean.loc[top_taxa, ['Sylph', 'MetaPhlAn']]
    data_to_plot = data_to_plot.sort_values('Sylph', ascending=False)
    
    # Create positions for the bars
    ind = np.arange(len(top_taxa))
    width = 0.35
    
    plt.bar(ind - width/2, data_to_plot['Sylph'], width, label='Sylph')
    plt.bar(ind + width/2, data_to_plot['MetaPhlAn'], width, label='MetaPhlAn')
    
    plt.ylabel('Mean Abundance (%)')
    plt.title('Top Taxa Abundance Comparison')
    plt.xticks(ind, data_to_plot.index, rotation=45, ha='right')
    plt.legend()
    plt.tight_layout()
    
    plt.savefig(os.path.join(figures_dir, 'top_taxa_abundance_comparison.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. Correlation heatmap
    # Filter to top 20 taxa for better visualization
    corr_data = []
    
    top20_taxa = combined_mean.sort_values('max', ascending=False).head(20).index.tolist()
    
    for taxon in top20_taxa:
        for sample in common_samples:
            if taxon in sylph_df.index and taxon in metaphlan_df.index:
                corr_data.append({
                    'Taxon': taxon,
                    'Sample': sample,
                    'Sylph': sylph_df.loc[taxon, sample],
                    'MetaPhlAn': metaphlan_df.loc[taxon, sample]
                })
    
    corr_df = pd.DataFrame(corr_data)
    
    # Reshape for heatmap (pivot table)
    pivot_df = pd.pivot_table(corr_df, values=['Sylph', 'MetaPhlAn'], index='Taxon', columns='Sample')
    
    # Calculate correlation between methods for each sample
    sample_correlations = {}
    for sample in common_samples:
        sylph_vals = corr_df[corr_df['Sample'] == sample]['Sylph'].values
        metaphlan_vals = corr_df[corr_df['Sample'] == sample]['MetaPhlAn'].values
        
        try:
            corr, _ = spearmanr(sylph_vals, metaphlan_vals)
            sample_correlations[sample] = corr
        except:
            sample_correlations[sample] = np.nan
    
    # Create a heatmap of correlations
    sample_corr_df = pd.DataFrame({'Sample': list(sample_correlations.keys()),
                                   'Correlation': list(sample_correlations.values())})
    
    plt.figure(figsize=(10, 6))
    sns.barplot(x='Sample', y='Correlation', data=sample_corr_df.sort_values('Correlation', ascending=False))
    plt.title('Spearman Correlation Between Sylph and MetaPhlAn by Sample')
    plt.xticks(rotation=90)
    plt.tight_layout()
    
    plt.savefig(os.path.join(figures_dir, 'sample_correlation_comparison.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # 4. Boxplot of correlation distribution
    plt.figure(figsize=(8, 6))
    
    from scipy.stats import iqr
    # Calculate the IQR of correlations
    correlation_iqr = iqr(sample_corr_df['Correlation'].dropna())
    # Calculate the bin width using Freedman-Diaconis rule
    bin_width = 2 * correlation_iqr / (len(sample_corr_df['Correlation'].dropna()) ** (1/3))
    # Calculate the number of bins
    num_bins = int(np.ceil((sample_corr_df['Correlation'].max() - sample_corr_df['Correlation'].min()) / bin_width))
    
    plt.hist(sample_corr_df['Correlation'].dropna(), bins=max(5, min(num_bins, 20)), alpha=0.7, color='skyblue')
    plt.axvline(sample_corr_df['Correlation'].mean(), color='red', linestyle='--', 
                label=f'Mean: {sample_corr_df["Correlation"].mean():.2f}')
    
    plt.title('Distribution of Sylph-MetaPhlAn Correlations Across Samples')
    plt.xlabel('Spearman Correlation')
    plt.ylabel('Number of Samples')
    plt.legend()
    plt.grid(alpha=0.3)
    
    plt.savefig(os.path.join(figures_dir, 'correlation_distribution.png'), dpi=300, bbox_inches='tight')
    plt.close()

def load_metadata(metadata_file, sample_id_column='SampleID'):
    """
    Load and process metadata file.
    
    Parameters:
    -----------
    metadata_file : str
        Path to metadata file
    sample_id_column : str
        Column name for sample IDs
        
    Returns:
    --------
    pandas.DataFrame
        Metadata DataFrame with sample IDs as index
    """
    try:
        metadata_df = pd.read_csv(metadata_file)
        metadata_df.set_index(sample_id_column, inplace=True)
        return metadata_df
    except Exception as e:
        print(f"Error loading metadata file: {str(e)}")
        return pd.DataFrame()

def create_venn_diagram(sylph_taxa, metaphlan_taxa, output_file):
    """
    Create a Venn diagram of taxa detected by each method.
    
    Parameters:
    -----------
    sylph_taxa : set
        Set of taxa detected by Sylph
    metaphlan_taxa : set
        Set of taxa detected by MetaPhlAn
    output_file : str
        Path to save the Venn diagram
        
    Returns:
    --------
    None
    """
    try:
        from matplotlib_venn import venn2
        
        plt.figure(figsize=(8, 6))
        venn = venn2([sylph_taxa, metaphlan_taxa], ('Sylph', 'MetaPhlAn'))
        
        # Add counts to labels
        venn.get_label_by_id('10').set_text(f'Sylph only\n{len(sylph_taxa - metaphlan_taxa)}')
        venn.get_label_by_id('01').set_text(f'MetaPhlAn only\n{len(metaphlan_taxa - sylph_taxa)}')
        venn.get_label_by_id('11').set_text(f'Both\n{len(sylph_taxa & metaphlan_taxa)}')
        
        plt.title('Taxa Detected by Each Method')
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
    except ImportError:
        print("Warning: matplotlib_venn not installed. Skipping Venn diagram.")
    except Exception as e:
        print(f"Error creating Venn diagram: {str(e)}")

def main():
    """Main function to compare Sylph and MetaPhlAn results."""
    # Parse arguments
    args = parse_args()
    
    # Set up directories
    os.makedirs(args.output_dir, exist_ok=True)
    figures_dir = os.path.join(args.output_dir, 'figures')
    os.makedirs(figures_dir, exist_ok=True)
    
    # Load Sylph data
    sylph_file = os.path.join(args.sylph_dir, 'sylph_bacteria_abundance.csv')
    sylph_df = load_abundance_table(sylph_file)
    
    if sylph_df.empty:
        print("Error: Could not load Sylph data")
        return
    
    # Load MetaPhlAn data
    metaphlan_file = os.path.join(args.metaphlan_dir, 'combined_abundance.csv')
    metaphlan_df = load_abundance_table(metaphlan_file)
    
    if metaphlan_df.empty:
        print("Error: Could not load MetaPhlAn data")
        return
    
    # Standardize taxonomy names
    print("Standardizing taxonomy names for comparison...")
    sylph_std = standardize_taxonomy(sylph_df, 'sylph')
    metaphlan_std = standardize_taxonomy(metaphlan_df, 'metaphlan')
    
    # Find common samples
    common_samples = find_common_samples(sylph_std, metaphlan_std)
    print(f"Found {len(common_samples)} samples common to both datasets")
    
    if not common_samples:
        print("Error: No common samples found")
        return
    
    # Calculate concordance
    print("Calculating concordance between methods...")
    concordance = calculate_concordance(sylph_std, metaphlan_std, common_samples)
    
    print(f"Total unique taxa: {concordance['total_taxa']}")
    print(f"Common to both methods: {concordance['common_taxa']} ({concordance['common_taxa']/concordance['total_taxa']*100:.1f}%)")
    print(f"Sylph only: {concordance['sylph_only_taxa']} ({concordance['sylph_only_taxa']/concordance['total_taxa']*100:.1f}%)")
    print(f"MetaPhlAn only: {concordance['metaphlan_only_taxa']} ({concordance['metaphlan_only_taxa']/concordance['total_taxa']*100:.1f}%)")
    print(f"Jaccard similarity index: {concordance['jaccard_index']:.3f}")
    
    # Create a Venn diagram
    venn_file = os.path.join(figures_dir, 'taxa_venn_diagram.png')
    create_venn_diagram(
        set(concordance['sylph_only_taxa_list'] + concordance['common_taxa_list']),
        set(concordance['metaphlan_only_taxa_list'] + concordance['common_taxa_list']),
        venn_file
    )
    
    # Compare abundance correlations for common taxa
    print("Comparing abundance correlations for common taxa...")
    correlations = compare_abundance_correlation(
        sylph_std, 
        metaphlan_std, 
        common_samples, 
        concordance['common_taxa_list']
    )
    
    print(f"Mean Pearson correlation: {correlations['mean_pearson']:.3f}")
    print(f"Mean Spearman correlation: {correlations['mean_spearman']:.3f}")
    print(f"Number of significantly correlated taxa: {correlations['significant_correlations']}")
    
    # Save correlation results
    correlation_file = os.path.join(args.output_dir, 'taxa_correlations.csv')
    correlations['taxa_correlations'].to_csv(correlation_file)
    print(f"Correlation results saved to {correlation_file}")
    
    # Create plots
    print("Creating comparison plots...")
    plot_abundance_comparison(
        sylph_std, 
        metaphlan_std, 
        common_samples, 
        concordance['common_taxa_list'],
        args.output_dir
    )
    
    # Save concordance results
    concordance_file = os.path.join(args.output_dir, 'method_concordance.txt')
    with open(concordance_file, 'w') as f:
        f.write("COMPARISON BETWEEN SYLPH AND METAPHLAN\n")
        f.write("=====================================\n\n")
        f.write(f"Total unique taxa: {concordance['total_taxa']}\n")
        f.write(f"Common to both methods: {concordance['common_taxa']} ({concordance['common_taxa']/concordance['total_taxa']*100:.1f}%)\n")
        f.write(f"Sylph only: {concordance['sylph_only_taxa']} ({concordance['sylph_only_taxa']/concordance['total_taxa']*100:.1f}%)\n")
        f.write(f"MetaPhlAn only: {concordance['metaphlan_only_taxa']} ({concordance['metaphlan_only_taxa']/concordance['total_taxa']*100:.1f}%)\n")
        f.write(f"Jaccard similarity index: {concordance['jaccard_index']:.3f}\n\n")
        
        f.write("CORRELATION METRICS\n")
        f.write("==================\n\n")
        f.write(f"Mean Pearson correlation: {correlations['mean_pearson']:.3f}\n")
        f.write(f"Mean Spearman correlation: {correlations['mean_spearman']:.3f}\n")
        f.write(f"Number of significantly correlated taxa: {correlations['significant_correlations']}\n\n")
        
        f.write("TOP 10 CORRELATED TAXA\n")
        f.write("====================\n\n")
        top_corr = correlations['taxa_correlations'].sort_values('spearman_r', ascending=False).head(10)
        for _, row in top_corr.iterrows():
            f.write(f"{row['taxon']}: Spearman r = {row['spearman_r']:.3f}, Pearson r = {row['pearson_r']:.3f}\n")
        
        f.write("\n\nTOP 10 MOST ABUNDANT TAXA IN SYLPH\n")
        f.write("===============================\n\n")
        sylph_top = sylph_std[common_samples].mean(axis=1).nlargest(10)
        for taxon, value in sylph_top.items():
            in_metaphlan = "    (also in MetaPhlAn)" if taxon in metaphlan_std.index else ""
            f.write(f"{taxon}: {value:.3f}%{in_metaphlan}\n")
        
        f.write("\n\nTOP 10 MOST ABUNDANT TAXA IN METAPHLAN\n")
        f.write("====================================\n\n")
        metaphlan_top = metaphlan_std[common_samples].mean(axis=1).nlargest(10)
        for taxon, value in metaphlan_top.items():
            in_sylph = "    (also in Sylph)" if taxon in sylph_std.index else ""
            f.write(f"{taxon}: {value:.3f}%{in_sylph}\n")
    
    print(f"Concordance results saved to {concordance_file}")
    
    # If metadata is available, create plots by groups
    if os.path.exists(args.metadata):
        print("Loading metadata for group analysis...")
        metadata_df = load_metadata(args.metadata, args.sample_id_column)
        
        if not metadata_df.empty:
            # Check for timing variable
            if 'Timing' in metadata_df.columns:
                print("Creating analysis by timing...")
                
                # Filter to samples with timing information
                common_with_timing = [s for s in common_samples if s in metadata_df.index and not pd.isna(metadata_df.loc[s, 'Timing'])]
                
                if common_with_timing:
                    # Group samples by timing
                    for timing, group in metadata_df.loc[common_with_timing].groupby('Timing'):
                        timing_samples = group.index.tolist()
                        
                        # Calculate correlation for this timing
                        sylph_timing = sylph_std[timing_samples].mean(axis=1)
                        metaphlan_timing = metaphlan_std[timing_samples].mean(axis=1)
                        
                        # Common taxa with non-zero abundance
                        common_taxa_timing = set(sylph_timing[sylph_timing > 0].index).intersection(
                            set(metaphlan_timing[metaphlan_timing > 0].index))
                        
                        if common_taxa_timing:
                            # Filter to common taxa
                            sylph_filtered = sylph_timing[sylph_timing.index.isin(common_taxa_timing)]
                            metaphlan_filtered = metaphlan_timing[metaphlan_timing.index.isin(common_taxa_timing)]
                            
                            # Ensure same order
                            common_taxa_list = sorted(list(common_taxa_timing))
                            sylph_filtered = sylph_filtered.loc[common_taxa_list]
                            metaphlan_filtered = metaphlan_filtered.loc[common_taxa_list]
                            
                            # Calculate correlation
                            try:
                                spearman, p_spearman = spearmanr(sylph_filtered, metaphlan_filtered)
                                
                                print(f"Timing {timing}: Spearman correlation = {spearman:.3f} (p = {p_spearman:.3f})")
                                
                                # Create scatterplot
                                plt.figure(figsize=(8, 6))
                                plt.scatter(sylph_filtered, metaphlan_filtered, alpha=0.7)
                                
                                # Add taxon labels for top 5
                                combined = pd.DataFrame({
                                    'Sylph': sylph_filtered,
                                    'MetaPhlAn': metaphlan_filtered
                                })
                                combined['max'] = combined.max(axis=1)
                                top5 = combined.sort_values('max', ascending=False).head(5).index
                                
                                for taxon in top5:
                                    x = combined.loc[taxon, 'Sylph']
                                    y = combined.loc[taxon, 'MetaPhlAn']
                                    plt.annotate(taxon, (x, y), textcoords="offset points", xytext=(0,10), ha='center')
                                
                                plt.title(f'Abundance Comparison for Timing "{timing}"')
                                plt.xlabel('Mean Abundance (%) - Sylph')
                                plt.ylabel('Mean Abundance (%) - MetaPhlAn')
                                plt.text(0.05, 0.95, f'Spearman r = {spearman:.3f}\np = {p_spearman:.3f}',
                                        transform=plt.gca().transAxes, verticalalignment='top',
                                        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
                                plt.grid(alpha=0.3)
                                
                                # Add identity line
                                max_val = max(sylph_filtered.max(), metaphlan_filtered.max())
                                plt.plot([0, max_val], [0, max_val], '--', color='gray')
                                
                                plt.savefig(os.path.join(figures_dir, f'abundance_comparison_timing_{timing}.png'), 
                                           dpi=300, bbox_inches='tight')
                                plt.close()
                                
                            except Exception as e:
                                print(f"Error calculating correlation for timing {timing}: {str(e)}")
            
            # Check for symptoms variable
            if 'Symptoms' in metadata_df.columns:
                print("Creating analysis by symptoms...")
                
                # Filter to samples with symptoms information
                common_with_symptoms = [s for s in common_samples if s in metadata_df.index and not pd.isna(metadata_df.loc[s, 'Symptoms'])]
                
                if common_with_symptoms:
                    # Similar analysis as above but for symptoms...
                    for symptom, group in metadata_df.loc[common_with_symptoms].groupby('Symptoms'):
                        symptom_samples = group.index.tolist()
                        
                        # Calculate correlation for this symptom group
                        sylph_symptom = sylph_std[symptom_samples].mean(axis=1)
                        metaphlan_symptom = metaphlan_std[symptom_samples].mean(axis=1)
                        
                        # Common taxa with non-zero abundance
                        common_taxa_symptom = set(sylph_symptom[sylph_symptom > 0].index).intersection(
                            set(metaphlan_symptom[metaphlan_symptom > 0].index))
                        
                        if common_taxa_symptom:
                            # Create similar plot as for timing
                            # (code similar to the timing section)
                            # Filter to common taxa
                            sylph_filtered = sylph_symptom[sylph_symptom.index.isin(common_taxa_symptom)]
                            metaphlan_filtered = metaphlan_symptom[metaphlan_symptom.index.isin(common_taxa_symptom)]
                            
                            # Ensure same order
                            common_taxa_list = sorted(list(common_taxa_symptom))
                            sylph_filtered = sylph_filtered.loc[common_taxa_list]
                            metaphlan_filtered = metaphlan_filtered.loc[common_taxa_list]
                            
                            # Calculate correlation
                            try:
                                spearman, p_spearman = spearmanr(sylph_filtered, metaphlan_filtered)
                                
                                print(f"Symptom {symptom}: Spearman correlation = {spearman:.3f} (p = {p_spearman:.3f})")
                                
                                # Create scatterplot
                                plt.figure(figsize=(8, 6))
                                plt.scatter(sylph_filtered, metaphlan_filtered, alpha=0.7)
                                
                                # Add taxon labels for top 5
                                combined = pd.DataFrame({
                                    'Sylph': sylph_filtered,
                                    'MetaPhlAn': metaphlan_filtered
                                })
                                combined['max'] = combined.max(axis=1)
                                top5 = combined.sort_values('max', ascending=False).head(5).index
                                
                                for taxon in top5:
                                    x = combined.loc[taxon, 'Sylph']
                                    y = combined.loc[taxon, 'MetaPhlAn']
                                    plt.annotate(taxon, (x, y), textcoords="offset points", xytext=(0,10), ha='center')
                                
                                plt.title(f'Abundance Comparison for Symptom "{symptom}"')
                                plt.xlabel('Mean Abundance (%) - Sylph')
                                plt.ylabel('Mean Abundance (%) - MetaPhlAn')
                                plt.text(0.05, 0.95, f'Spearman r = {spearman:.3f}\np = {p_spearman:.3f}',
                                        transform=plt.gca().transAxes, verticalalignment='top',
                                        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
                                plt.grid(alpha=0.3)
                                
                                # Add identity line
                                max_val = max(sylph_filtered.max(), metaphlan_filtered.max())
                                plt.plot([0, max_val], [0, max_val], '--', color='gray')
                                
                                plt.savefig(os.path.join(figures_dir, f'abundance_comparison_symptom_{symptom}.png'), 
                                           dpi=300, bbox_inches='tight')
                                plt.close()
                                
                            except Exception as e:
                                print(f"Error calculating correlation for symptom {symptom}: {str(e)}")
    
    print("\nComparison complete!")

if __name__ == "__main__":
    main()
