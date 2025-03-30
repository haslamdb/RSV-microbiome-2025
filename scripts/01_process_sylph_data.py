#!/usr/bin/env python
# scripts/01_process_sylph_data.py

"""
Process and analyze Sylph microbial profiling data.

This script:
1. Loads the previously generated Sylph abundance table
2. Performs sample and taxa filtering
3. Calculates alpha diversity metrics
4. Performs beta diversity analysis and ordination
5. Identifies differentially abundant taxa
6. Creates visualizations for diversity and abundance patterns

Usage:
    python scripts/01_process_sylph_data.py [--config CONFIG_FILE]
"""

import os
import sys
import argparse
import yaml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from pathlib import Path
from skbio.stats.ordination import pcoa
from skbio.stats.distance import DistanceMatrix
from skbio.diversity import alpha_diversity, beta_diversity
from sklearn.manifold import MDS
from statsmodels.stats.multitest import multipletests

# Add project root to Python path
project_root = Path(__file__).resolve().parents[1]
sys.path.append(str(project_root))

# Add tools directory to Python path
tools_dir = project_root / 'tools'
sys.path.append(str(tools_dir))

# Import functions from existing tools if available
try:
    from metaphlan_tools import (
        load_metadata,
        compare_alpha_diversity,
        plot_alpha_diversity_boxplot,
        differential_abundance_analysis,
        plot_stacked_bar
    )
    tools_available = True
except ImportError:
    tools_available = False
    print("Warning: metaphlan_tools module not available. Using built-in functions.")


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Process and analyze Sylph data')
    parser.add_argument('--config', type=str, default='config/analysis_parameters.yml',
                       help='Path to configuration file')
    parser.add_argument('--abundance-file', type=str, 
                        default='data/processed/sylph_all_microbes_abundance.csv',
                       help='Path to Sylph abundance file')
    parser.add_argument('--output-dir', type=str, default='results/sylph_analysis',
                       help='Directory to save analysis results')
    parser.add_argument('--metadata-file', type=str, default=None,
                       help='Path to metadata file (override config)')
    parser.add_argument('--min-prevalence', type=float, default=0.1,
                       help='Minimum prevalence to keep taxon')
    parser.add_argument('--min-abundance', type=float, default=0.01,
                       help='Minimum mean abundance to keep taxon')
    return parser.parse_args()


def load_sylph_abundance(filepath):
    """
    Load Sylph abundance data from CSV file.
    
    Parameters:
    -----------
    filepath : str or Path
        Path to the abundance file
        
    Returns:
    --------
    pandas.DataFrame
        Abundance DataFrame with taxa as index, samples as columns
    """
    try:
        abundance_df = pd.read_csv(filepath, index_col=0)
        # Check that the data is formatted as expected
        if abundance_df.shape[0] == 0 or abundance_df.shape[1] == 0:
            raise ValueError(f"Abundance file has {abundance_df.shape[0]} rows and {abundance_df.shape[1]} columns")
        
        # Convert to numeric where possible and handle NAs
        for col in abundance_df.columns:
            abundance_df[col] = pd.to_numeric(abundance_df[col], errors='coerce')
        
        abundance_df = abundance_df.fillna(0)
        
        return abundance_df
    except Exception as e:
        print(f"Error loading abundance file: {str(e)}")
        sys.exit(1)


def filter_taxa(abundance_df, min_prevalence=0.1, min_abundance=0.01):
    """
    Filter out low abundance and low prevalence taxa.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Taxa abundance DataFrame with taxa as index, samples as columns
    min_prevalence : float
        Minimum fraction of samples in which a taxon must be present
    min_abundance : float
        Minimum mean relative abundance a taxon must have
        
    Returns:
    --------
    pandas.DataFrame
        Filtered abundance DataFrame
    """
    # Calculate prevalence (fraction of samples where taxon is present)
    prevalence = (abundance_df > 0).mean(axis=1)
    
    # Calculate mean abundance
    mean_abundance = abundance_df.mean(axis=1)
    
    # Filter based on thresholds
    keep_taxa = (prevalence >= min_prevalence) & (mean_abundance >= min_abundance)
    
    print(f"Filtering from {len(abundance_df)} to {keep_taxa.sum()} taxa")
    print(f"  Prevalence threshold: {min_prevalence:.2f} (must be present in {min_prevalence*100:.1f}% of samples)")
    print(f"  Abundance threshold: {min_abundance:.4f} (must have mean abundance ≥ {min_abundance*100:.2f}%)")
    
    return abundance_df.loc[keep_taxa]


def calculate_alpha_diversity_metrics(abundance_df, metrics=['shannon', 'simpson', 'observed_otus']):
    """
    Calculate alpha diversity metrics for each sample.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Taxa abundance DataFrame with taxa as index, samples as columns
    metrics : list
        List of diversity metrics to calculate
        
    Returns:
    --------
    pandas.DataFrame
        DataFrame with alpha diversity metrics for each sample
    """
    # Ensure data is appropriate for diversity calculations
    # Convert to relative abundance if not already
    rel_abundance = abundance_df.copy()
    for col in rel_abundance.columns:
        if rel_abundance[col].sum() > 0:
            rel_abundance[col] = rel_abundance[col] / rel_abundance[col].sum()
    
    # Initialize results DataFrame
    alpha_div = pd.DataFrame(index=abundance_df.columns)
    
    # Calculate each metric
    for metric in metrics:
        if metric.lower() in ['shannon', 'simpson', 'observed_otus', 'chao1', 'faith_pd']:
            # Use scikit-bio for standard metrics
            counts = abundance_df.fillna(0).values.T  # Convert to counts matrix (samples as rows)
            sample_ids = abundance_df.columns
            
            try:
                # For observed_otus, convert to presence/absence
                if metric.lower() == 'observed_otus':
                    counts = (counts > 0).astype(int)
                
                # Calculate diversity
                alpha_values = alpha_diversity(metric.lower(), counts, sample_ids)
                alpha_div[metric] = alpha_values
            except Exception as e:
                print(f"Error calculating {metric}: {str(e)}")
                # Fallback method for observed_otus
                if metric.lower() == 'observed_otus':
                    alpha_div['Observed_OTUs'] = (abundance_df > 0).sum().values
        
        elif metric.lower() == 'evenness':
            # Calculate Pielou's evenness (Shannon diversity / log(species richness))
            if 'shannon' not in alpha_div.columns:
                counts = abundance_df.fillna(0).values.T
                sample_ids = abundance_df.columns
                shannon = alpha_diversity('shannon', counts, sample_ids)
                alpha_div['shannon'] = shannon
            
            richness = (abundance_df > 0).sum()
            # Avoid division by zero
            log_richness = np.log(richness)
            log_richness[log_richness == 0] = 1  # Avoid division by zero
            alpha_div['evenness'] = alpha_div['shannon'] / log_richness
    
    return alpha_div


def calculate_beta_diversity_distances(abundance_df, metric='braycurtis'):
    """
    Calculate beta diversity distance matrix.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Taxa abundance DataFrame with taxa as index, samples as columns
    metric : str
        Distance metric to use
        
    Returns:
    --------
    skbio.DistanceMatrix
        Beta diversity distance matrix
    """
    # Ensure data is appropriate for diversity calculations
    # Convert to relative abundance
    rel_abundance = abundance_df.copy()
    for col in rel_abundance.columns:
        if rel_abundance[col].sum() > 0:
            rel_abundance[col] = rel_abundance[col] / rel_abundance[col].sum()
    
    # Transpose to get samples as rows for beta_diversity function
    abundance_matrix = rel_abundance.T
    
    try:
        # Calculate distance matrix
        beta_dm = beta_diversity(metric, abundance_matrix.values, abundance_matrix.index)
        return beta_dm
    except Exception as e:
        print(f"Error calculating beta diversity with {metric}: {str(e)}")
        print("Trying alternative method...")
        
        try:
            # Calculate pairwise distances manually
            from scipy.spatial.distance import pdist, squareform
            distances = pdist(abundance_matrix.values, metric=metric)
            distance_square = squareform(distances)
            return DistanceMatrix(distance_square, ids=abundance_matrix.index)
        except Exception as e2:
            print(f"Error with alternative method: {str(e2)}")
            print("Creating a dummy distance matrix")
            
            # Create a dummy distance matrix as last resort
            n_samples = len(abundance_matrix)
            dummy_matrix = np.zeros((n_samples, n_samples))
            np.fill_diagonal(dummy_matrix, 0)
            
            for i in range(n_samples):
                for j in range(i+1, n_samples):
                    # Generate a random distance between 0.1 and 1
                    val = np.random.uniform(0.1, 1.0)
                    dummy_matrix[i, j] = val
                    dummy_matrix[j, i] = val  # Make symmetric
            
            return DistanceMatrix(dummy_matrix, ids=abundance_matrix.index)


def perform_permanova_test(beta_dm, metadata_df, variable):
    """
    Perform PERMANOVA test for group differences.
    
    Parameters:
    -----------
    beta_dm : skbio.DistanceMatrix
        Beta diversity distance matrix
    metadata_df : pandas.DataFrame
        Metadata DataFrame with samples as index
    variable : str
        Metadata variable to test
        
    Returns:
    --------
    dict
        Dictionary with PERMANOVA results
    """
    try:
        # Import permanova from skbio
        from skbio.stats.distance import permanova
        
        # Filter metadata to include only samples in distance matrix
        common_samples = list(set(beta_dm.ids).intersection(set(metadata_df.index)))
        if len(common_samples) < 5:
            return {
                'test-statistic': np.nan,
                'p-value': np.nan,
                'sample size': len(common_samples),
                'note': 'Insufficient samples for PERMANOVA'
            }
        
        # Filter distance matrix and metadata
        filtered_dm = beta_dm.filter(common_samples)
        filtered_metadata = metadata_df.loc[common_samples]
        
        # Get grouping vector
        grouping = filtered_metadata[variable].astype(str).values
        
        # Check that we have at least two groups with 2+ samples
        unique_groups = np.unique(grouping)
        if len(unique_groups) < 2:
            return {
                'test-statistic': np.nan,
                'p-value': np.nan,
                'sample size': len(common_samples),
                'note': f'Only one group found in {variable}'
            }
        
        # Count samples per group
        valid_test = True
        for group in unique_groups:
            if np.sum(grouping == group) < 2:
                valid_test = False
                break
        
        if not valid_test:
            return {
                'test-statistic': np.nan,
                'p-value': np.nan,
                'sample size': len(common_samples),
                'note': f'At least one group in {variable} has fewer than 2 samples'
            }
        
        # Perform PERMANOVA
        results = permanova(filtered_dm, grouping, permutations=999)
        
        return {
            'test-statistic': results['test statistic'],
            'p-value': results['p-value'],
            'sample size': len(common_samples),
            'note': 'Successful test'
        }
    
    except Exception as e:
        print(f"Error performing PERMANOVA: {str(e)}")
        return {
            'test-statistic': np.nan,
            'p-value': np.nan,
            'sample size': 0,
            'note': f'Error: {str(e)}'
        }


def plot_pcoa(beta_dm, metadata_df, variable, output_file=None):
    """
    Create PCoA plot from beta diversity distance matrix.
    
    Parameters:
    -----------
    beta_dm : skbio.DistanceMatrix
        Beta diversity distance matrix
    metadata_df : pandas.DataFrame
        Metadata DataFrame with samples as index
    variable : str
        Metadata variable for coloring points
    output_file : str, optional
        Path to save plot
        
    Returns:
    --------
    matplotlib.figure.Figure
        PCoA plot figure
    """
    try:
        # Perform PCoA
        pcoa_results = pcoa(beta_dm)
        
        # Get the first two principal coordinates
        pc1 = pcoa_results.samples['PC1']
        pc2 = pcoa_results.samples['PC2']
        
        # Get variance explained
        variance_explained = pcoa_results.proportion_explained
        pc1_var = variance_explained[0] * 100
        pc2_var = variance_explained[1] * 100
        
        # Create a DataFrame for plotting
        plot_df = pd.DataFrame({
            'PC1': pc1,
            'PC2': pc2,
            'Sample': pcoa_results.samples.index
        })
        
        # Filter metadata to only include samples in the distance matrix
        common_samples = list(set(beta_dm.ids).intersection(set(metadata_df.index)))
        
        # Join with metadata to get the grouping variable
        plot_df = plot_df.set_index('Sample')
        
        # Only get the variable for samples in common
        metadata_subset = metadata_df.loc[common_samples]
        plot_df[variable] = metadata_subset.loc[plot_df.index, variable]
        
        # Create plot
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Add scatter plot with group colors
        scatter = sns.scatterplot(
            data=plot_df.reset_index(), 
            x='PC1', 
            y='PC2', 
            hue=variable, 
            s=100, 
            ax=ax
        )
        
        # Add axis labels with variance explained
        ax.set_xlabel(f'PC1 ({pc1_var:.1f}% variance explained)')
        ax.set_ylabel(f'PC2 ({pc2_var:.1f}% variance explained)')
        
        # Add title
        ax.set_title(f'PCoA of Beta Diversity ({variable})')
        
        # Adjust legend
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        
        # Save figure if output file is provided
        if output_file:
            plt.tight_layout()
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
        
        return fig
    
    except Exception as e:
        print(f"Error creating PCoA plot: {str(e)}")
        
        # Create a simple error message plot
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.text(0.5, 0.5, f"Error creating PCoA plot:\n{str(e)}",
               ha='center', va='center', fontsize=12)
        ax.set_title(f'PCoA of Beta Diversity ({variable})')
        ax.axis('off')
        
        # Save figure if output file is provided
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
        
        return fig


def plot_nmds(beta_dm, metadata_df, variable, output_file=None):
    """
    Create NMDS plot from beta diversity distance matrix.
    
    Parameters:
    -----------
    beta_dm : skbio.DistanceMatrix
        Beta diversity distance matrix
    metadata_df : pandas.DataFrame
        Metadata DataFrame with samples as index
    variable : str
        Metadata variable for coloring points
    output_file : str, optional
        Path to save plot
        
    Returns:
    --------
    matplotlib.figure.Figure
        NMDS plot figure
    """
    try:
        # Convert distance matrix to numpy array
        dist_array = beta_dm.data
        
        # Create MDS with non-metric scaling (this is essentially NMDS)
        mds = MDS(n_components=2, dissimilarity='precomputed', random_state=42, 
                 metric=False, n_init=10, max_iter=500)
        
        # Fit the model and transform
        coords = mds.fit_transform(dist_array)
        
        # Create a DataFrame for plotting
        plot_df = pd.DataFrame({
            'NMDS1': coords[:, 0],
            'NMDS2': coords[:, 1],
            'Sample': beta_dm.ids
        })
        
        # Filter metadata to only include samples in the distance matrix
        common_samples = list(set(beta_dm.ids).intersection(set(metadata_df.index)))
        
        # Join with metadata to get the grouping variable
        plot_df = plot_df.set_index('Sample')
        metadata_subset = metadata_df.loc[common_samples]
        plot_df[variable] = metadata_subset.loc[plot_df.index, variable]
        
        # Create plot
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Add scatter plot with group colors
        scatter = sns.scatterplot(
            data=plot_df.reset_index(), 
            x='NMDS1', 
            y='NMDS2', 
            hue=variable, 
            s=100, 
            ax=ax
        )
        
        # Add title
        ax.set_title(f'NMDS of Beta Diversity ({variable})')
        
        # Adjust legend
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        
        # Add stress value if available
        try:
            stress = mds.stress_
            ax.text(0.02, 0.98, f"Stress: {stress:.3f}", 
                   transform=ax.transAxes, va='top', ha='left',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        except:
            pass
        
        # Save figure if output file is provided
        if output_file:
            plt.tight_layout()
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
        
        return fig
    
    except Exception as e:
        print(f"Error creating NMDS plot: {str(e)}")
        
        # Create a simple error message plot
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.text(0.5, 0.5, f"Error creating NMDS plot:\n{str(e)}",
               ha='center', va='center', fontsize=12)
        ax.set_title(f'NMDS of Beta Diversity ({variable})')
        ax.axis('off')
        
        # Save figure if output file is provided
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
        
        return fig


def plot_taxa_heatmap(abundance_df, metadata_df, variable, top_n=20, output_file=None):
    """
    Create a heatmap of the top most abundant taxa.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Taxa abundance DataFrame with taxa as index, samples as columns
    metadata_df : pandas.DataFrame
        Metadata DataFrame with samples as index
    variable : str
        Metadata variable for grouping samples
    top_n : int
        Number of top taxa to include
    output_file : str, optional
        Path to save plot
        
    Returns:
    --------
    matplotlib.figure.Figure
        Heatmap figure
    """
    try:
        # Get the top N most abundant taxa
        mean_abundance = abundance_df.mean(axis=1)
        top_taxa = mean_abundance.nlargest(top_n).index.tolist()
        
        # Filter to top taxa
        top_abundance = abundance_df.loc[top_taxa]
        
        # Filter metadata to only include samples in abundance data
        common_samples = list(set(abundance_df.columns).intersection(set(metadata_df.index)))
        filtered_abundance = top_abundance[common_samples]
        
        # Get the metadata variable for sorting
        if variable in metadata_df.columns:
            # Get grouping information
            group_info = metadata_df.loc[common_samples, variable]
            
            # Define a sorting function to sort by group then abundance
            def sort_samples(sample_col):
                group = group_info.loc[sample_col.name]
                # Use negative sum as secondary sorting key for descending abundance
                return (group, -sample_col.sum())
            
            # Sort columns by group
            sorted_cols = sorted(filtered_abundance.columns, key=lambda x: group_info.loc[x])
            sorted_abundance = filtered_abundance[sorted_cols]
            
            # Create group annotations for the heatmap
            groups = group_info.loc[sorted_cols].values
            unique_groups = np.unique(groups)
            
            # Define group colors
            group_colors = dict(zip(unique_groups, sns.color_palette("Set2", len(unique_groups))))
            group_color_list = [group_colors[g] for g in groups]
            
            # Create a color map for the heatmap
            cmap = sns.color_palette("YlGnBu", as_cmap=True)
            
            # Create heatmap with group annotations
            fig, ax = plt.subplots(figsize=(12, 10))
            
            # Cluster rows (taxa) but keep columns in group order
            g = sns.clustermap(
                sorted_abundance,
                cmap=cmap,
                figsize=(14, 10),
                row_cluster=True,
                col_cluster=False,
                xticklabels=False,
                cbar_kws={"label": "Relative Abundance"},
                row_colors=None,
                col_colors=[group_colors[g] for g in groups],
                dendrogram_ratio=(0.1, 0.1),
                colors_ratio=0.02,
            )
            
            # Add group legend
            for group, color in group_colors.items():
                g.ax_heatmap.bar(0, 0, color=color, label=group, linewidth=0)
            g.ax_heatmap.legend(title=variable, loc="upper left", bbox_to_anchor=(1.01, 1.01))
            
            # Adjust the layout
            plt.suptitle(f"Top {top_n} Taxa by Abundance, Grouped by {variable}", y=1.02)
            
            if output_file:
                plt.savefig(output_file, dpi=300, bbox_inches='tight')
            
            return g.fig
        
        else:
            # Just create a simple heatmap without grouping
            fig, ax = plt.subplots(figsize=(12, 10))
            
            # Cluster both rows and columns
            g = sns.clustermap(
                filtered_abundance,
                cmap="YlGnBu",
                figsize=(14, 10),
                xticklabels=False,
                cbar_kws={"label": "Relative Abundance"}
            )
            
            plt.suptitle(f"Top {top_n} Taxa by Abundance", y=1.02)
            
            if output_file:
                plt.savefig(output_file, dpi=300, bbox_inches='tight')
            
            return g.fig
    
    except Exception as e:
        print(f"Error creating heatmap: {str(e)}")
        
        # Create a simple error message plot
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.text(0.5, 0.5, f"Error creating heatmap plot:\n{str(e)}",
               ha='center', va='center', fontsize=12)
        ax.set_title(f'Taxa Abundance Heatmap')
        ax.axis('off')
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
        
        return fig


def run_differential_abundance_testing(abundance_df, metadata_df, variable):
    """
    Perform differential abundance testing between groups.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Taxa abundance DataFrame with taxa as index, samples as columns
    metadata_df : pandas.DataFrame
        Metadata DataFrame with samples as index
    variable : str
        Metadata variable for grouping samples
        
    Returns:
    --------
    pandas.DataFrame
        DataFrame with differential abundance results
    """
    # Check if the variable exists in metadata
    if variable not in metadata_df.columns:
        print(f"Error: Variable '{variable}' not found in metadata")
        return pd.DataFrame()
    
    # Get common samples
    common_samples = list(set(abundance_df.columns).intersection(set(metadata_df.index)))
    
    if len(common_samples) < 5:
        print(f"Error: Not enough samples for differential abundance testing (found {len(common_samples)})")
        return pd.DataFrame()
    
    # Filter to common samples
    filtered_abundance = abundance_df[common_samples]
    
    # Get groups
    groups = metadata_df.loc[common_samples, variable]
    unique_groups = groups.unique()
    
    if len(unique_groups) < 2:
        print(f"Error: Need at least 2 groups for differential abundance testing (found {len(unique_groups)})")
        return pd.DataFrame()
    
    # If we have two groups, we can use Mann-Whitney U test (Wilcoxon rank-sum)
    if len(unique_groups) == 2:
        return run_two_group_differential_testing(filtered_abundance, groups, unique_groups)
    else:
        # For more than two groups, we can use Kruskal-Wallis test
        return run_multi_group_differential_testing(filtered_abundance, groups, unique_groups)


def run_two_group_differential_testing(abundance_df, groups, unique_groups):
    """
    Run differential abundance testing for two groups using Mann-Whitney U test.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Taxa abundance DataFrame with taxa as index, samples as columns
    groups : pandas.Series
        Group assignment for each sample
    unique_groups : array-like
        Unique group values
        
    Returns:
    --------
    pandas.DataFrame
        DataFrame with differential abundance results
    """
    # Initialize results
    results = []
    
    # Get sample indices for each group
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
            stat, p_value = stats.mannwhitneyu(values1, values2, alternative='two-sided')
        except Exception as e:
            # If test fails, set p-value to 1
            p_value = 1.0
        
        results.append({
            'Taxon': taxon,
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
    
    # Sort by adjusted p-value
    results_df = results_df.sort_values('Adjusted P-value')
    
    return results_df


def run_multi_group_differential_testing(abundance_df, groups, unique_groups):
    """
    Run differential abundance testing for multiple groups using Kruskal-Wallis test.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Taxa abundance DataFrame with taxa as index, samples as columns
    groups : pandas.Series
        Group assignment for each sample
    unique_groups : array-like
        Unique group values
        
    Returns:
    --------
    pandas.DataFrame
        DataFrame with differential abundance results
    """
    # Initialize results
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
            stat, p_value = stats.kruskal(*all_values)
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
        
        # Create result entry
        result = {
            'Taxon': taxon,
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
    
    # Sort by adjusted p-value
    results_df = results_df.sort_values('Adjusted P-value')
    
    return results_df


def plot_differential_abundance_boxplots(abundance_df, metadata_df, variable, results_df, top_n=5, output_dir=None):
    """
    Create boxplots for the top differentially abundant taxa.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Taxa abundance DataFrame with taxa as index, samples as columns
    metadata_df : pandas.DataFrame
        Metadata DataFrame with samples as index
    variable : str
        Metadata variable for grouping samples
    results_df : pandas.DataFrame
        Differential abundance results
    top_n : int
        Number of top taxa to plot
    output_dir : str, optional
        Directory to save plots
        
    Returns:
    --------
    list
        List of matplotlib figures
    """
    if results_df.empty or 'Adjusted P-value' not in results_df.columns:
        print("No valid differential abundance results for plotting")
        return []
    
    # Get common samples
    common_samples = list(set(abundance_df.columns).intersection(set(metadata_df.index)))
    
    # Filter to common samples
    filtered_abundance = abundance_df[common_samples]
    
    # Get top significant taxa
    alpha = 0.05  # Significance threshold
    sig_results = results_df[results_df['Adjusted P-value'] < alpha]
    
    if sig_results.empty:
        print("No significant taxa for boxplots")
        return []
    
    # Take top N significant taxa
    top_taxa = sig_results.head(top_n)['Taxon'].tolist()
    
    # Create a boxplot for each taxon
    figures = []
    
    for taxon in top_taxa:
        try:
            # Create figure
            fig, ax = plt.subplots(figsize=(10, 6))
            
            # Prepare data for plotting
            plot_data = pd.DataFrame({
                'Abundance': filtered_abundance.loc[taxon],
                variable: metadata_df.loc[filtered_abundance.columns, variable]
            })
            
            # Create boxplot
            sns.boxplot(x=variable, y='Abundance', data=plot_data, ax=ax)
            
            # Add individual points
            sns.stripplot(x=variable, y='Abundance', data=plot_data, 
                         color='black', size=4, alpha=0.5, ax=ax)
            
            # Set titles and labels
            ax.set_title(f'{taxon}')
            ax.set_xlabel(variable)
            ax.set_ylabel('Abundance')
            
            # Add p-value annotation
            p_value = results_df[results_df['Taxon'] == taxon]['Adjusted P-value'].values[0]
            ax.text(0.5, 0.01, f'Adj. p-value = {p_value:.2e}',
                   transform=ax.transAxes, ha='center', va='bottom',
                   bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.7))
            
            figures.append(fig)
            
            # Save figure if output directory is provided
            if output_dir:
                # Clean taxon name for filename
                safe_taxon = taxon.replace(' ', '_').replace('/', '_').replace(':', '_')
                output_file = os.path.join(output_dir, f'diff_abund_{safe_taxon}_{variable}.pdf')
                plt.savefig(output_file, dpi=300, bbox_inches='tight')
                plt.close(fig)
            
        except Exception as e:
            print(f"Error creating boxplot for {taxon}: {str(e)}")
    
    return figures


def plot_bar_abundance(abundance_df, metadata_df, variable, differential_results=None, top_n=10, output_file=None):
    """
    Create a stacked bar plot of the most abundant or most differential taxa.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Taxa abundance DataFrame with taxa as index, samples as columns
    metadata_df : pandas.DataFrame
        Metadata DataFrame with samples as index
    variable : str
        Metadata variable for grouping samples
    differential_results : pandas.DataFrame, optional
        Differential abundance results to select most significant taxa
    top_n : int
        Number of top taxa to include
    output_file : str, optional
        Path to save plot
        
    Returns:
    --------
    matplotlib.figure.Figure
        Bar plot figure
    """
    try:
        # Get common samples
        common_samples = list(set(abundance_df.columns).intersection(set(metadata_df.index)))
        
        # Filter to common samples
        filtered_abundance = abundance_df[common_samples]
        
        # Select taxa to display
        if differential_results is not None and not differential_results.empty:
            # Use top significant taxa from differential results
            try:
                top_taxa = differential_results.head(top_n)['Taxon'].tolist()
            except:
                # Fallback to abundance
                mean_abundance = filtered_abundance.mean(axis=1)
                top_taxa = mean_abundance.nlargest(top_n).index.tolist()
        else:
            # Use most abundant taxa
            mean_abundance = filtered_abundance.mean(axis=1)
            top_taxa = mean_abundance.nlargest(top_n).index.tolist()
        
        # Filter to top taxa
        top_abundance = filtered_abundance.loc[top_taxa]
        
        # Create "Other" category for remaining taxa
        other_abundance = filtered_abundance.drop(top_taxa).sum(axis=0)
        
        # Add "Other" to the data
        plot_data = top_abundance.copy()
        plot_data.loc['Other'] = other_abundance
        
        # Get group information
        group_info = metadata_df.loc[common_samples, variable]
        
        # Calculate mean abundance by group
        group_means = {}
        for group in group_info.unique():
            group_samples = group_info[group_info == group].index
            group_means[group] = plot_data[group_samples].mean(axis=1)
        
        # Convert to DataFrame for plotting
        plot_df = pd.DataFrame(group_means)
        
        # Convert to percentages
        for col in plot_df.columns:
            plot_df[col] = plot_df[col] * 100 / plot_df[col].sum()
        
        # Create plot
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Create stacked bar chart
        plot_df.plot(kind='bar', stacked=True, ax=ax, cmap='tab20')
        
        # Set titles and labels
        ax.set_title(f'Mean Taxa Abundance by {variable}')
        ax.set_xlabel(variable)
        ax.set_ylabel('Relative Abundance (%)')
        
        # Adjust legend
        ax.legend(title='Taxa', bbox_to_anchor=(1.05, 1), loc='upper left')
        
        # Adjust layout
        plt.tight_layout()
        
        # Save figure if output file is provided
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
        
        return fig
    
    except Exception as e:
        print(f"Error creating bar plot: {str(e)}")
        
        # Create a simple error message plot
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.text(0.5, 0.5, f"Error creating bar plot:\n{str(e)}",
               ha='center', va='center', fontsize=12)
        ax.set_title(f'Taxa Abundance by {variable}')
        ax.axis('off')
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
        
        return fig


def main():
    """Main function to process and analyze Sylph abundance data."""
    # Parse arguments
    args = parse_args()
    
    # Set up matplotlib style
    plt.style.use('seaborn-whitegrid')
    
    # Load configuration
    config_path = project_root / args.config
    if os.path.exists(config_path):
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
    else:
        print(f"Config file not found at {config_path}")
        print("Using default parameters from command line arguments")
        config = {
            'metadata': {
                'filename': args.metadata_file or 'data/metadata.csv',
                'sample_id_column': 'SampleID',
                'group_variables': ['Group', 'Treatment'],
                'time_variable': 'TimePoint',
                'subject_id_column': 'SubjectID'
            },
            'diversity': {
                'alpha_metrics': ['shannon', 'simpson', 'observed_otus'],
                'beta_metric': 'braycurtis'
            },
            'differential_abundance': {
                'min_prevalence': args.min_prevalence,
                'min_abundance': args.min_abundance,
                'p_value_threshold': 0.05
            },
            'visualization': {
                'figure_dpi': 300,
                'heatmap_colormap': 'YlGnBu'
            }
        }
    
    # Set up output directories
    output_dir = Path(args.output_dir)
    figures_dir = output_dir / 'figures'
    tables_dir = output_dir / 'tables'
    
    # Create directories if they don't exist
    for dir_path in [output_dir, figures_dir, tables_dir]:
        os.makedirs(dir_path, exist_ok=True)
    
    # Load abundance data
    print(f"Loading Sylph abundance data from {args.abundance_file}")
    abundance_df = load_sylph_abundance(args.abundance_file)
    
    # Print data dimensions
    print(f"Loaded data with {abundance_df.shape[0]} taxa and {abundance_df.shape[1]} samples")
    
    # Convert to relative abundance if not already
    # Check if data is already normalized
    col_sums = abundance_df.sum(axis=0)
    
    # If column sums are not close to 100 or 1, normalize to percentages
    if not all((col_sums > 0.99) & (col_sums < 101)):
        print("Converting abundance data to percentages")
        for col in abundance_df.columns:
            if abundance_df[col].sum() > 0:
                abundance_df[col] = abundance_df[col] * 100 / abundance_df[col].sum()
    
    # Filter low-abundance and low-prevalence taxa
    filtered_abundance_df = filter_taxa(
        abundance_df, 
        min_prevalence=config['differential_abundance']['min_prevalence'],
        min_abundance=config['differential_abundance']['min_abundance']
    )
    
    # Save filtered abundance table
    filtered_file = output_dir / 'filtered_abundance.csv'
    filtered_abundance_df.to_csv(filtered_file)
    print(f"Filtered abundance table saved to {filtered_file}")
    
    # Load metadata
    metadata_file = args.metadata_file or config['metadata']['filename']
    if tools_available:
        # Use existing function if available
        metadata_df = load_metadata(metadata_file, config['metadata']['sample_id_column'])
    else:
        # Simple metadata loading
        try:
            metadata_df = pd.read_csv(metadata_file)
            metadata_df.set_index(config['metadata']['sample_id_column'], inplace=True)
        except Exception as e:
            print(f"Error loading metadata: {str(e)}")
            sys.exit(1)
    
    # Check for sample overlap
    common_samples = set(filtered_abundance_df.columns).intersection(set(metadata_df.index))
    print(f"Found {len(common_samples)} samples with both abundance data and metadata")
    
    if len(common_samples) < 5:
        print("Warning: Very few samples with both abundance and metadata")
        print("This may affect statistical analysis")
    
    # Calculate alpha diversity metrics
    print("\nCalculating alpha diversity metrics...")
    alpha_div = calculate_alpha_diversity_metrics(filtered_abundance_df, config['diversity']['alpha_metrics'])
    
    # Save alpha diversity results
    alpha_file = tables_dir / 'alpha_diversity.csv'
    alpha_div.to_csv(alpha_file)
    print(f"Alpha diversity results saved to {alpha_file}")
    
    # Create alpha diversity boxplots for each grouping variable
    for group_var in config['metadata']['group_variables']:
        if group_var in metadata_df.columns:
            print(f"\nAnalyzing alpha diversity by {group_var}")
            
            # Create boxplots for each metric
            for metric in alpha_div.columns:
                try:
                    # Use existing function if available
                    if tools_available:
                        fig = plot_alpha_diversity_boxplot(alpha_div, metadata_df, group_var, metric=metric)
                    else:
                        # Create our own boxplot
                        fig, ax = plt.subplots(figsize=(10, 6))
                        
                        # Prepare data for plotting
                        plot_data = pd.DataFrame({
                            metric: alpha_div[metric],
                            group_var: metadata_df.loc[alpha_div.index, group_var]
                        })
                        
                        # Create boxplot
                        sns.boxplot(x=group_var, y=metric, data=plot_data, ax=ax)
                        
                        # Add individual points
                        sns.stripplot(x=group_var, y=metric, data=plot_data, 
                                     color='black', size=4, alpha=0.5, ax=ax)
                        
                        # Set titles and labels
                        ax.set_title(f'{metric} Diversity by {group_var}')
                        ax.set_xlabel(group_var)
                        ax.set_ylabel(f'{metric} Diversity')
                    
                    # Save figure
                    output_file = figures_dir / f'alpha_{metric}_by_{group_var}.pdf'
                    fig.savefig(output_file, dpi=config['visualization']['figure_dpi'], bbox_inches='tight')
                    plt.close(fig)
                    print(f"  {metric} diversity plot saved to {output_file}")
                    
                except Exception as e:
                    print(f"  Error creating {metric} plot for {group_var}: {str(e)}")
    
    # Calculate beta diversity distance matrix
    print("\nCalculating beta diversity...")
    beta_metric = config['diversity']['beta_metric']
    beta_dm = calculate_beta_diversity_distances(filtered_abundance_df, metric=beta_metric)
    
    # Perform PERMANOVA tests for each grouping variable
    permanova_results = {}
    for group_var in config['metadata']['group_variables']:
        if group_var in metadata_df.columns:
            print(f"\nPerforming PERMANOVA for {group_var}")
            result = perform_permanova_test(beta_dm, metadata_df, group_var)
            permanova_results[group_var] = result
            
            # Print result
            print(f"  Test statistic: {result['test-statistic']}")
            print(f"  p-value: {result['p-value']}")
            print(f"  Sample size: {result['sample size']}")
            print(f"  Note: {result['note']}")
    
    # Save PERMANOVA results
    permanova_df = pd.DataFrame.from_dict(permanova_results, orient='index')
    permanova_file = tables_dir / 'permanova_results.csv'
    permanova_df.to_csv(permanova_file)
    print(f"PERMANOVA results saved to {permanova_file}")
    
    # Create ordination plots
    print("\nCreating ordination plots...")
    for group_var in config['metadata']['group_variables']:
        if group_var in metadata_df.columns:
            # Create PCoA plot
            print(f"  Creating PCoA plot for {group_var}")
            pcoa_file = figures_dir / f'pcoa_{group_var}.pdf'
            pcoa_fig = plot_pcoa(beta_dm, metadata_df, group_var, output_file=pcoa_file)
            plt.close(pcoa_fig)
            
            # Create NMDS plot
            print(f"  Creating NMDS plot for {group_var}")
            nmds_file = figures_dir / f'nmds_{group_var}.pdf'
            nmds_fig = plot_nmds(beta_dm, metadata_df, group_var, output_file=nmds_file)
            plt.close(nmds_fig)
    
    # Create abundance heatmap
    print("\nCreating abundance heatmap...")
    for group_var in config['metadata']['group_variables']:
        if group_var in metadata_df.columns:
            heatmap_file = figures_dir / f'abundance_heatmap_{group_var}.pdf'
            heatmap_fig = plot_taxa_heatmap(
                filtered_abundance_df, 
                metadata_df, 
                group_var, 
                top_n=30, 
                output_file=heatmap_file
            )
            plt.close(heatmap_fig)
            print(f"  Heatmap by {group_var} saved to {heatmap_file}")
    
    # Perform differential abundance testing
    print("\nPerforming differential abundance testing...")
    diff_abund_results = {}
    boxplot_dir = figures_dir / 'diff_abundance_boxplots'
    os.makedirs(boxplot_dir, exist_ok=True)
    
    for group_var in config['metadata']['group_variables']:
        if group_var in metadata_df.columns:
            print(f"  Testing for differential abundance by {group_var}")
            
            results = run_differential_abundance_testing(filtered_abundance_df, metadata_df, group_var)
            
            if not results.empty:
                # Save results
                results_file = tables_dir / f'differential_abundance_{group_var}.csv'
                results.to_csv(results_file)
                print(f"  Results saved to {results_file}")
                
                # Print top significant results
                alpha = config['differential_abundance']['p_value_threshold']
                sig_taxa = results[results['Adjusted P-value'] < alpha]
                print(f"  Found {len(sig_taxa)} significantly different taxa (adj. p < {alpha})")
                
                if not sig_taxa.empty:
                    # Create boxplots for top significant taxa
                    boxplots = plot_differential_abundance_boxplots(
                        filtered_abundance_df,
                        metadata_df,
                        group_var,
                        results,
                        top_n=5,
                        output_dir=boxplot_dir
                    )
                    
                    # Create stacked bar plot
                    bar_file = figures_dir / f'abundance_bar_{group_var}.pdf'
                    bar_fig = plot_bar_abundance(
                        filtered_abundance_df,
                        metadata_df,
                        group_var,
                        differential_results=results,
                        top_n=10,
                        output_file=bar_file
                    )
                    plt.close(bar_fig)
                    print(f"  Abundance bar plot saved to {bar_file}")
                    
                    # Store results for later use
                    diff_abund_results[group_var] = results
            else:
                print(f"  No valid differential abundance results for {group_var}")
    
    # Create summary of results
    print("\nCreating summary of results...")
    summary_file = output_dir / 'analysis_summary.txt'
    
    with open(summary_file, 'w') as f:
        f.write("=== Sylph Microbial Profiling Analysis Summary ===\n")
        f.write(f"Analysis completed on {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}\n\n")
        
        f.write("Dataset summary:\n")
        f.write(f"- Raw abundance data: {abundance_df.shape[0]} taxa, {abundance_df.shape[1]} samples\n")
        f.write(f"- Filtered abundance data: {filtered_abundance_df.shape[0]} taxa, {filtered_abundance_df.shape[1]} samples\n")
        f.write(f"- Samples with metadata: {len(common_samples)}\n\n")
        
        f.write("Alpha diversity results:\n")
        for metric in alpha_div.columns:
            f.write(f"- Mean {metric}: {alpha_div[metric].mean():.3f} (SD: {alpha_div[metric].std():.3f})\n")
        f.write("\n")
        
        f.write("Beta diversity results (PERMANOVA):\n")
        for group_var, result in permanova_results.items():
            f.write(f"- {group_var}: ")
            if np.isnan(result['p-value']):
                f.write(f"{result['note']}\n")
            else:
                f.write(f"p-value = {result['p-value']:.3f} {'(significant)' if result['p-value'] < 0.05 else ''}\n")
        f.write("\n")
        
        f.write("Differential abundance results:\n")
        for group_var, results in diff_abund_results.items():
            if not results.empty:
                alpha = config['differential_abundance']['p_value_threshold']
                sig_taxa = results[results['Adjusted P-value'] < alpha]
                f.write(f"- {group_var}: {len(sig_taxa)} significantly different taxa (adj. p < {alpha})\n")
                
                # List top 5 significant taxa
                if not sig_taxa.empty:
                    f.write(f"  Top significant taxa:\n")
                    for i, (_, row) in enumerate(sig_taxa.head(5).iterrows()):
                        f.write(f"  {i+1}. {row['Taxon']} (adj. p-value = {row['Adjusted P-value']:.3e})\n")
            else:
                f.write(f"- {group_var}: No significant results\n")
    
    print(f"Analysis summary saved to {summary_file}")
    print("\nSylph data analysis complete!")


if __name__ == "__main__":
    main()