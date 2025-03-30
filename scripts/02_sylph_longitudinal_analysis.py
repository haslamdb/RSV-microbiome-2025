#!/usr/bin/env python3
"""
Analyze longitudinal changes in Sylph microbial profiling data.

This script:
1. Loads the processed Sylph abundance table and metadata
2. Tracks changes in taxa abundance across time points
3. Analyzes changes in diversity metrics over time
4. Creates visualizations of temporal microbiome changes
5. Identifies taxa with significant temporal patterns
6. Examines time-dependent effects of clinical variables

Usage:
    python scripts/02_sylph_longitudinal_analysis.py [--config CONFIG_FILE]
"""

import os
import sys
import argparse
import yaml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats
from statsmodels.stats.multitest import multipletests

# Add project root to Python path
project_root = Path(__file__).resolve().parents[1]
sys.path.append(str(project_root))

# Add tools/sylph_tools directory to Python path
tools_dir = project_root / 'tools' / 'sylph_tools'
sys.path.append(str(tools_dir))

# Try to import from sylph_tools if available
try:
    from sylph_tools import (
        load_metadata,
        calculate_alpha_diversity,
        filter_low_abundance
    )
    tools_available = True
except ImportError:
    tools_available = False
    print("Warning: sylph_tools module not available. Using built-in functions.")


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Analyze longitudinal Sylph data')
    parser.add_argument('--config', type=str, default='config/analysis_parameters.yml',
                       help='Path to configuration file')
    parser.add_argument('--abundance-file', type=str, 
                        default='results/sylph_analysis/filtered_abundance.csv',
                       help='Path to filtered Sylph abundance file')
    parser.add_argument('--output-dir', type=str, default='results/sylph_longitudinal',
                       help='Directory to save analysis results')
    parser.add_argument('--metadata-file', type=str, default=None,
                       help='Path to metadata file (override config)')
    return parser.parse_args()


def load_abundance_data(filepath):
    """
    Load abundance data from CSV file.
    
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


def identify_temporal_patterns(abundance_df, metadata_df, time_var, subject_var, top_n=10):
    """
    Identify taxa with significant temporal patterns.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Taxa abundance DataFrame with taxa as index, samples as columns
    metadata_df : pandas.DataFrame
        Metadata DataFrame with samples as index
    time_var : str
        Metadata variable for time point
    subject_var : str
        Metadata variable for subject ID
    top_n : int
        Number of top taxa to return
        
    Returns:
    --------
    list
        List of taxa with significant temporal patterns
    dict
        Dictionary with statistical results for each taxon
    """
    # Get mean abundance across all samples for each taxon
    mean_abundance = abundance_df.mean(axis=1)
    
    # Sort taxa by mean abundance and get top 50 to test
    abundant_taxa = mean_abundance.sort_values(ascending=False).head(50).index.tolist()
    
    # Initialize results
    results = {}
    significant_taxa = []
    
    # For each abundant taxon, test for temporal patterns
    for taxon in abundant_taxa:
        taxon_data = abundance_df.loc[taxon]
        
        # Create DataFrame for analysis
        df = pd.DataFrame({
            'Abundance': taxon_data,
            time_var: metadata_df.loc[taxon_data.index, time_var],
            subject_var: metadata_df.loc[taxon_data.index, subject_var]
        })
        
        # Check if we have enough time points and subjects
        if df[time_var].nunique() < 2 or df[subject_var].nunique() < 3:
            continue
        
        try:
            # Perform non-parametric test (Friedman test for repeated measures)
            # Group by subject and time, then reshape for Friedman test
            pivot_df = df.pivot_table(index=subject_var, columns=time_var, values='Abundance')
            
            # Drop subjects with missing time points
            pivot_df = pivot_df.dropna()
            
            # Need at least 3 subjects with all time points for Friedman test
            if pivot_df.shape[0] < 3:
                continue
                
            # Apply Friedman test
            statistic, p_value = stats.friedmanchisquare(*[pivot_df[col].values for col in pivot_df.columns])
            
            results[taxon] = {
                'test': 'Friedman',
                'statistic': statistic,
                'p-value': p_value,
                'mean_abundance': mean_abundance[taxon],
                'n_subjects': pivot_df.shape[0]
            }
            
            # Check for significance
            if p_value < 0.05:
                significant_taxa.append(taxon)
        except Exception as e:
            print(f"Error testing {taxon}: {str(e)}")
            continue
    
    # Sort significant taxa by abundance
    significant_taxa.sort(key=lambda s: -mean_abundance[s])
    
    # Return top N taxa
    return significant_taxa[:min(top_n, len(significant_taxa))], results


def plot_longitudinal_changes(abundance_df, metadata_df, taxon, time_var, subject_var, group_var=None):
    """
    Plot longitudinal changes in taxon abundance.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Taxa abundance DataFrame with taxa as index, samples as columns
    metadata_df : pandas.DataFrame
        Metadata DataFrame with samples as index
    taxon : str
        Taxon name to plot
    time_var : str
        Metadata variable for time point
    subject_var : str
        Metadata variable for subject ID
    group_var : str, optional
        Metadata variable for grouping subjects
        
    Returns:
    --------
    matplotlib.figure.Figure
        Longitudinal plot figure
    """
    # Get common samples
    common_samples = list(set(abundance_df.columns).intersection(set(metadata_df.index)))
    
    # Get taxon abundance
    if taxon in abundance_df.index:
        taxon_abundance = abundance_df.loc[taxon, common_samples]
    else:
        raise ValueError(f"Taxon '{taxon}' not found in abundance data")
    
    # Get metadata information
    time_info = metadata_df.loc[common_samples, time_var]
    subject_info = metadata_df.loc[common_samples, subject_var]
    
    # Create data frame for plotting
    plot_data = pd.DataFrame({
        'Abundance': taxon_abundance,
        time_var: time_info,
        subject_var: subject_info
    })
    
    # Add group information if provided
    if group_var is not None and group_var in metadata_df.columns:
        plot_data[group_var] = metadata_df.loc[common_samples, group_var]
    
    # Create figure
    if group_var is not None and group_var in plot_data.columns:
        # Create faceted plot by group
        groups = plot_data[group_var].unique()
        
        if len(groups) <= 1:
            # If only one group, fall back to non-grouped plot
            return plot_longitudinal_changes(abundance_df, metadata_df, taxon, time_var, subject_var)
        
        # Create subplot for each group
        fig, axes = plt.subplots(1, len(groups), figsize=(5*len(groups), 6), sharey=True)
        
        # Handle case with only one group (axes not in array)
        if len(groups) == 1:
            axes = [axes]
        
        # Plot each group
        for i, group in enumerate(groups):
            ax = axes[i]
            group_data = plot_data[plot_data[group_var] == group]
            
            # Plot individual subject trajectories
            for subject, subject_data in group_data.groupby(subject_var):
                # Sort by time
                subject_data = subject_data.sort_values(time_var)
                ax.plot(subject_data[time_var], subject_data['Abundance'], 'o-', alpha=0.3, label=subject)
            
            # Plot group mean
            mean_data = group_data.groupby(time_var)['Abundance'].mean().reset_index()
            ax.plot(mean_data[time_var], mean_data['Abundance'], 'o-', color='red', linewidth=3, label='Mean')
            
            # Set titles and labels
            ax.set_title(f"{group_var}: {group}")
            ax.set_xlabel(time_var)
            
            # Only add y-label on first subplot
            if i == 0:
                ax.set_ylabel('Abundance')
        
        # Add overall title
        plt.suptitle(f"{taxon} Abundance Over Time by {group_var}", fontsize=14)
        
    else:
        # Create simple plot without grouping
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Plot individual subject trajectories
        for subject, subject_data in plot_data.groupby(subject_var):
            # Sort by time
            subject_data = subject_data.sort_values(time_var)
            ax.plot(subject_data[time_var], subject_data['Abundance'], 'o-', alpha=0.3, label=subject)
        
        # Plot mean trajectory
        mean_data = plot_data.groupby(time_var)['Abundance'].mean().reset_index()
        ax.plot(mean_data[time_var], mean_data['Abundance'], 'o-', color='red', linewidth=3, label='Mean')
        
        # Set titles and labels
        ax.set_title(f"{taxon} Abundance Over Time")
        ax.set_xlabel(time_var)
        ax.set_ylabel('Abundance')
    
    # Adjust layout
    plt.tight_layout()
    
    return fig


def analyze_paired_timepoints(abundance_df, metadata_df, time_var, subject_var, time1, time2):
    """
    Analyze changes between paired time points using Wilcoxon signed-rank test.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Taxa abundance DataFrame with taxa as index, samples as columns
    metadata_df : pandas.DataFrame
        Metadata DataFrame with samples as index
    time_var : str
        Metadata variable for time point
    subject_var : str
        Metadata variable for subject ID
    time1 : str
        First time point
    time2 : str
        Second time point
        
    Returns:
    --------
    pandas.DataFrame
        DataFrame with paired test results
    """
    # Get common samples
    common_samples = list(set(abundance_df.columns).intersection(set(metadata_df.index)))
    
    # Filter metadata
    filtered_metadata = metadata_df.loc[common_samples]
    
    # Get subjects with both time points
    time1_subjects = set(filtered_metadata[filtered_metadata[time_var] == time1][subject_var])
    time2_subjects = set(filtered_metadata[filtered_metadata[time_var] == time2][subject_var])
    paired_subjects = time1_subjects.intersection(time2_subjects)
    
    if len(paired_subjects) < 3:
        print(f"Not enough subjects with both time points {time1} and {time2} (found {len(paired_subjects)})")
        return pd.DataFrame()
    
    print(f"Found {len(paired_subjects)} subjects with both time points {time1} and {time2}")
    
    # Initialize results
    results = []
    
    # Test each taxon
    for taxon in abundance_df.index:
        # Get samples for each time point
        time1_samples = filtered_metadata[(filtered_metadata[time_var] == time1) & 
                                        (filtered_metadata[subject_var].isin(paired_subjects))].index
        time2_samples = filtered_metadata[(filtered_metadata[time_var] == time2) & 
                                        (filtered_metadata[subject_var].isin(paired_subjects))].index
        
        # Get abundance values
        time1_abundance = abundance_df.loc[taxon, time1_samples].values
        time2_abundance = abundance_df.loc[taxon, time2_samples].values
        
        # Order samples by subject
        paired_values = []
        for subject in paired_subjects:
            subject_time1 = filtered_metadata[(filtered_metadata[time_var] == time1) & 
                                            (filtered_metadata[subject_var] == subject)].index
            subject_time2 = filtered_metadata[(filtered_metadata[time_var] == time2) & 
                                            (filtered_metadata[subject_var] == subject)].index
            
            if len(subject_time1) == 1 and len(subject_time2) == 1:
                value1 = abundance_df.loc[taxon, subject_time1[0]]
                value2 = abundance_df.loc[taxon, subject_time2[0]]
                paired_values.append((value1, value2))
        
        # Convert to arrays
        if paired_values:
            values1 = np.array([p[0] for p in paired_values])
            values2 = np.array([p[1] for p in paired_values])
            
            # Calculate mean and change
            mean1 = values1.mean()
            mean2 = values2.mean()
            mean_change = mean2 - mean1
            
            # Calculate fold change
            if mean1 > 0:
                fold_change = mean2 / mean1
            else:
                fold_change = np.nan
            
            # Perform Wilcoxon signed-rank test
            try:
                statistic, p_value = stats.wilcoxon(values1, values2)
            except Exception as e:
                statistic = np.nan
                p_value = 1.0
            
            results.append({
                'Taxon': taxon,
                'Time1': time1,
                'Time2': time2,
                'Mean at Time1': mean1,
                'Mean at Time2': mean2,
                'Mean Change': mean_change,
                'Fold Change': fold_change,
                'P-value': p_value,
                'Test Statistic': statistic,
                'N Pairs': len(paired_values)
            })
    
    # Create DataFrame from results
    results_df = pd.DataFrame(results)
    
    # Apply multiple testing correction if needed
    if not results_df.empty and len(results_df) > 1:
        try:
            results_df['Adjusted P-value'] = multipletests(
                results_df['P-value'], method='fdr_bh'
            )[1]
        except Exception as e:
            print(f"Error applying multiple testing correction: {str(e)}")
            results_df['Adjusted P-value'] = results_df['P-value']
    else:
        if not results_df.empty:
            results_df['Adjusted P-value'] = results_df['P-value']
    
    # Sort by adjusted p-value
    if not results_df.empty:
        results_df = results_df.sort_values('Adjusted P-value')
    
    return results_df


def plot_paired_changes(abundance_df, metadata_df, taxon, time_var, subject_var, time1, time2, output_file=None):
    """
    Create paired changes plot for a taxon between two time points.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Taxa abundance DataFrame with taxa as index, samples as columns
    metadata_df : pandas.DataFrame
        Metadata DataFrame with samples as index
    taxon : str
        Taxon name to plot
    time_var : str
        Metadata variable for time point
    subject_var : str
        Metadata variable for subject ID
    time1 : str
        First time point
    time2 : str
        Second time point
    output_file : str, optional
        Path to save plot
        
    Returns:
    --------
    matplotlib.figure.Figure
        Paired changes plot figure
    """
    # Get common samples
    common_samples = list(set(abundance_df.columns).intersection(set(metadata_df.index)))
    
    # Filter metadata
    filtered_metadata = metadata_df.loc[common_samples]
    
    # Get subjects with both time points
    time1_subjects = set(filtered_metadata[filtered_metadata[time_var] == time1][subject_var])
    time2_subjects = set(filtered_metadata[filtered_metadata[time_var] == time2][subject_var])
    paired_subjects = time1_subjects.intersection(time2_subjects)
    
    if len(paired_subjects) < 3:
        print(f"Not enough subjects with both time points {time1} and {time2} (found {len(paired_subjects)})")
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, f"Not enough subjects with both time points {time1} and {time2}",
               ha='center', va='center', fontsize=12)
        ax.set_title(f"{taxon}: {time1} vs {time2}")
        ax.axis('off')
        return fig
    
    # Get taxon abundance
    if taxon not in abundance_df.index:
        print(f"Taxon '{taxon}' not found in abundance data")
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, f"Taxon '{taxon}' not found in abundance data",
               ha='center', va='center', fontsize=12)
        ax.set_title(f"{taxon}: {time1} vs {time2}")
        ax.axis('off')
        return fig
    
    # Collect paired values
    paired_values = []
    labels = []
    
    for subject in paired_subjects:
        subject_time1 = filtered_metadata[(filtered_metadata[time_var] == time1) & 
                                        (filtered_metadata[subject_var] == subject)].index
        subject_time2 = filtered_metadata[(filtered_metadata[time_var] == time2) & 
                                        (filtered_metadata[subject_var] == subject)].index
        
        if len(subject_time1) == 1 and len(subject_time2) == 1:
            value1 = abundance_df.loc[taxon, subject_time1[0]]
            value2 = abundance_df.loc[taxon, subject_time2[0]]
            paired_values.append((value1, value2))
            labels.append(subject)
    
    # Create plot
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Convert to arrays
    values1 = np.array([p[0] for p in paired_values])
    values2 = np.array([p[1] for p in paired_values])
    
    # Plot paired lines
    for i in range(len(paired_values)):
        ax.plot([0, 1], [values1[i], values2[i]], 'o-', alpha=0.5)
    
    # Plot mean values
    mean1 = values1.mean()
    mean2 = values2.mean()
    ax.plot([0, 1], [mean1, mean2], 'o-', color='red', linewidth=3, label='Mean')
    
    # Add statistical test
    try:
        statistic, p_value = stats.wilcoxon(values1, values2)
        ax.text(0.5, 0.95, f"Wilcoxon signed-rank test\np-value = {p_value:.4f}",
               ha='center', va='top', transform=ax.transAxes,
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    except Exception as e:
        pass
    
    # Set labels and title
    ax.set_xticks([0, 1])
    ax.set_xticklabels([time1, time2])
    ax.set_ylabel('Abundance')
    ax.set_title(f"{taxon}: {time1} vs {time2}")
    
    # Add legend
    ax.legend()
    
    # Save figure if output file is provided
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
    
    return fig


def plot_diversity_changes(alpha_div, metadata_df, time_var, subject_var, group_var=None):
    """
    Plot changes in alpha diversity over time.
    
    Parameters:
    -----------
    alpha_div : pandas.DataFrame
        Alpha diversity DataFrame with samples as index
    metadata_df : pandas.DataFrame
        Metadata DataFrame with samples as index
    time_var : str
        Metadata variable for time point
    subject_var : str
        Metadata variable for subject ID
    group_var : str, optional
        Metadata variable for grouping subjects
        
    Returns:
    --------
    dict
        Dictionary of plot figures for each metric
    """
    # Get common samples
    common_samples = list(set(alpha_div.index).intersection(set(metadata_df.index)))
    
    # Filter to common samples
    alpha_subset = alpha_div.loc[common_samples]
    metadata_subset = metadata_df.loc[common_samples]
    
    # Add time and subject information
    alpha_meta = alpha_subset.copy()
    alpha_meta[time_var] = metadata_subset[time_var]
    alpha_meta[subject_var] = metadata_subset[subject_var]
    
    # Add group information if provided
    if group_var is not None and group_var in metadata_subset.columns:
        alpha_meta[group_var] = metadata_subset[group_var]
    
    # Create plots for each metric
    figures = {}
    
    for metric in alpha_subset.columns:
        if group_var is not None and group_var in alpha_meta.columns:
            # Create faceted plot by group
            groups = alpha_meta[group_var].unique()
            
            if len(groups) <= 1:
                # If only one group, fall back to non-grouped plot
                fig = _plot_diversity_time_series(alpha_meta, time_var, subject_var, metric)
                figures[metric] = fig
                continue
            
            # Create subplot for each group
            fig, axes = plt.subplots(1, len(groups), figsize=(5*len(groups), 6), sharey=True)
            
            # Handle case with only one group (axes not in array)
            if len(groups) == 1:
                axes = [axes]
            
            # Plot each group
            for i, group in enumerate(groups):
                ax = axes[i]
                group_data = alpha_meta[alpha_meta[group_var] == group]
                
                # Plot individual subject trajectories
                for subject, subject_data in group_data.groupby(subject_var):
                    # Sort by time
                    subject_data = subject_data.sort_values(time_var)
                    ax.plot(subject_data[time_var], subject_data[metric], 'o-', alpha=0.3)
                
                # Add boxplot
                sns.boxplot(x=time_var, y=metric, data=group_data, ax=ax)
                
                # Set titles and labels
                ax.set_title(f"{group_var}: {group}")
                ax.set_xlabel(time_var)
                
                # Only add y-label on first subplot
                if i == 0:
                    ax.set_ylabel(f'{metric} Diversity')
                
                # Rotate x-axis labels if needed
                plt.setp(ax.get_xticklabels(), rotation=45)
            
            # Add overall title
            plt.suptitle(f"{metric} Diversity Over Time by {group_var}", fontsize=14)
            
        else:
            # Create simple plot without grouping
            fig = _plot_diversity_time_series(alpha_meta, time_var, subject_var, metric)
        
        figures[metric] = fig
    
    return figures


def _plot_diversity_time_series(data, time_var, subject_var, metric):
    """Helper function to create a diversity time series plot."""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Plot individual subject trajectories
    for subject, subject_data in data.groupby(subject_var):
        # Sort by time
        subject_data = subject_data.sort_values(time_var)
        ax.plot(subject_data[time_var], subject_data[metric], 'o-', alpha=0.3)
    
    # Add boxplot
    sns.boxplot(x=time_var, y=metric, data=data, ax=ax)
    
    # Set titles and labels
    ax.set_title(f"{metric} Diversity Over Time")
    ax.set_xlabel(time_var)
    ax.set_ylabel(f'{metric} Diversity')
    
    # Rotate x-axis labels if needed
    plt.setp(ax.get_xticklabels(), rotation=45)
    
    # Adjust layout
    plt.tight_layout()
    
    return fig


def main():
    """Main function to analyze longitudinal Sylph data."""
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
        print("Using default parameters")
        config = {
            'metadata': {
                'filename': args.metadata_file or 'data/metadata.csv',
                'sample_id_column': 'SampleID',
                'time_variable': 'TimePoint',
                'subject_id_column': 'SubjectID',
                'group_variables': ['Group', 'Treatment']
            },
            'diversity': {
                'alpha_metrics': ['shannon', 'simpson', 'observed_otus']
            },
            'visualization': {
                'figure_dpi': 300
            },
            'longitudinal': {
                'analyze_top_n_species': 10
            }
        }
    
    # Set up output directories
    output_dir = Path(args.output_dir)
    figures_dir = output_dir / 'figures'
    tables_dir = output_dir / 'tables'
    paired_dir = figures_dir / 'paired_plots'
    
    # Create directories if they don't exist
    for dir_path in [output_dir, figures_dir, tables_dir, paired_dir]:
        os.makedirs(dir_path, exist_ok=True)
    
    # Load abundance data
    print(f"Loading Sylph abundance data from {args.abundance_file}")
    abundance_df = load_abundance_data(args.abundance_file)
    
    # Load metadata
    metadata_file = args.metadata_file or config['metadata']['filename']
    if tools_available:
        metadata_df = load_metadata(metadata_file, config['metadata']['sample_id_column'])
    else:
        try:
            metadata_df = pd.read_csv(metadata_file)
            metadata_df.set_index(config['metadata']['sample_id_column'], inplace=True)
        except Exception as e:
            print(f"Error loading metadata: {str(e)}")
            sys.exit(1)
    
    # Check for required metadata variables
    time_var = config['metadata']['time_variable']
    subject_var = config['metadata']['subject_id_column']
    
    if time_var not in metadata_df.columns:
        print(f"Error: Time variable '{time_var}' not found in metadata")
        sys.exit(1)
    
    if subject_var not in metadata_df.columns:
        print(f"Error: Subject ID variable '{subject_var}' not found in metadata")
        sys.exit(1)
    
    # Get sample overlap
    common_samples = list(set(abundance_df.columns).intersection(set(metadata_df.index)))
    print(f"Found {len(common_samples)} samples with both abundance data and metadata")
    
    # Filter to common samples
    filtered_abundance = abundance_df[common_samples]
    
    # Calculate alpha diversity
    print("\nCalculating alpha diversity...")
    if tools_available:
        alpha_div = calculate_alpha_diversity(filtered_abundance, config['diversity']['alpha_metrics'])
    else:
        from skbio.diversity import alpha_diversity
        alpha_div = pd.DataFrame(index=filtered_abundance.columns)
        
        for metric in config['diversity']['alpha_metrics']:
            # Convert to numpy array (samples as rows)
            counts = filtered_abundance.values.T
            sample_ids = filtered_abundance.columns
            
            # Calculate diversity
            if metric.lower() == 'observed_otus':
                # Count non-zero elements for each sample
                alpha_div[metric] = (filtered_abundance > 0).sum().values
            else:
                # Use scikit-bio
                alpha_values = alpha_diversity(metric.lower(), counts, sample_ids)
                alpha_div[metric] = alpha_values
    
    # Save alpha diversity results
    alpha_file = tables_dir / 'alpha_diversity.csv'
    alpha_div.to_csv(alpha_file)
    print(f"Alpha diversity results saved to {alpha_file}")
    
    # Plot alpha diversity changes over time
    print("\nAnalyzing alpha diversity changes over time...")
    diversity_figures = plot_diversity_changes(
        alpha_div, 
        metadata_df, 
        time_var, 
        subject_var, 
        group_var=config['metadata']['group_variables'][0] if config['metadata']['group_variables'] else None
    )
    
    # Save diversity figures
    for metric, fig in diversity_figures.items():
        fig_file = figures_dir / f'alpha_{metric}_over_time.pdf'
        fig.savefig(fig_file, dpi=config['visualization']['figure_dpi'], bbox_inches='tight')
        plt.close(fig)
        print(f"  {metric} diversity time plot saved to {fig_file}")
    
    # Identify taxa with significant temporal patterns
    print("\nIdentifying taxa with significant temporal patterns...")
    top_n = config['longitudinal']['analyze_top_n_species']
    temporal_taxa, stats_results = identify_temporal_patterns(
        filtered_abundance, metadata_df, time_var, subject_var, top_n
    )
    
    # Save statistical results
    stats_df = pd.DataFrame.from_dict(stats_results, orient='index')
    stats_file = tables_dir / 'temporal_pattern_stats.csv'
    stats_df.to_csv(stats_file)
    print(f"Temporal pattern statistics saved to {stats_file}")
    
    # If no taxa with significant temporal patterns were found, use top abundant taxa
    if not temporal_taxa:
        print("No taxa with significant temporal patterns found.")
        print("Using top taxa by abundance instead.")
        mean_abundance = filtered_abundance.mean(axis=1)
        temporal_taxa = mean_abundance.nlargest(top_n).index.tolist()
    else:
        print(f"Found {len(temporal_taxa)} taxa with significant temporal patterns:")
        for taxon in temporal_taxa:
            p_value = stats_results[taxon]['p-value']
            print(f"  - {taxon}: p-value = {p_value:.4f}")
    
    # Create longitudinal plots for each taxon
    print(f"\nCreating longitudinal plots for {len(temporal_taxa)} taxa...")
    for taxon in temporal_taxa:
        print(f"  Plotting {taxon}...")
        
        # Plot without grouping
        try:
            fig = plot_longitudinal_changes(
                filtered_abundance, 
                metadata_df, 
                taxon, 
                time_var, 
                subject_var
            )
            
            # Add p-value if available
            if taxon in stats_results:
                p_value = stats_results[taxon]['p-value']
                plt.annotate(f"p-value = {p_value:.4f}", xy=(0.02, 0.98), xycoords='axes fraction',
                           ha='left', va='top', fontsize=10,
                           bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=0.7))
            
            # Save figure
            safe_taxon = taxon.replace(' ', '_').replace('/', '_').replace(':', '_')
            fig_file = figures_dir / f"taxon_{safe_taxon}_over_time.pdf"
            fig.savefig(fig_file, dpi=config['visualization']['figure_dpi'], bbox_inches='tight')
            plt.close(fig)
            print(f"    Plot saved to {fig_file}")
        except Exception as e:
            print(f"    Error plotting {taxon}: {str(e)}")
        
        # Plot with grouping if available
        for group_var in config['metadata']['group_variables']:
            if group_var in metadata_df.columns:
                try:
                    fig = plot_longitudinal_changes(
                        filtered_abundance, 
                        metadata_df, 
                        taxon, 
                        time_var, 
                        subject_var, 
                        group_var
                    )
                    
                    # Save figure
                    safe_taxon = taxon.replace(' ', '_').replace('/', '_').replace(':', '_')
                    fig_file = figures_dir / f"taxon_{safe_taxon}_over_time_by_{group_var}.pdf"
                    fig.savefig(fig_file, dpi=config['visualization']['figure_dpi'], bbox_inches='tight')
                    plt.close(fig)
                    print(f"    Plot with {group_var} saved to {fig_file}")
                except Exception as e:
                    print(f"    Error plotting {taxon} with {group_var}: {str(e)}")
    
    # Analyze consecutive time points if more than one time point exists
    time_points = sorted(metadata_df[time_var].unique())
    
    if len(time_points) >= 2:
        print("\nAnalyzing changes between consecutive time points...")
        
        for i in range(len(time_points) - 1):
            time1 = time_points[i]
            time2 = time_points[i + 1]
            
            print(f"\nComparing {time1} vs {time2}...")
            
            # Perform paired time point analysis
            paired_results = analyze_paired_timepoints(
                filtered_abundance,
                metadata_df,
                time_var,
                subject_var,
                time1,
                time2
            )
            
            if not paired_results.empty:
                # Save results
                paired_file = tables_dir / f"paired_changes_{time1}_vs_{time2}.csv"
                paired_results.to_csv(paired_file)
                print(f"Paired analysis results saved to {paired_file}")
                
                # Print top significant results
                alpha = 0.05
                sig_taxa = paired_results[paired_results['Adjusted P-value'] < alpha]
                print(f"Found {len(sig_taxa)} significantly changing taxa (adj. p < {alpha})")
                
                if not sig_taxa.empty:
                    # Create paired plots for top significant taxa
                    for _, row in sig_taxa.head(5).iterrows():
                        taxon = row['Taxon']
                        
                        try:
                            # Create paired plot
                            safe_taxon = taxon.replace(' ', '_').replace('/', '_').replace(':', '_')
                            plot_file = paired_dir / f"paired_{safe_taxon}_{time1}_vs_{time2}.pdf"
                            
                            fig = plot_paired_changes(
                                filtered_abundance,
                                metadata_df,
                                taxon,
                                time_var,
                                subject_var,
                                time1,
                                time2,
                                output_file=plot_file
                            )
                            
                            plt.close(fig)
                            print(f"  Paired plot for {taxon} saved to {plot_file}")
                        except Exception as e:
                            print(f"  Error creating paired plot for {taxon}: {str(e)}")
    
    # Create visualization showing abundance changes over time for multiple taxa
    print("\nCreating combined visualization of top taxa changes...")
    try:
        # Get data for top taxa
        taxa_data = filtered_abundance.loc[temporal_taxa]
        
        # Calculate mean abundance at each time point
        time_points = sorted(metadata_df[time_var].unique())
        mean_by_time = {time_point: [] for time_point in time_points}
        
        for taxon in temporal_taxa:
            taxon_series = filtered_abundance.loc[taxon]
            for time_point in time_points:
                # Get samples for this time point
                time_samples = metadata_df[metadata_df[time_var] == time_point].index
                time_samples = [s for s in time_samples if s in taxon_series.index]
                
                # Calculate mean for this taxon at this time point
                if time_samples:
                    mean = taxon_series[time_samples].mean()
                    mean_by_time[time_point].append(mean)
                else:
                    mean_by_time[time_point].append(0)
        
        # Create a stacked bar chart
        fig, ax = plt.subplots(figsize=(10, 8))
        
        bottom = np.zeros(len(time_points))
        
        for i, taxon in enumerate(temporal_taxa):
            values = [mean_by_time[t][i] for t in time_points]
            ax.bar(time_points, values, bottom=bottom, label=taxon)
            bottom += values
        
        ax.set_title(f"Mean Abundance of Top Taxa Across {time_var}")
        ax.set_xlabel(time_var)
        ax.set_ylabel("Mean Relative Abundance (%)")
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        
        # Save figure
        fig_file = figures_dir / f"top_taxa_stacked_by_{time_var}.pdf"
        fig.savefig(fig_file, dpi=config['visualization']['figure_dpi'], bbox_inches='tight')
        plt.close(fig)
        print(f"Combined visualization saved to {fig_file}")
    except Exception as e:
        print(f"Error creating combined visualization: {str(e)}")
    
    # Create a summary of results
    print("\nCreating summary of results...")
    summary_file = output_dir / 'longitudinal_summary.txt'
    
    with open(summary_file, 'w') as f:
        f.write("=== Sylph Longitudinal Analysis Summary ===\n")
        f.write(f"Analysis completed on {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}\n\n")
        
        f.write("Dataset summary:\n")
        f.write(f"- Abundance data: {filtered_abundance.shape[0]} taxa, {filtered_abundance.shape[1]} samples\n")
        f.write(f"- Time points: {', '.join(str(t) for t in time_points)}\n")
        f.write(f"- Subjects: {metadata_df[subject_var].nunique()}\n\n")
        
        f.write("Alpha diversity changes over time:\n")
        for metric in diversity_figures.keys():
            try:
                # Calculate mean for each time point
                means = {}
                for time_point in time_points:
                    time_samples = metadata_df[metadata_df[time_var] == time_point].index
                    time_samples = list(set(time_samples).intersection(set(alpha_div.index)))
                    if time_samples:
                        means[time_point] = alpha_div.loc[time_samples, metric].mean()
                
                f.write(f"- {metric}: ")
                for time_point, mean in means.items():
                    f.write(f"{time_point}={mean:.3f} ")
                f.write("\n")
            except Exception as e:
                f.write(f"- {metric}: Error calculating means - {str(e)}\n")
        f.write("\n")
        
        f.write("Taxa with significant temporal patterns:\n")
        for i, taxon in enumerate(temporal_taxa):
            if taxon in stats_results:
                p_value = stats_results[taxon]['p-value']
                f.write(f"{i+1}. {taxon}: p-value = {p_value:.4f}\n")
            else:
                f.write(f"{i+1}. {taxon}\n")
        f.write("\n")
        
        if len(time_points) >= 2:
            f.write("Significant changes between consecutive time points:\n")
            for i in range(len(time_points) - 1):
                time1 = time_points[i]
                time2 = time_points[i + 1]
                
                paired_file = tables_dir / f"paired_changes_{time1}_vs_{time2}.csv"
                if os.path.exists(paired_file):
                    paired_results = pd.read_csv(paired_file)
                    sig_taxa = paired_results[paired_results['Adjusted P-value'] < 0.05]
                    
                    f.write(f"\n{time1} vs {time2}: {len(sig_taxa)} significant changes\n")
                    
                    if not sig_taxa.empty:
                        for j, (_, row) in enumerate(sig_taxa.head(5).iterrows()):
                            taxon = row['Taxon']
                            adj_p = row['Adjusted P-value']
                            fold_change = row['Fold Change']
                            direction = "increase" if row['Mean Change'] > 0 else "decrease"
                            
                            f.write(f"  {j+1}. {taxon}: {direction}, adj. p-value = {adj_p:.4e}, fold change = {fold_change:.2f}\n")
    
    print(f"Summary saved to {summary_file}")
    print("\nLongitudinal analysis complete!")