#!/usr/bin/env python3
"""
Analyze longitudinal changes in the nasal microbiome over time.

This script:
1. Loads the combined abundance table and metadata
2. Tracks changes in species abundance across time points (Prior, Acute, Post)
3. Analyzes changes in diversity metrics over time
4. Creates visualizations of temporal microbiome changes
5. Identifies species with significant temporal patterns
6. Examines time-dependent effects of clinical variables

Usage:
    python scripts/04_longitudinal_analysis.py [--config CONFIG_FILE]
"""

import os
import sys
import argparse
import yaml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import stats

# Add project root to Python path
project_root = Path(__file__).resolve().parents[1]
sys.path.append(str(project_root))

# Add tools directory to Python path
tools_dir = project_root / 'tools'
sys.path.append(str(tools_dir))

# Import functions from metaphlan_tools
from metaphlan_tools import (
    load_metadata,
    calculate_alpha_diversity,
    plot_longitudinal_changes
)


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Analyze longitudinal microbiome changes')
    parser.add_argument('--config', type=str, default='config/analysis_parameters.yml',
                       help='Path to configuration file')
    return parser.parse_args()


def identify_temporal_patterns(abundance_df, metadata_df, time_var, subject_var, top_n=10):
    """
    Identify species with significant temporal patterns.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Species abundance DataFrame with species as index, samples as columns
    metadata_df : pandas.DataFrame
        Metadata DataFrame with samples as index
    time_var : str
        Metadata variable for time point
    subject_var : str
        Metadata variable for subject ID
    top_n : int
        Number of top species to return
        
    Returns:
    --------
    list
        List of species with significant temporal patterns
    dict
        Dictionary with statistical results for each species
    """
    # Get mean abundance across all samples for each species
    mean_abundance = abundance_df.mean(axis=1)
    
    # Sort species by mean abundance and get top 50 to test
    abundant_species = mean_abundance.sort_values(ascending=False).head(50).index.tolist()
    
    # Initialize results
    results = {}
    significant_species = []
    
    # For each abundant species, test for temporal patterns
    for species in abundant_species:
        species_data = abundance_df.loc[species]
        
        # Create DataFrame for analysis
        df = pd.DataFrame({
            'Abundance': species_data,
            time_var: metadata_df.loc[species_data.index, time_var],
            subject_var: metadata_df.loc[species_data.index, subject_var]
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
            
            results[species] = {
                'test': 'Friedman',
                'statistic': statistic,
                'p-value': p_value,
                'mean_abundance': mean_abundance[species],
                'n_subjects': pivot_df.shape[0]
            }
            
            # Check for significance
            if p_value < 0.05:
                significant_species.append(species)
        except Exception as e:
            print(f"Error testing {species}: {str(e)}")
            continue
    
    # Sort significant species by abundance
    significant_species.sort(key=lambda s: -mean_abundance[s])
    
    # Return top N species
    return significant_species[:min(top_n, len(significant_species))], results


def main():
    """Main function to analyze longitudinal changes."""
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
    longitudinal_dir = figures_dir / 'longitudinal_plots'
    
    # Create directories if they don't exist
    figures_dir.mkdir(exist_ok=True, parents=True)
    tables_dir.mkdir(exist_ok=True, parents=True)
    longitudinal_dir.mkdir(exist_ok=True, parents=True)
    
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
    
    # Define required variables
    time_var = config['metadata']['time_variable']
    subject_var = config['metadata']['subject_id_column']
    group_vars = config['metadata']['group_variables']
    
    # Check if required variables exist
    if time_var not in metadata_df.columns:
        print(f"Error: Time variable '{time_var}' not found in metadata")
        sys.exit(1)
    
    if subject_var not in metadata_df.columns:
        print(f"Error: Subject ID variable '{subject_var}' not found in metadata")
        sys.exit(1)
    
    # Calculate alpha diversity for longitudinal analysis
    print("\nCalculating alpha diversity for longitudinal analysis...")
    alpha_metrics = config['diversity']['alpha_metrics']
    alpha_df = calculate_alpha_diversity(abundance_df, metrics=alpha_metrics)
    
    # Analyze alpha diversity changes over time
    print(f"\nAnalyzing alpha diversity changes over {time_var}")
    
    # Join alpha diversity with metadata
    alpha_meta_df = alpha_df.copy()
    alpha_meta_df[time_var] = metadata_df.loc[alpha_meta_df.index, time_var]
    alpha_meta_df[subject_var] = metadata_df.loc[alpha_meta_df.index, subject_var]
    
    # Add group variables if present
    for group_var in group_vars:
        if group_var in metadata_df.columns:
            alpha_meta_df[group_var] = metadata_df.loc[alpha_meta_df.index, group_var]
    
    # Plot alpha diversity changes over time
    for metric in alpha_metrics:
        # Skip metrics that aren't in the alpha_df
        if metric not in alpha_df.columns:
            print(f"Warning: Metric '{metric}' not found in alpha diversity results. Skipping.")
            continue
            
        # Create summary boxplot
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Add individual subject trajectories as light lines
        for subject, data in alpha_meta_df.groupby(subject_var):
            data_sorted = data.sort_values(time_var)
            ax.plot(data_sorted[time_var], data_sorted[metric], 'o-', alpha=0.3, color='gray')
        
        # Add boxplot on top
        sns.boxplot(x=time_var, y=metric, data=alpha_meta_df, ax=ax)
        
        ax.set_title(f"{metric} Diversity Over {time_var}")
        ax.set_xlabel(time_var)
        ax.set_ylabel(f"{metric} Diversity")
        
        # Save figure
        fig_file = longitudinal_dir / f"alpha_{metric}_by_{time_var}.pdf"
        fig.savefig(fig_file, dpi=config['visualization']['figure_dpi'], bbox_inches='tight')
        plt.close(fig)
        print(f"  Alpha diversity time plot saved to {fig_file}")
        
        # If we have a grouping variable, create faceted plot
        for group_var in group_vars:
            if group_var in metadata_df.columns:
                # Create faceted plot
                fig, axes = plt.subplots(nrows=1, ncols=metadata_df[group_var].nunique(), 
                                         figsize=(4*metadata_df[group_var].nunique(), 6), 
                                         sharey=True)
                
                # Handle single subplot case
                if metadata_df[group_var].nunique() == 1:
                    axes = [axes]
                
                # Plot for each group
                for i, (group, group_data) in enumerate(alpha_meta_df.groupby(group_var)):
                    ax = axes[i]
                    
                    # Add individual subject trajectories
                    for subject, subject_data in group_data.groupby(subject_var):
                        subject_data_sorted = subject_data.sort_values(time_var)
                        ax.plot(subject_data_sorted[time_var], subject_data_sorted[metric], 
                               'o-', alpha=0.3, color='gray')
                    
                    # Add boxplot
                    sns.boxplot(x=time_var, y=metric, data=group_data, ax=ax)
                    
                    ax.set_title(f"{group_var}: {group}")
                    ax.set_xlabel(time_var)
                    
                    # Only add y-label on first subplot
                    if i == 0:
                        ax.set_ylabel(f"{metric} Diversity")
                    else:
                        ax.set_ylabel("")
                
                # Add overall title
                fig.suptitle(f"{metric} Diversity Over {time_var} by {group_var}", fontsize=14)
                plt.tight_layout(rect=[0, 0, 1, 0.95])  # Make room for suptitle
                
                # Save figure
                fig_file = longitudinal_dir / f"alpha_{metric}_by_{time_var}_and_{group_var}.pdf"
                fig.savefig(fig_file, dpi=config['visualization']['figure_dpi'], bbox_inches='tight')
                plt.close(fig)
                print(f"  Faceted time plot saved to {fig_file}")
    
    # Identify species with significant temporal patterns
    print("\nIdentifying species with significant temporal patterns...")
    top_n = config['longitudinal']['analyze_top_n_species']
    temporal_species, stats_results = identify_temporal_patterns(
        abundance_df, metadata_df, time_var, subject_var, top_n
    )
    
    # Save statistical results
    stats_df = pd.DataFrame.from_dict(stats_results, orient='index')
    stats_file = tables_dir / 'temporal_pattern_stats.csv'
    stats_df.to_csv(stats_file)
    print(f"Statistical results saved to {stats_file}")
    
    # If no species with significant temporal patterns were found, use top abundant species
    if not temporal_species:
        print("No species with significant temporal patterns found.")
        print("Using top species by abundance instead.")
        mean_abundance = abundance_df.mean(axis=1)
        temporal_species = mean_abundance.nlargest(top_n).index.tolist()
    else:
        print(f"Found {len(temporal_species)} species with significant temporal patterns:")
        for species in temporal_species:
            p_value = stats_results[species]['p-value']
            print(f"  - {species}: p-value = {p_value:.4f}")
    
    # Create longitudinal plots for each species
    print(f"\nCreating longitudinal plots for {len(temporal_species)} species...")
    for species in temporal_species:
        print(f"  Plotting {species}...")
        
        # Plot without grouping
        try:
            fig = plot_longitudinal_changes(
                abundance_df, 
                metadata_df, 
                species, 
                time_var, 
                subject_var
            )
            
            # Add p-value if available
            if species in stats_results:
                p_value = stats_results[species]['p-value']
                plt.annotate(f"p-value = {p_value:.4f}", xy=(0.02, 0.98), xycoords='axes fraction',
                           ha='left', va='top', fontsize=10,
                           bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=0.7))
            
            # Save figure
            fig_file = longitudinal_dir / f"species_{species.replace(' ', '_')}_by_{time_var}.pdf"
            fig.savefig(fig_file, dpi=config['visualization']['figure_dpi'], bbox_inches='tight')
            plt.close(fig)
            print(f"    Plot saved to {fig_file}")
        except Exception as e:
            print(f"    Error plotting {species}: {str(e)}")
        
        # Plot with grouping if available
        for group_var in group_vars:
            if group_var in metadata_df.columns:
                try:
                    fig = plot_longitudinal_changes(
                        abundance_df, 
                        metadata_df, 
                        species, 
                        time_var, 
                        subject_var, 
                        group_var
                    )
                    
                    # Save figure
                    fig_file = longitudinal_dir / f"species_{species.replace(' ', '_')}_by_{time_var}_and_{group_var}.pdf"
                    fig.savefig(fig_file, dpi=config['visualization']['figure_dpi'], bbox_inches='tight')
                    plt.close(fig)
                    print(f"    Plot with {group_var} saved to {fig_file}")
                except Exception as e:
                    print(f"    Error plotting {species} with {group_var}: {str(e)}")
    
    # Create visualization showing abundance changes over time for multiple species
    print("\nCreating combined visualization of top species changes...")
    try:
        # Get data for top species
        species_data = abundance_df.loc[temporal_species]
        
        # Calculate mean abundance at each time point
        time_points = sorted(metadata_df[time_var].unique())
        mean_by_time = {time_point: [] for time_point in time_points}
        
        for species in temporal_species:
            species_series = abundance_df.loc[species]
            for time_point in time_points:
                # Get samples for this time point
                time_samples = metadata_df[metadata_df[time_var] == time_point].index
                time_samples = [s for s in time_samples if s in species_series.index]
                
                # Calculate mean for this species at this time point
                if time_samples:
                    mean = species_series[time_samples].mean()
                    mean_by_time[time_point].append(mean)
                else:
                    mean_by_time[time_point].append(0)
        
        # Create a stacked bar chart
        fig, ax = plt.subplots(figsize=(10, 8))
        
        bottom = np.zeros(len(time_points))
        
        for i, species in enumerate(temporal_species):
            values = [mean_by_time[t][i] for t in time_points]
            ax.bar(time_points, values, bottom=bottom, label=species)
            bottom += values
        
        ax.set_title(f"Mean Abundance of Top Species Across {time_var}")
        ax.set_xlabel(time_var)
        ax.set_ylabel("Mean Relative Abundance (%)")
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        
        # Save figure
        fig_file = figures_dir / f"top_species_stacked_by_{time_var}.pdf"
        fig.savefig(fig_file, dpi=config['visualization']['figure_dpi'], bbox_inches='tight')
        plt.close(fig)
        print(f"Combined visualization saved to {fig_file}")
    except Exception as e:
        print(f"Error creating combined visualization: {str(e)}")
    
    print("\nLongitudinal analysis complete!")


if __name__ == "__main__":
    # Import seaborn here to avoid import error in function
    import seaborn as sns
    sns.set(style="whitegrid")
    
    main()