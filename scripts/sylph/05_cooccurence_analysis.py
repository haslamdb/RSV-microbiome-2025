#!/usr/bin/env python3
"""
Analyze Streptococcus pneumoniae and Haemophilus influenzae in RSV microbiome data.

This script:
1. Loads the previously generated abundance data
2. Extracts S. pneumoniae and H. influenzae abundance data
3. Analyzes their abundance patterns:
   - Across time points (Prior, Acute, Post)
   - Between severity groups at each time point
   - Between symptom groups at each time point
4. Examines co-occurrence patterns of both species
5. Generates visualizations regardless of statistical significance

Usage:
    python analyze_strep_haemophilus.py [--input-file FILTERED_ABUNDANCE_FILE] [--metadata METADATA_FILE]
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import mannwhitneyu, kruskal, spearmanr
from statsmodels.stats.multitest import multipletests
import traceback

# Set plotting style
plt.style.use('seaborn-v0_8-whitegrid')
sns.set(font_scale=1.1)

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Analyze S. pneumoniae and H. influenzae in RSV microbiome data')
    parser.add_argument('--input-file', type=str, default='results/rsv_analysis/abundance_tables/filtered_abundance.csv',
                       help='Path to filtered abundance table')
    parser.add_argument('--metadata', type=str, default='metadata.csv',
                       help='Path to metadata file')
    parser.add_argument('--output-dir', type=str, default='results/strep_haemophilus_analysis',
                       help='Directory to save analysis results')
    parser.add_argument('--sample-id-column', type=str, default='SampleID',
                       help='Column name for sample IDs in metadata')
    return parser.parse_args()

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
        
        # Check if the sample ID column exists
        if sample_id_column not in metadata_df.columns:
            raise ValueError(f"Sample ID column '{sample_id_column}' not found in metadata")
            
        # Set index and remove any duplicate sample IDs
        metadata_df = metadata_df.set_index(sample_id_column)
        if metadata_df.index.duplicated().any():
            print(f"Warning: Found {metadata_df.index.duplicated().sum()} duplicate sample IDs in metadata")
            metadata_df = metadata_df[~metadata_df.index.duplicated(keep='first')]
        
        # Convert categorical variables to string
        for col in metadata_df.columns:
            if metadata_df[col].dtype == 'object' or metadata_df[col].dtype.name == 'category':
                metadata_df[col] = metadata_df[col].astype(str)
        
        return metadata_df
    
    except Exception as e:
        print(f"Error loading metadata file: {str(e)}")
        return pd.DataFrame()

def find_species_in_data(abundance_df, species_names):
    """
    Find specified species in the abundance data.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Taxa abundance DataFrame with taxa as index, samples as columns
    species_names : list
        List of species names or patterns to search for
        
    Returns:
    --------
    dict
        Dictionary mapping species names to their exact index names in the DataFrame
    """
    species_indices = {}
    
    for species in species_names:
        matches = [idx for idx in abundance_df.index if species.lower() in idx.lower()]
        if matches:
            # Use the first match as the index
            species_indices[species] = matches[0]
            print(f"Found '{species}' as '{matches[0]}'")
            if len(matches) > 1:
                print(f"  Note: Multiple matches found: {matches}")
        else:
            print(f"Warning: Could not find '{species}' in the data")
    
    return species_indices

def perform_statistical_test(abundance_df, metadata_df, taxon, variable):
    """
    Perform statistical test to compare abundance across groups.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Taxa abundance DataFrame with taxa as index, samples as columns
    metadata_df : pandas.DataFrame
        Metadata DataFrame with samples as index
    taxon : str
        Taxon name to analyze
    variable : str
        Grouping variable to test (e.g., 'Timing', 'Severity', 'Symptoms')
        
    Returns:
    --------
    dict
        Dictionary with test results
    """
    # Get common samples
    common_samples = list(set(abundance_df.columns).intersection(set(metadata_df.index)))
    
    # Get abundance values for the taxon
    taxon_abundance = abundance_df.loc[taxon, common_samples]
    
    # Get group information
    groups = metadata_df.loc[common_samples, variable]
    unique_groups = groups.unique()
    
    # If there's only one group, no test can be performed
    if len(unique_groups) < 2:
        return {
            'taxon': taxon,
            'variable': variable,
            'test': 'None',
            'statistic': np.nan,
            'p-value': np.nan,
            'note': f'Only one group found in {variable}'
        }
    
    # For two groups, use Mann-Whitney U test
    if len(unique_groups) == 2:
        # Get values for each group
        group1_data = taxon_abundance[groups == unique_groups[0]]
        group2_data = taxon_abundance[groups == unique_groups[1]]
        
        # Calculate mean abundance in each group
        mean1 = group1_data.mean()
        mean2 = group2_data.mean()
        
        # Perform Mann-Whitney U test
        try:
            stat, p_value = mannwhitneyu(group1_data, group2_data, alternative='two-sided')
            
            return {
                'taxon': taxon,
                'variable': variable,
                'test': 'Mann-Whitney U',
                'statistic': stat,
                'p-value': p_value,
                'group1': unique_groups[0],
                'group2': unique_groups[1],
                'mean1': mean1,
                'mean2': mean2,
                'fold_change': (mean2 / mean1) if mean1 > 0 else np.nan
            }
        except Exception as e:
            return {
                'taxon': taxon,
                'variable': variable,
                'test': 'Mann-Whitney U',
                'statistic': np.nan,
                'p-value': np.nan,
                'note': f'Test failed: {str(e)}'
            }
    
    # For more than two groups, use Kruskal-Wallis test
    else:
        # Get values for each group
        group_data = [taxon_abundance[groups == group] for group in unique_groups]
        
        # Calculate mean abundance in each group
        group_means = {group: taxon_abundance[groups == group].mean() for group in unique_groups}
        
        # Perform Kruskal-Wallis test
        try:
            stat, p_value = kruskal(*group_data)
            
            result = {
                'taxon': taxon,
                'variable': variable,
                'test': 'Kruskal-Wallis',
                'statistic': stat,
                'p-value': p_value
            }
            
            # Add mean for each group
            for group, mean in group_means.items():
                result[f'mean_{group}'] = mean
            
            return result
        except Exception as e:
            return {
                'taxon': taxon,
                'variable': variable,
                'test': 'Kruskal-Wallis',
                'statistic': np.nan,
                'p-value': np.nan,
                'note': f'Test failed: {str(e)}'
            }

def analyze_species_by_variable(abundance_df, metadata_df, species_indices, variable):
    """
    Analyze species abundance patterns by a specific variable.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Taxa abundance DataFrame with taxa as index, samples as columns
    metadata_df : pandas.DataFrame
        Metadata DataFrame with samples as index
    species_indices : dict
        Dictionary mapping species names to their indices in the abundance DataFrame
    variable : str
        Metadata variable to analyze by
        
    Returns:
    --------
    pandas.DataFrame
        Results of statistical tests
    """
    results = []
    
    for species_name, taxon in species_indices.items():
        # Perform statistical test
        test_result = perform_statistical_test(abundance_df, metadata_df, taxon, variable)
        test_result['species_name'] = species_name
        results.append(test_result)
    
    # Convert to DataFrame
    results_df = pd.DataFrame(results)
    
    return results_df

def analyze_species_by_timepoint_and_variable(abundance_df, metadata_df, species_indices, time_var, group_var):
    """
    Analyze species abundance patterns by variable at each time point.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Taxa abundance DataFrame with taxa as index, samples as columns
    metadata_df : pandas.DataFrame
        Metadata DataFrame with samples as index
    species_indices : dict
        Dictionary mapping species names to their indices in the abundance DataFrame
    time_var : str
        Time variable column name (e.g., 'Timing')
    group_var : str
        Grouping variable column name (e.g., 'Severity', 'Symptoms')
        
    Returns:
    --------
    dict
        Dictionary with results for each time point
    """
    # Get time points
    time_points = sorted(metadata_df[time_var].unique())
    
    # Initialize results
    time_results = {}
    
    for time_point in time_points:
        print(f"Analyzing {time_point} by {group_var}")
        
        # Get samples for this time point
        time_samples = metadata_df[metadata_df[time_var] == time_point].index
        
        # Get common samples with abundance data
        common_samples = list(set(time_samples).intersection(set(abundance_df.columns)))
        
        if len(common_samples) < 3:
            print(f"  Not enough samples for time point {time_point} (found {len(common_samples)}), skipping")
            continue
            
        # Filter abundance and metadata to this time point
        time_abundance = abundance_df[common_samples]
        time_metadata = metadata_df.loc[common_samples]
        
        # Analyze by group variable
        results = []
        
        for species_name, taxon in species_indices.items():
            # Perform statistical test
            test_result = perform_statistical_test(time_abundance, time_metadata, taxon, group_var)
            test_result['species_name'] = species_name
            test_result['time_point'] = time_point
            results.append(test_result)
        
        # Convert to DataFrame
        results_df = pd.DataFrame(results)
        
        # Store results
        time_results[time_point] = results_df
    
    return time_results

def analyze_co_occurrence(abundance_df, metadata_df, species_indices):
    """
    Analyze co-occurrence patterns of the two species.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Taxa abundance DataFrame with taxa as index, samples as columns
    metadata_df : pandas.DataFrame
        Metadata DataFrame with samples as index
    species_indices : dict
        Dictionary mapping species names to their indices in the abundance DataFrame
        
    Returns:
    --------
    dict
        Dictionary with co-occurrence analysis results
    """
    # Check if we have both species
    if len(species_indices) < 2:
        print("Need both species for co-occurrence analysis")
        return {}
    
    # Get common samples
    common_samples = list(set(abundance_df.columns).intersection(set(metadata_df.index)))
    
    # Extract abundance for both species
    species_names = list(species_indices.keys())
    species1 = species_names[0]
    species2 = species_names[1]
    
    taxon1 = species_indices[species1]
    taxon2 = species_indices[species2]
    
    abundance1 = abundance_df.loc[taxon1, common_samples]
    abundance2 = abundance_df.loc[taxon2, common_samples]
    
    # Calculate presence/absence
    presence1 = (abundance1 > 0).astype(int)
    presence2 = (abundance2 > 0).astype(int)
    
    # Create contingency table
    contingency = pd.crosstab(presence1, presence2)
    
    # Calculate correlation coefficient
    corr, p_value = spearmanr(abundance1, abundance2)
    
    # Calculate co-occurrence rates by timepoint and clinical variables
    time_var = 'Timing'
    variables = ['Severity', 'Symptoms']
    
    co_occurrence_rates = {}
    
    # Overall co-occurrence
    co_occurrence_rates['overall'] = {
        'both_present': ((presence1 == 1) & (presence2 == 1)).mean() * 100,
        'species1_only': ((presence1 == 1) & (presence2 == 0)).mean() * 100,
        'species2_only': ((presence1 == 0) & (presence2 == 1)).mean() * 100,
        'none_present': ((presence1 == 0) & (presence2 == 0)).mean() * 100,
        'correlation': corr,
        'p_value': p_value
    }
    
    # By timepoint
    if time_var in metadata_df.columns:
        for time_point in metadata_df[time_var].unique():
            time_samples = metadata_df[metadata_df[time_var] == time_point].index
            time_samples = [s for s in time_samples if s in common_samples]
            
            if len(time_samples) < 3:
                continue
                
            time_presence1 = presence1.loc[time_samples]
            time_presence2 = presence2.loc[time_samples]
            
            time_abundance1 = abundance1.loc[time_samples]
            time_abundance2 = abundance2.loc[time_samples]
            
            # Calculate correlation
            try:
                time_corr, time_p = spearmanr(time_abundance1, time_abundance2)
            except:
                time_corr, time_p = np.nan, np.nan
            
            co_occurrence_rates[time_point] = {
                'both_present': ((time_presence1 == 1) & (time_presence2 == 1)).mean() * 100,
                'species1_only': ((time_presence1 == 1) & (time_presence2 == 0)).mean() * 100,
                'species2_only': ((time_presence1 == 0) & (time_presence2 == 1)).mean() * 100,
                'none_present': ((time_presence1 == 0) & (time_presence2 == 0)).mean() * 100,
                'correlation': time_corr,
                'p_value': time_p
            }
    
    # By clinical variables
    for variable in variables:
        if variable in metadata_df.columns:
            for group in metadata_df[variable].unique():
                group_samples = metadata_df[metadata_df[variable] == group].index
                group_samples = [s for s in group_samples if s in common_samples]
                
                if len(group_samples) < 3:
                    continue
                    
                group_presence1 = presence1.loc[group_samples]
                group_presence2 = presence2.loc[group_samples]
                
                group_abundance1 = abundance1.loc[group_samples]
                group_abundance2 = abundance2.loc[group_samples]
                
                # Calculate correlation
                try:
                    group_corr, group_p = spearmanr(group_abundance1, group_abundance2)
                except:
                    group_corr, group_p = np.nan, np.nan
                
                co_occurrence_rates[f"{variable}_{group}"] = {
                    'both_present': ((group_presence1 == 1) & (group_presence2 == 1)).mean() * 100,
                    'species1_only': ((group_presence1 == 1) & (group_presence2 == 0)).mean() * 100,
                    'species2_only': ((group_presence1 == 0) & (group_presence2 == 1)).mean() * 100,
                    'none_present': ((group_presence1 == 0) & (group_presence2 == 0)).mean() * 100,
                    'correlation': group_corr,
                    'p_value': group_p
                }
    
    return {
        'contingency': contingency,
        'correlation': corr,
        'p_value': p_value,
        'co_occurrence_rates': co_occurrence_rates
    }

def plot_species_by_time(abundance_df, metadata_df, species_indices, time_var, output_dir):
    """
    Create boxplots showing species abundance by time.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Taxa abundance DataFrame with taxa as index, samples as columns
    metadata_df : pandas.DataFrame
        Metadata DataFrame with samples as index
    species_indices : dict
        Dictionary mapping species names to their indices in the abundance DataFrame
    time_var : str
        Time variable column name (e.g., 'Timing')
    output_dir : str
        Directory to save plots
        
    Returns:
    --------
    None
    """
    # Get common samples
    common_samples = list(set(abundance_df.columns).intersection(set(metadata_df.index)))
    
    # Define time order
    time_order = ['Prior', 'Acute', 'Post']
    
    # Create figure and axes for each species
    for species_name, taxon in species_indices.items():
        # Get abundance values
        taxon_abundance = abundance_df.loc[taxon, common_samples]
        
        # Get time information
        time_info = metadata_df.loc[common_samples, time_var].copy()
        
        # Create DataFrame for plotting
        plot_data = pd.DataFrame({
            'Abundance': taxon_abundance,
            time_var: time_info
        })
        
        # Create categorical variable with correct order
        plot_data[time_var] = pd.Categorical(
            plot_data[time_var],
            categories=time_order,
            ordered=True
        )
        
        # Create figure
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Create boxplot
        sns.boxplot(x=time_var, y='Abundance', data=plot_data, ax=ax)
        
        # Add individual points
        sns.stripplot(x=time_var, y='Abundance', data=plot_data, 
                     color='black', size=4, alpha=0.5, ax=ax)
        
        # Perform Kruskal-Wallis test
        try:
            # Group data by time point
            groups = []
            for t in time_order:
                time_group = plot_data[plot_data[time_var] == t]['Abundance'].values
                if len(time_group) > 0:
                    groups.append(time_group)
            
            # Only perform test if we have at least 2 groups with data
            if len(groups) >= 2 and all(len(g) > 0 for g in groups):
                stat, p_value = kruskal(*groups)
                ax.text(0.5, 0.01, f'Kruskal-Wallis p = {p_value:.3f}', ha='center', va='bottom', 
                      transform=ax.transAxes,
                      bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.5'))
        except Exception as e:
            print(f"Error calculating p-value: {str(e)}")
        
        # Set title and labels
        ax.set_title(f'{species_name} Abundance Across Time Points', fontsize=14, fontweight='bold')
        ax.set_xlabel(time_var, fontsize=12)
        ax.set_ylabel('Abundance', fontsize=12)
        
        # Adjust layout
        plt.tight_layout()
        
        # Save figure
        output_file = os.path.join(output_dir, f"{species_name.replace(' ', '_')}_by_time.pdf")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Saved plot to {output_file}")
        plt.close(fig)

def plot_taxa_facet(abundance_df, metadata_df, taxon, time_var, group_var, output_file=None):
    """
    Create a faceted boxplot showing taxon abundance across time points by group.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Taxa abundance DataFrame with taxa as index, samples as columns
    metadata_df : pandas.DataFrame
        Metadata DataFrame with samples as index
    taxon : str
        Taxon name to plot
    time_var : str
        Time variable column name (e.g., 'Timing')
    group_var : str
        Grouping variable column name (e.g., 'Severity', 'Symptoms')
    output_file : str, optional
        Path to save the plot
        
    Returns:
    --------
    matplotlib.figure.Figure
        Boxplot figure
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    import numpy as np
    from scipy.stats import mannwhitneyu
    
    # Get common samples
    common_samples = list(set(abundance_df.columns).intersection(set(metadata_df.index)))
    
    # Get taxon abundance
    if taxon not in abundance_df.index:
        print(f"Error: Taxon '{taxon}' not found in abundance data")
        return None
    
    # Extract abundance data for the taxon
    taxon_abundance = abundance_df.loc[taxon, common_samples]
    
    # Get relevant metadata
    meta_subset = metadata_df.loc[common_samples, [time_var, group_var]].copy()
    
    # Define the correct order for time points (vertically from top to bottom)
    time_order = ['Prior', 'Acute', 'Post']
    
    # Create a DataFrame for plotting by merging abundance with metadata
    plot_data = pd.DataFrame()
    plot_data['SampleID'] = common_samples
    plot_data['Abundance'] = [taxon_abundance[sample] for sample in common_samples]
    plot_data[time_var] = [meta_subset.loc[sample, time_var] for sample in common_samples]
    plot_data[group_var] = [meta_subset.loc[sample, group_var] for sample in common_samples]
    
    # Create a categorical variable with the correct order
    plot_data[time_var] = pd.Categorical(
        plot_data[time_var],
        categories=time_order,
        ordered=True
    )
    
    # Create a vertical facet grid with time points as rows
    fig, axes = plt.subplots(3, 1, figsize=(10, 12), sharex=True)
    
    # Map time points to row indices
    time_to_row = {time: i for i, time in enumerate(time_order)}
    
    # Process each time point
    for time_point in time_order:
        # Get the corresponding row index
        row_idx = time_to_row[time_point]
        
        # Filter data for this time point
        time_data = plot_data[plot_data[time_var] == time_point]
        
        if time_data.empty:
            # If no data for this time point, add text to the axis
            axes[row_idx].text(0.5, 0.5, f'No data for {time_point}', 
                             ha='center', va='center', fontsize=12)
            axes[row_idx].set_title(time_point)
            continue
        
        # Create boxplot for this time point
        sns.boxplot(x=group_var, y='Abundance', data=time_data, ax=axes[row_idx])
        
        # Add individual points
        sns.stripplot(x=group_var, y='Abundance', data=time_data, 
                    color='black', size=4, alpha=0.5, ax=axes[row_idx])
        
        # Calculate p-value for this time point
        groups = time_data[group_var].unique()
        if len(groups) == 2:
            group1_data = time_data[time_data[group_var] == groups[0]]['Abundance']
            group2_data = time_data[time_data[group_var] == groups[1]]['Abundance']
            
            if len(group1_data) > 1 and len(group2_data) > 1:
                try:
                    _, p_value = mannwhitneyu(group1_data, group2_data, alternative='two-sided')
                    axes[row_idx].text(0.5, 0.01, f'p = {p_value:.3f}', ha='center', va='bottom', 
                                transform=axes[row_idx].transAxes,
                                bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.5'))
                except Exception as e:
                    print(f"Error calculating p-value: {str(e)}")
                    pass
        
        # Set title and labels
        axes[row_idx].set_title(f'{time_point}', fontsize=12, fontweight='bold')
        
        # Add y-label only to the middle subplot
        if row_idx == 1:
            axes[row_idx].set_ylabel('Abundance', fontsize=12)
        else:
            axes[row_idx].set_ylabel('')
            
        # Only add x-label to the bottom subplot
        if row_idx == len(time_order) - 1:
            axes[row_idx].set_xlabel(group_var, fontsize=12)
        else:
            axes[row_idx].set_xlabel('')
    
    # Add overall title
    plt.suptitle(f'{taxon} Abundance by {group_var} Across Time Points', 
                fontsize=14, fontweight='bold', y=0.98)
    
    # Adjust layout
    plt.tight_layout()
    fig.subplots_adjust(top=0.92)
    
    # Save figure if output file is provided
    if output_file is not None:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {output_file}")
    
    return fig

def plot_co_occurrence_heatmap(abundance_df, metadata_df, species_indices, output_dir):
    """
    Create a heatmap of co-occurrence patterns.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Taxa abundance DataFrame with taxa as index, samples as columns
    metadata_df : pandas.DataFrame
        Metadata DataFrame with samples as index
    species_indices : dict
        Dictionary mapping species names to their indices in the abundance DataFrame
    output_dir : str
        Directory to save plots
        
    Returns:
    --------
    None
    """
    # Check if we have both species
    if len(species_indices) < 2:
        print("Need both species for co-occurrence heatmap")
        return
    
    # Get common samples
    common_samples = list(set(abundance_df.columns).intersection(set(metadata_df.index)))
    
    # Extract abundance for both species
    species_names = list(species_indices.keys())
    species1 = species_names[0]
    species2 = species_names[1]
    
    taxon1 = species_indices[species1]
    taxon2 = species_indices[species2]
    
    abundance1 = abundance_df.loc[taxon1, common_samples]
    abundance2 = abundance_df.loc[taxon2, common_samples]
    
    # Create DataFrame for heatmap
    heatmap_data = pd.DataFrame({
        'Sample': common_samples,
        species1: abundance1.values,
        species2: abundance2.values
    })
    
    # Get timing and clinical variables
    time_var = 'Timing'
    clinical_vars = ['Severity', 'Symptoms']
    
    for var in [time_var] + clinical_vars:
        if var in metadata_df.columns:
            heatmap_data[var] = [metadata_df.loc[s, var] if s in metadata_df.index else 'Unknown' for s in common_samples]
    
    # Sort by variables
    sort_vars = [var for var in [time_var] + clinical_vars if var in heatmap_data.columns]
    if sort_vars:
        heatmap_data = heatmap_data.sort_values(sort_vars)
    
    # Convert to numeric for heatmap
    plot_data = heatmap_data[[species1, species2]].apply(lambda x: np.log1p(x) if (x > 0).any() else x)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Create heatmap
    ax = sns.heatmap(
        plot_data.T,
        cmap='YlGnBu',
        cbar_kws={'label': 'Log-transformed Abundance'},
        ax=ax
    )
    
    # Add row labels
    ax.set_yticklabels(species_names, rotation=0)
    
    # Add sample annotations if we have them
    if sort_vars:
        # Create a color list for each variable
        colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
        
        # Add colored bars for each variable
        for i, var in enumerate(sort_vars):
            unique_values = sorted(heatmap_data[var].unique())
            value_to_int = {val: j for j, val in enumerate(unique_values)}
            
            # Create color map
            color_map = {val: colors[j % len(colors)] for j, val in enumerate(unique_values)}
            
            # Create row colors
            row_colors = [color_map[val] for val in heatmap_data[var]]
            
            # Add colored bar above heatmap
            for j, sample in enumerate(heatmap_data.index):
                ax.add_patch(plt.Rectangle(
                    (j, plot_data.shape[1] + i*0.2),
                    1, 0.2,
                    color=color_map[heatmap_data.loc[sample, var]],
                    ec='none'
                ))
            
# Add legend
            from matplotlib.patches import Patch
            legend_elements = [Patch(facecolor=color_map[val], label=f'{var}: {val}') 
                              for val in unique_values]
            ax.legend(handles=legend_elements, loc='upper center', 
                     bbox_to_anchor=(0.5, 1.15), ncol=len(unique_values))
    
    # Set title
    plt.title('Co-occurrence of S. pneumoniae and H. influenzae', fontsize=14)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save figure
    output_file = os.path.join(output_dir, 'co_occurrence_heatmap.pdf')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved co-occurrence heatmap to {output_file}")
    plt.close(fig)

def plot_scatter_correlation(abundance_df, metadata_df, species_indices, output_dir):
    """
    Create a scatter plot showing correlation between the two species.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Taxa abundance DataFrame with taxa as index, samples as columns
    metadata_df : pandas.DataFrame
        Metadata DataFrame with samples as index
    species_indices : dict
        Dictionary mapping species names to their indices in the abundance DataFrame
    output_dir : str
        Directory to save plots
        
    Returns:
    --------
    None
    """
    # Check if we have both species
    if len(species_indices) < 2:
        print("Need both species for correlation plot")
        return
    
    # Get common samples
    common_samples = list(set(abundance_df.columns).intersection(set(metadata_df.index)))
    
    # Extract abundance for both species
    species_names = list(species_indices.keys())
    species1 = species_names[0]
    species2 = species_names[1]
    
    taxon1 = species_indices[species1]
    taxon2 = species_indices[species2]
    
    abundance1 = abundance_df.loc[taxon1, common_samples]
    abundance2 = abundance_df.loc[taxon2, common_samples]
    
    # Create DataFrame for scatter plot
    plot_data = pd.DataFrame({
        'Sample': common_samples,
        species1: abundance1.values,
        species2: abundance2.values
    })
    
    # Add timing information if available
    time_var = 'Timing'
    if time_var in metadata_df.columns:
        plot_data[time_var] = [metadata_df.loc[s, time_var] if s in metadata_df.index else 'Unknown' for s in common_samples]
    
    # Calculate correlation
    corr, p_value = spearmanr(abundance1, abundance2)
    
    # Create figure
    plt.figure(figsize=(10, 8))
    
    # Create scatter plot
    if time_var in plot_data.columns:
        # Create scatter plot colored by time point
        sns.scatterplot(
            data=plot_data, 
            x=species1, 
            y=species2, 
            hue=time_var,
            s=80,
            alpha=0.7
        )
    else:
        # Create simple scatter plot
        sns.scatterplot(
            data=plot_data, 
            x=species1, 
            y=species2,
            s=80,
            alpha=0.7
        )
    
    # Add regression line
    sns.regplot(
        data=plot_data, 
        x=species1, 
        y=species2,
        scatter=False,
        line_kws={'color': 'red', 'linestyle': '--'}
    )
    
    # Add correlation text
    plt.text(
        0.05, 0.95, 
        f'Spearman r = {corr:.3f}\np-value = {p_value:.3e}',
        transform=plt.gca().transAxes,
        bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.5')
    )
    
    # Set labels and title
    plt.xlabel(f'{species1} Abundance', fontsize=12)
    plt.ylabel(f'{species2} Abundance', fontsize=12)
    plt.title(f'Correlation between {species1} and {species2}', fontsize=14)
    
    # Save figure
    output_file = os.path.join(output_dir, 'species_correlation.pdf')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved correlation plot to {output_file}")
    plt.close()

def main():
    """Main function for analyzing S. pneumoniae and H. influenzae."""
    # Parse arguments
    args = parse_args()
    
    # Set up output directories
    output_dir = Path(args.output_dir)
    tables_dir = output_dir / "tables"
    figures_dir = output_dir / "figures"
    
    # Create output directories
    for dir_path in [output_dir, tables_dir, figures_dir]:
        os.makedirs(dir_path, exist_ok=True)
    
    print(f"Output will be saved to {output_dir}")
    
    # Load abundance data
    print("\n1. Loading abundance data")
    if os.path.exists(args.input_file):
        abundance_df = pd.read_csv(args.input_file, index_col=0)
        print(f"Loaded abundance data with {len(abundance_df.index)} taxa and {len(abundance_df.columns)} samples")
    else:
        print(f"Error: Input file not found at {args.input_file}")
        return
    
    # Load metadata
    print("\n2. Loading metadata")
    if os.path.exists(args.metadata):
        metadata_df = load_metadata(args.metadata, args.sample_id_column)
        
        if metadata_df.empty:
            print("Error: Could not load metadata")
            return
            
        print(f"Loaded metadata for {len(metadata_df)} samples")
        
        # Define required variables
        time_var = 'Timing'
        severity_var = 'Severity'
        symptoms_var = 'Symptoms'
        
        # Check if required variables exist
        missing_vars = []
        for var in [time_var, severity_var, symptoms_var]:
            if var not in metadata_df.columns:
                missing_vars.append(var)
                
        if missing_vars:
            print(f"Warning: The following variables are missing from metadata: {missing_vars}")
    else:
        print(f"Error: Metadata file not found at {args.metadata}")
        return
    
    # Find the two species in the data
    print("\n3. Finding target species in the data")
    target_species = ['Streptococcus pneumoniae', 'Haemophilus influenzae']
    species_indices = find_species_in_data(abundance_df, target_species)
    
    if not species_indices:
        print("Error: Could not find any target species in the data")
        return
    
    # Perform analyses
    print("\n4. Performing statistical analyses")
    
    # Analysis by timing
    print("\n4.1 Analyzing species abundance across timing")
    timing_results = analyze_species_by_variable(abundance_df, metadata_df, species_indices, time_var)
    
    if not timing_results.empty:
        # Save results
        timing_file = tables_dir / 'species_by_timing.csv'
        timing_results.to_csv(timing_file)
        print(f"Timing analysis results saved to {timing_file}")
        
        # Create plots
        print("Creating boxplots for species abundance by timing")
        plot_species_by_time(abundance_df, metadata_df, species_indices, time_var, figures_dir)
    
    # Analysis by severity
    if severity_var in metadata_df.columns:
        print(f"\n4.2 Analyzing species abundance by {severity_var}")
        severity_results = analyze_species_by_variable(abundance_df, metadata_df, species_indices, severity_var)
        
        if not severity_results.empty:
            # Save results
            severity_file = tables_dir / 'species_by_severity.csv'
            severity_results.to_csv(severity_file)
            print(f"Severity analysis results saved to {severity_file}")
            
            # Create plots for each species by severity
            for species_name, taxon in species_indices.items():
                try:
                    output_file = figures_dir / f"{species_name.replace(' ', '_')}_by_{severity_var}.pdf"
                    plot_taxa_facet(abundance_df, metadata_df, taxon, time_var, severity_var, output_file)
                except Exception as e:
                    print(f"Error creating plot for {species_name} by {severity_var}: {str(e)}")
    
    # Analysis by symptoms
    if symptoms_var in metadata_df.columns:
        print(f"\n4.3 Analyzing species abundance by {symptoms_var}")
        symptoms_results = analyze_species_by_variable(abundance_df, metadata_df, species_indices, symptoms_var)
        
        if not symptoms_results.empty:
            # Save results
            symptoms_file = tables_dir / 'species_by_symptoms.csv'
            symptoms_results.to_csv(symptoms_file)
            print(f"Symptoms analysis results saved to {symptoms_file}")
            
            # Create plots for each species by symptoms
            for species_name, taxon in species_indices.items():
                try:
                    output_file = figures_dir / f"{species_name.replace(' ', '_')}_by_{symptoms_var}.pdf"
                    plot_taxa_facet(abundance_df, metadata_df, taxon, time_var, symptoms_var, output_file)
                except Exception as e:
                    print(f"Error creating plot for {species_name} by {symptoms_var}: {str(e)}")
    
    # Co-occurrence analysis
    if len(species_indices) >= 2:
        print("\n4.4 Analyzing co-occurrence patterns")
        co_occurrence_results = analyze_co_occurrence(abundance_df, metadata_df, species_indices)
        
        if co_occurrence_results:
            # Save results
            result_file = tables_dir / 'co_occurrence_results.csv'
            pd.DataFrame(co_occurrence_results['co_occurrence_rates']).T.to_csv(result_file)
            print(f"Co-occurrence results saved to {result_file}")
            
            # Create visualizations
            print("Creating co-occurrence visualizations")
            plot_co_occurrence_heatmap(abundance_df, metadata_df, species_indices, figures_dir)
            plot_scatter_correlation(abundance_df, metadata_df, species_indices, figures_dir)
    
    # Analysis by timepoint and severity/symptoms
    print("\n4.5 Analyzing species abundance by severity and symptoms at each time point")
    
    # By severity
    if severity_var in metadata_df.columns:
        severity_time_results = analyze_species_by_timepoint_and_variable(
            abundance_df, metadata_df, species_indices, time_var, severity_var
        )
        
        for time_point, results_df in severity_time_results.items():
            # Save results
            result_file = tables_dir / f"{time_point}_{severity_var}_results.csv"
            results_df.to_csv(result_file)
            print(f"Results for {time_point} by {severity_var} saved to {result_file}")
    
    # By symptoms
    if symptoms_var in metadata_df.columns:
        symptoms_time_results = analyze_species_by_timepoint_and_variable(
            abundance_df, metadata_df, species_indices, time_var, symptoms_var
        )
        
        for time_point, results_df in symptoms_time_results.items():
            # Save results
            result_file = tables_dir / f"{time_point}_{symptoms_var}_results.csv"
            results_df.to_csv(result_file)
            print(f"Results for {time_point} by {symptoms_var} saved to {result_file}")
    
    print("\nAnalysis complete!")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"Error in main execution: {str(e)}")
        traceback.print_exc()