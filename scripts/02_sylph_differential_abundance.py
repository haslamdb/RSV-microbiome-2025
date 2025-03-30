#!/usr/bin/env python3
"""
Analyze RSV nasal microbiome data across timing points based on severity and symptoms.

This script:
1. Parses Sylph profiling output files (bacterial, viral, fungal)
2. Creates combined abundance tables
3. Performs differential abundance testing:
   - Between timing points (Prior, Acute, Post)
   - Between severity groups at each time point
   - Between symptom groups at each time point
4. Generates faceted boxplots for differentially abundant species

Usage:
    python 02_sylph_differential_abundance.py [--input-dir INPUT_DIR] [--output-dir OUTPUT_DIR]
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

# Set plotting style
plt.style.use('seaborn-v0_8-whitegrid')
sns.set(font_scale=1.1)

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Analyze RSV microbiome data')
    parser.add_argument('--input-dir', type=str, default='data/SylphProfiles',
                       help='Directory containing Sylph output files')
    parser.add_argument('--output-dir', type=str, default='results/rsv_analysis',
                       help='Directory to save analysis results')
    parser.add_argument('--metadata', type=str, default='metadata.csv',
                       help='Path to metadata file')
    parser.add_argument('--sample-id-column', type=str, default='SampleID',
                       help='Column name for sample IDs in metadata')
    parser.add_argument('--abundance-type', type=str, default='Taxonomic_abundance',
                       choices=['Taxonomic_abundance', 'Sequence_abundance'],
                       help='Which abundance measure to use')
    parser.add_argument('--min-ani', type=float, default=95.0,
                       help='Minimum Adjusted ANI to include in analysis')
    parser.add_argument('--min-coverage', type=float, default=0.05,
                       help='Minimum effective coverage to include in analysis')
    parser.add_argument('--min-prevalence', type=float, default=0.1,
                       help='Minimum prevalence to keep taxon')
    parser.add_argument('--min-abundance', type=float, default=0.01,
                       help='Minimum mean abundance to keep taxon')
    return parser.parse_args()

def extract_taxonomic_info(contig_name):
    """
    Extract taxonomic information from the contig name field in Sylph output.
    
    Parameters:
    -----------
    contig_name : str
        Contig name string from Sylph output
        
    Returns:
    --------
    dict
        Dictionary with taxonomic information (genus, species)
    """
    # Initialize with empty values
    taxonomy = {
        'genus': '',
        'species': '',
        'strain': '',
        'full_name': ''
    }
    
    if pd.isna(contig_name) or not contig_name:
        return taxonomy
    
    # Try to extract taxonomic information
    # First, keep everything after the first space (which should contain the organism name)
    parts = contig_name.split(' ', 1)
    if len(parts) < 2:
        return taxonomy
        
    organism_part = parts[1]
    
    # Look for typical binomial nomenclature pattern
    species_match = re.search(r'([A-Z][a-z]+)\s+([a-z]+)', organism_part)
    if species_match:
        taxonomy['genus'] = species_match.group(1)
        taxonomy['species'] = species_match.group(2)
        taxonomy['full_name'] = f"{taxonomy['genus']} {taxonomy['species']}"
        
        # Try to extract strain if present
        strain_match = re.search(r'strain\s+([^\s,]+)', organism_part)
        if strain_match:
            taxonomy['strain'] = strain_match.group(1)
    else:
        # Fall back to using the full text
        taxonomy['full_name'] = organism_part.split(',')[0].strip()
    
    return taxonomy

def parse_sylph_file(file_path, abundance_type='Taxonomic_abundance', min_ani=95.0, min_coverage=0.05):
    """
    Parse a Sylph output file and extract abundance data.
    
    Parameters:
    -----------
    file_path : str
        Path to the Sylph output file
    abundance_type : str
        Which abundance column to use ('Taxonomic_abundance' or 'Sequence_abundance')
    min_ani : float
        Minimum Adjusted ANI to include
    min_coverage : float
        Minimum effective coverage to include
        
    Returns:
    --------
    pandas.DataFrame
        DataFrame with species as index and abundance as values
    """
    # Check if file exists
    if not os.path.exists(file_path):
        print(f"Warning: File not found: {file_path}")
        return pd.DataFrame()
    
    try:
        # Read the Sylph file
        # The file doesn't have a clean header - need to handle specially
        with open(file_path, 'r') as f:
            header_line = f.readline().strip()
        
        # Clean up header line and split into column names
        header_columns = re.split(r'\s+', header_line.strip())
        
        # Read the file with the cleaned column names
        try:
            df = pd.read_csv(file_path, sep='\t', skiprows=1, names=header_columns)
            
            # Check if required columns exist
            required_cols = ['Adjusted_ANI', 'Eff_cov', 'Contig_name', abundance_type]
            missing_cols = [col for col in required_cols if col not in df.columns]
            if missing_cols:
                print(f"  Warning: Missing required columns: {missing_cols}")
                print(f"  Available columns: {df.columns.tolist()}")
                return pd.DataFrame()
            
        except Exception as e:
            print(f"  Error reading CSV: {str(e)}")
            return pd.DataFrame()
        
        # Filter based on ANI and coverage
        df = df[df['Adjusted_ANI'] >= min_ani]
        df = df[df['Eff_cov'] >= min_coverage]
        
        # Extract taxonomic information
        taxonomic_info = df['Contig_name'].apply(extract_taxonomic_info)
        tax_df = pd.DataFrame.from_records(taxonomic_info.tolist())
        
        # Join with the original dataframe
        df = pd.concat([df, tax_df], axis=1)
        
        # Use full species names for more accurate identification
        # Aggregate by species name and sum abundances
        abundance_data = {}
        
        # For each unique species, sum the abundance
        for species_name, group in df.groupby('full_name'):
            if species_name and not pd.isna(species_name):
                abundance_data[species_name] = group[abundance_type].sum()
        
        # Convert to DataFrame
        abundance_df = pd.DataFrame(list(abundance_data.items()), columns=['Taxon', 'abundance'])
        abundance_df.set_index('Taxon', inplace=True)
        
        return abundance_df
    
    except Exception as e:
        print(f"Error parsing {file_path}: {str(e)}")
        return pd.DataFrame()

def combine_sylph_samples(files, sample_ids=None, abundance_type='Taxonomic_abundance', 
                          min_ani=95.0, min_coverage=0.05):
    """
    Combine multiple Sylph output files into a single abundance table.
    
    Parameters:
    -----------
    files : list
        List of file paths to Sylph output files
    sample_ids : list, optional
        List of sample IDs corresponding to each file
    abundance_type : str
        Which abundance column to use
    min_ani : float
        Minimum Adjusted ANI to include
    min_coverage : float
        Minimum effective coverage to include
        
    Returns:
    --------
    pandas.DataFrame
        Combined abundance table with species as rows and samples as columns
    """
    dfs = []
    successful_files = []
    successful_sample_ids = []
    
    if sample_ids is None:
        # Extract sample IDs from file names
        sample_ids = []
        for f in files:
            # Extract sample ID from the Sylph output file name pattern
            match = re.search(r'(\d+-DNA|\w+)_profiled_', os.path.basename(f))
            if match:
                sample_ids.append(match.group(1))
            else:
                # Fallback to file basename without extension
                sample_ids.append(os.path.splitext(os.path.basename(f))[0])
    
    if len(files) != len(sample_ids):
        raise ValueError("Number of files must match number of sample IDs")
    
    for i, file_path in enumerate(files):
        try:
            # Parse the Sylph file
            print(f"Processing file {i+1}/{len(files)}: {os.path.basename(file_path)}")
            df = parse_sylph_file(file_path, abundance_type, min_ani, min_coverage)
            
            if df.empty:
                print(f"Warning: No data extracted from {file_path}")
                continue
                
            # Set column name to sample ID
            df.columns = [sample_ids[i]]
            
            # Check for duplicate indices and handle them
            if df.index.duplicated().any():
                print(f"Warning: Found duplicate taxa in {sample_ids[i]}, keeping first occurrence")
                df = df[~df.index.duplicated(keep='first')]
            
            dfs.append(df)
            successful_files.append(file_path)
            successful_sample_ids.append(sample_ids[i])
            print(f"Successfully processed {os.path.basename(file_path)}")
        except Exception as e:
            print(f"Error processing {file_path}: {str(e)}")
            continue
    
    if not dfs:
        # Return an empty DataFrame with a warning
        print("\nWarning: No valid data frames to combine. Returning empty DataFrame.")
        return pd.DataFrame()
    
    print(f"\nSuccessfully processed {len(dfs)}/{len(files)} files")
    
    # Combine along columns (each sample is a column)
    combined_df = pd.concat(dfs, axis=1)
    
    # Fill missing values with zeros
    combined_df = combined_df.fillna(0)
    
    return combined_df

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

def filter_low_abundance(abundance_df, min_prevalence=0.1, min_abundance=0.01):
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

def differential_abundance_analysis(abundance_df, metadata_df, group_var, adjusted_p_threshold=0.05):
    """
    Identify differentially abundant taxa between groups.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Taxa abundance DataFrame with taxa as index, samples as columns
    metadata_df : pandas.DataFrame
        Metadata DataFrame with samples as index
    group_var : str
        Grouping variable in metadata (e.g., 'Severity', 'Symptoms')
    adjusted_p_threshold : float
        P-value threshold after multiple testing correction
        
    Returns:
    --------
    pandas.DataFrame
        Results of differential abundance testing
    """
    # Get common samples
    common_samples = list(set(abundance_df.columns).intersection(set(metadata_df.index)))
    
    if len(common_samples) < 5:
        print(f"Error: Not enough samples for differential abundance testing (found {len(common_samples)})")
        return pd.DataFrame()
    
    # Filter to common samples
    filtered_abundance = abundance_df[common_samples]
    
    # Get groups
    groups = metadata_df.loc[common_samples, group_var]
    unique_groups = groups.unique()
    
    if len(unique_groups) < 2:
        print(f"Error: Need at least 2 groups for differential abundance testing (found {len(unique_groups)})")
        return pd.DataFrame()
    
    # Determine which method to use based on number of groups
    if len(unique_groups) == 2:
        # Use Mann-Whitney U test for two groups
        results_df = _two_group_differential_testing(filtered_abundance, groups, unique_groups)
    else:
        # Use Kruskal-Wallis test for multiple groups
        results_df = _multi_group_differential_testing(filtered_abundance, groups, unique_groups)
    
    # Add variable name for reference
    results_df['Variable'] = group_var
    
    # Filter to significant results
    significant_df = results_df[results_df['Adjusted P-value'] < adjusted_p_threshold]
    
    print(f"Found {len(significant_df)} significantly different taxa (adj. p < {adjusted_p_threshold})")
    
    return results_df

def plot_timing_boxplot(abundance_df, metadata_df, taxon, time_var, output_file=None):
    """
    Create a boxplot showing taxon abundance changes across time points.
    
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
    output_file : str, optional
        Path to save the plot
        
    Returns:
    --------
    matplotlib.figure.Figure
        Boxplot figure
    """
    # Get common samples
    common_samples = list(set(abundance_df.columns).intersection(set(metadata_df.index)))
    
    # Get taxon abundance
    if taxon not in abundance_df.index:
        print(f"Error: Taxon '{taxon}' not found in abundance data")
        return None
    
    taxon_abundance = abundance_df.loc[taxon, common_samples]
    
    # Get time information
    time_info = metadata_df.loc[common_samples, time_var].copy()
    
    # Create a DataFrame for plotting
    plot_data = pd.DataFrame({
        'Abundance': taxon_abundance,
        time_var: time_info
    })
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))

    # Define the correct order for time points
    time_order = ['Prior', 'Acute', 'Post']
    
    # Ensure plot_data is initialized
    plot_data = metadata_df.copy()


    # Create a categorical variable with the correct order
    plot_data[time_var] = pd.Categorical(
        plot_data[time_var],
        categories=time_order,
        ordered=True
    )
    
    
    # Create boxplot
    sns.boxplot(x=time_var, y='Abundance', data=plot_data, ax=ax)
    
    # Add individual points
    sns.stripplot(x=time_var, y='Abundance', data=plot_data, 
                color='black', size=4, alpha=0.5, ax=ax)
    
    # Perform Kruskal-Wallis test
    try:
        groups = [plot_data[plot_data[time_var] == t]['Abundance'].values for t in plot_data[time_var].unique()]
        if all(len(g) > 0 for g in groups):
            stat, p_value = kruskal(*groups)
            ax.text(0.5, 0.01, f'Kruskal-Wallis p = {p_value:.3f}', ha='center', va='bottom', 
                  transform=ax.transAxes,
                  bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.5'))
    except:
        pass
    
    # Set title and labels
    ax.set_title(f'{taxon} Abundance Across Time Points')
    ax.set_xlabel(time_var)
    ax.set_ylabel('Abundance')
    
    # Adjust layout
    plt.tight_layout()
    
    # Save figure if output file is provided
    if output_file is not None:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {output_file}")
    
    return fig

def _two_group_differential_testing(abundance_df, groups, unique_groups):
    """Helper function for differential testing with two groups using Mann-Whitney U test."""
    # Initialize results
    results = []
    
    # Get sample indices for each group
    group1_samples = groups[groups == unique_groups[0]].index
    group2_samples = groups[groups == unique_groups[1]].index
    
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
            stat = np.nan
        
        results.append({
            'Taxon': taxon,
            'P-value': p_value,
            'Group1': unique_groups[0],
            'Group2': unique_groups[1],
            'Mean in Group1': mean1,
            'Mean in Group2': mean2,
            'Log2 Fold Change': fold_change,
            'Fold Change': 2**fold_change,
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
        if not results_df.empty:
            results_df['Adjusted P-value'] = results_df['P-value']
    
    # Sort by adjusted p-value
    if not results_df.empty:
        results_df = results_df.sort_values('Adjusted P-value')
    
    return results_df

def _multi_group_differential_testing(abundance_df, groups, unique_groups):
    """Helper function for differential testing with multiple groups using Kruskal-Wallis test."""
    # Initialize results
    results = []
    
    # Loop through each taxon
    for taxon in abundance_df.index:
        # Initialize group values and means
        group_values = {}
        group_means = {}
        
        # Get values and calculate means for each group
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
            stat = np.nan
        
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
            'Test': 'Kruskal-Wallis',
            'Max Log2 Fold Change': max_fold_change,
        }
        
        # Add mean for each group
        for group in unique_groups:
            result[f'Mean in {group}'] = group_means[group]
        
        # Add fold change information
        if max_group_pair:
            result['Max Fold Change Groups'] = f'{max_group_pair[0]} vs {max_group_pair[1]}'
            result['Fold Change'] = 2**max_fold_change
        
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
        if not results_df.empty:
            results_df['Adjusted P-value'] = results_df['P-value']
    
    # Sort by adjusted p-value
    if not results_df.empty:
        results_df = results_df.sort_values('Adjusted P-value')
    
    return results_df

def plot_taxa_boxplot(abundance_df, metadata_df, taxon, time_var, group_var, output_file=None):
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
    # Get common samples
    common_samples = list(set(abundance_df.columns).intersection(set(metadata_df.index)))
    
    # Get taxon abundance
    if taxon not in abundance_df.index:
        print(f"Error: Taxon '{taxon}' not found in abundance data")
        return None
    
    taxon_abundance = abundance_df.loc[taxon, common_samples]
    
    # Get metadata information
    meta_subset = metadata_df.loc[common_samples, [time_var, group_var]].copy()

    # Define the correct order for time points
    time_order = ['Prior', 'Acute', 'Post']
    
    # Create a categorical variable with the correct order
    plot_data[time_var] = pd.Categorical(
        plot_data[time_var],
        categories=time_order,
        ordered=True
    )
    
    
    # Create a DataFrame for plotting
    plot_data = pd.DataFrame({
        'Abundance': taxon_abundance,
        time_var: meta_subset[time_var],
        group_var: meta_subset[group_var]
    })
    
    # Create figure with appropriate size for facets
    time_points = plot_data[time_var].unique()
    fig, axes = plt.subplots(1, len(time_points), figsize=(4*len(time_points), 5), sharey=True)
    
    # Handle the case with only one time point
    if len(time_points) == 1:
        axes = [axes]
    
    # Create boxplot for each time point
    for i, time_point in enumerate(sorted(time_points)):
        time_data = plot_data[plot_data[time_var] == time_point]
        
        # Create boxplot for this time point
        sns.boxplot(x=group_var, y='Abundance', data=time_data, ax=axes[i])
        
        # Add individual points
        sns.stripplot(x=group_var, y='Abundance', data=time_data, 
                     color='black', size=4, alpha=0.5, ax=axes[i])
        
        # Calculate p-value for this time point
        groups = time_data[group_var].unique()
        if len(groups) == 2:
            group1_data = time_data[time_data[group_var] == groups[0]]['Abundance']
            group2_data = time_data[time_data[group_var] == groups[1]]['Abundance']
            
            if len(group1_data) > 1 and len(group2_data) > 1:
                try:
                    _, p_value = mannwhitneyu(group1_data, group2_data, alternative='two-sided')
                    axes[i].text(0.5, 0.01, f'p = {p_value:.3f}', ha='center', va='bottom', 
                                transform=axes[i].transAxes,
                                bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.5'))
                except:
                    pass
        
        # Set title and labels
        axes[i].set_title(f'{time_point}')
        axes[i].set_xlabel(group_var)
        
        # Only add y-label to the first subplot
        if i == 0:
            axes[i].set_ylabel('Abundance')
        else:
            axes[i].set_ylabel('')
    
    # Add overall title
    plt.suptitle(f'{taxon} Abundance by {group_var} Across Time Points', fontsize=14, y=1.05)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save figure if output file is provided
    if output_file is not None:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {output_file}")
    
    return fig

def analyze_by_timepoint(abundance_df, metadata_df, time_var, group_vars, adjusted_p_threshold=0.05):
    """
    Perform differential abundance analysis at each time point for each grouping variable.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Taxa abundance DataFrame with taxa as index, samples as columns
    metadata_df : pandas.DataFrame
        Metadata DataFrame with samples as index
    time_var : str
        Time variable column name (e.g., 'Timing')
    group_vars : list
        List of grouping variables to analyze (e.g., ['Severity', 'Symptoms'])
    adjusted_p_threshold : float
        P-value threshold after multiple testing correction
        
    Returns:
    --------
    dict
        Dictionary of differential abundance results by time point and group variable
    """
    # Check if time variable exists in metadata
    if time_var not in metadata_df.columns:
        print(f"Error: Time variable '{time_var}' not found in metadata")
        return {}
    
        # Define the correct order for time points
    time_order = ['Prior', 'Acute', 'Post']
    
    # Create a categorical variable with the correct order
    metadata_df[time_var] = pd.Categorical(
        metadata_df[time_var],
        categories=time_order,
        ordered=True
    )
    
    
    # Get all time points
    time_points = metadata_df[time_var].unique()
    
    # Initialize results dictionary
    results = {}
    
    # For each time point, analyze by each group variable
    for time_point in time_points:
        print(f"\nAnalyzing time point: {time_point}")
        
        # Get samples for this time point
        time_samples = metadata_df[metadata_df[time_var] == time_point].index
        
        # Get common samples with abundance data
        common_samples = list(set(time_samples).intersection(set(abundance_df.columns)))
        
        if len(common_samples) < 5:
            print(f"  Not enough samples for time point {time_point} (found {len(common_samples)}), skipping")
            continue
            
        # Filter abundance and metadata to this time point
        time_abundance = abundance_df[common_samples]
        time_metadata = metadata_df.loc[common_samples]
        
        # Initialize results for this time point
        results[time_point] = {}
        
        # Analyze by each group variable
        for group_var in group_vars:
            if group_var not in metadata_df.columns:
                print(f"  Warning: Group variable '{group_var}' not found in metadata, skipping")
                continue
                
            print(f"  Analyzing by {group_var}")
            
            # Perform differential abundance analysis
            diff_results = differential_abundance_analysis(
                time_abundance, 
                time_metadata, 
                group_var,
                adjusted_p_threshold
            )
            
            # Store results
            results[time_point][group_var] = diff_results
            
            # Print top significant results
            sig_results = diff_results[diff_results['Adjusted P-value'] < adjusted_p_threshold]
            
            if not sig_results.empty:
                print(f"    Top significant taxa:")
                for i, (_, row) in enumerate(sig_results.head(5).iterrows()):
                    if i < 5:  # Only show top 5
                        print(f"      {i+1}. {row['Taxon']}: adj. p-value = {row['Adjusted P-value']:.4f}")
    
    return results

def analyze_by_timing(abundance_df, metadata_df, time_var, adjusted_p_threshold=0.05):
    """
    Perform analysis to identify taxa that change significantly across time points.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Taxa abundance DataFrame with taxa as index, samples as columns
    metadata_df : pandas.DataFrame
        Metadata DataFrame with samples as index
    time_var : str
        Time variable column name (e.g., 'Timing')
    adjusted_p_threshold : float
        P-value threshold after multiple testing correction
        
    Returns:
    --------
    pandas.DataFrame
        Differential abundance results across time points
    """
    # Check if time variable exists in metadata
    if time_var not in metadata_df.columns:
        print(f"Error: Time variable '{time_var}' not found in metadata")
        return pd.DataFrame()
    
        # Define the correct order for time points
    time_order = ['Prior', 'Acute', 'Post']
    
    # Create a categorical variable with the correct order
    metadata_df[time_var] = pd.Categorical(
        metadata_df[time_var],
        categories=time_order,
        ordered=True
    )
    
    
    # Get all time points
    time_points = sorted(metadata_df[time_var].unique())
    
    if len(time_points) < 2:
        print(f"Error: Need at least 2 time points for temporal analysis")
        return pd.DataFrame()
    
    # Initialize results
    results = []
    
    # For each taxon, test for differences across time points
    for taxon in abundance_df.index:
        print(f"  Testing taxon: {taxon}")
        
        # Create DataFrame for analysis
        taxon_data = []
        
        for time_point in time_points:
            # Get samples for this time point
            time_samples = metadata_df[metadata_df[time_var] == time_point].index
            
            # Get common samples with abundance data
            common_samples = list(set(time_samples).intersection(set(abundance_df.columns)))
            
            if len(common_samples) < 3:
                print(f"    Not enough samples for time point {time_point} (found {len(common_samples)}), skipping taxon")
                continue
                
            # Get abundance values for this time point
            time_abundance = abundance_df.loc[taxon, common_samples]
            
            # Add to data
            for sample, abundance in time_abundance.items():
                taxon_data.append({
                    'Sample': sample,
                    'Abundance': abundance,
                    time_var: time_point
                })
        
        # Convert to DataFrame
        df = pd.DataFrame(taxon_data)
        
        # Skip if not enough data
        if len(df) < 5:
            continue
            
        # Perform Kruskal-Wallis test across time points
        try:
            # Group data by time point
            groups = [df[df[time_var] == time_point]['Abundance'].values for time_point in time_points
                     if len(df[df[time_var] == time_point]) > 0]
            
            # Skip if any group is empty
            if any(len(g) == 0 for g in groups):
                continue
                
            # Perform Kruskal-Wallis test
            stat, p_value = kruskal(*groups)
            
            # Calculate mean abundance for each time point
            time_means = {}
            for time_point in time_points:
                time_group = df[df[time_var] == time_point]
                if not time_group.empty:
                    time_means[time_point] = time_group['Abundance'].mean()
                else:
                    time_means[time_point] = np.nan
            
            # Calculate maximum fold change between any two time points
            max_fold_change = 0
            max_time_pair = None
            
            for i, time1 in enumerate(time_points):
                for time2 in time_points[i+1:]:
                    # Skip if either mean is NaN
                    if np.isnan(time_means[time1]) or np.isnan(time_means[time2]):
                        continue
                        
                    # Add small pseudocount to avoid division by zero
                    pseudocount = 1e-5
                    fold_change = np.abs(np.log2((time_means[time2] + pseudocount) / (time_means[time1] + pseudocount)))
                    
                    if fold_change > max_fold_change:
                        max_fold_change = fold_change
                        max_time_pair = (time1, time2)
            
            # Add results
            result = {
                'Taxon': taxon,
                'P-value': p_value,
                'Test': 'Kruskal-Wallis',
                'Max Log2 Fold Change': max_fold_change,
            }
            
            # Add mean for each time point
            for time_point in time_points:
                result[f'Mean at {time_point}'] = time_means.get(time_point, np.nan)
            
            # Add fold change information
            if max_time_pair:
                result['Max Fold Change Times'] = f'{max_time_pair[0]} vs {max_time_pair[1]}'
                result['Fold Change'] = 2**max_fold_change
            
            results.append(result)
            
        except Exception as e:
            print(f"    Error testing {taxon}: {str(e)}")
            continue
    
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
        if not results_df.empty:
            results_df['Adjusted P-value'] = results_df['P-value']
    
    # Sort by adjusted p-value
    if not results_df.empty:
        results_df = results_df.sort_values('Adjusted P-value')
        
        # Print top significant results
        sig_results = results_df[results_df['Adjusted P-value'] < adjusted_p_threshold]
        
        if not sig_results.empty:
            print(f"\nFound {len(sig_results)} taxa with significant changes across time points:")
            for i, (_, row) in enumerate(sig_results.head(10).iterrows()):
                if i < 10:  # Only show top 10
                    print(f"  {i+1}. {row['Taxon']}: adj. p-value = {row['Adjusted P-value']:.4f}")
    
    return results_df

def main():
    """Main function for RSV microbiome analysis."""
    # Parse arguments
    args = parse_args()
    
    # Set up output directories
    output_dir = Path(args.output_dir)
    abundance_dir = output_dir / "abundance_tables"
    results_dir = output_dir / "results"
    figures_dir = output_dir / "figures"
    boxplots_dir = figures_dir / "boxplots"
    
    # Create output directories
    for dir_path in [output_dir, abundance_dir, results_dir, figures_dir, boxplots_dir]:
        os.makedirs(dir_path, exist_ok=True)
    
    print(f"Output will be saved to {output_dir}")
    
    # Process Sylph files for bacteria, viruses, and fungi
    print("\n1. Processing Sylph output files")
    
    # Find bacterial profiles
    bacteria_pattern = os.path.join(args.input_dir, "*profiled_bacteria.tsv")
    bacteria_files = glob.glob(bacteria_pattern)
    
    if bacteria_files:
        print(f"\nProcessing {len(bacteria_files)} bacterial profile files...")
        bacteria_df = combine_sylph_samples(
            bacteria_files, 
            abundance_type=args.abundance_type,
            min_ani=args.min_ani,
            min_coverage=args.min_coverage
        )
        
        if not bacteria_df.empty:
            # Save combined table
            bacteria_file = abundance_dir / 'bacteria_abundance.csv'
            bacteria_df.to_csv(bacteria_file)
            print(f"Bacterial abundance table saved to {bacteria_file}")
            print(f"Table contains {len(bacteria_df.index)} taxa across {len(bacteria_df.columns)} samples")
    else:
        print("No bacterial profile files found")
        bacteria_df = pd.DataFrame()
    
    # Find viral profiles
    virus_pattern = os.path.join(args.input_dir, "*profiled_viruses.tsv")
    virus_files = glob.glob(virus_pattern)
    
    if virus_files:
        print(f"\nProcessing {len(virus_files)} viral profile files...")
        virus_df = combine_sylph_samples(
            virus_files, 
            abundance_type=args.abundance_type,
            min_ani=args.min_ani,
            min_coverage=args.min_coverage
        )
        
        if not virus_df.empty:
            # Save combined table
            virus_file = abundance_dir / 'virus_abundance.csv'
            virus_df.to_csv(virus_file)
            print(f"Viral abundance table saved to {virus_file}")
            print(f"Table contains {len(virus_df.index)} taxa across {len(virus_df.columns)} samples")
    else:
        print("No viral profile files found")
        virus_df = pd.DataFrame()
    
    # Find fungal profiles
    fungi_pattern = os.path.join(args.input_dir, "*profiled_fungi.tsv")
    fungi_files = glob.glob(fungi_pattern)
    
    if fungi_files:
        print(f"\nProcessing {len(fungi_files)} fungal profile files...")
        fungi_df = combine_sylph_samples(
            fungi_files, 
            abundance_type=args.abundance_type,
            min_ani=args.min_ani,
            min_coverage=args.min_coverage
        )
        
        if not fungi_df.empty:
            # Save combined table
            fungi_file = abundance_dir / 'fungi_abundance.csv'
            fungi_df.to_csv(fungi_file)
            print(f"Fungal abundance table saved to {fungi_file}")
            print(f"Table contains {len(fungi_df.index)} taxa across {len(fungi_df.columns)} samples")
    else:
        print("No fungal profile files found")
        fungi_df = pd.DataFrame()
    
    # Create a combined table with all microbes
    print("\nCreating combined microbial abundance table...")
    combined_tables = []
    
    if not bacteria_df.empty:
        bacteria_df.index = ['Bacteria: ' + idx for idx in bacteria_df.index]
        combined_tables.append(bacteria_df)
        
    if not virus_df.empty:
        virus_df.index = ['Virus: ' + idx for idx in virus_df.index]
        combined_tables.append(virus_df)
        
    if not fungi_df.empty:
        fungi_df.index = ['Fungus: ' + idx for idx in fungi_df.index]
        combined_tables.append(fungi_df)
    
    if combined_tables:
        all_microbes_df = pd.concat(combined_tables)
        all_microbes_file = abundance_dir / 'all_microbes_abundance.csv'
        all_microbes_df.to_csv(all_microbes_file)
        print(f"Combined microbial abundance table saved to {all_microbes_file}")
        print(f"Table contains {len(all_microbes_df.index)} taxa across {len(all_microbes_df.columns)} samples")
    else:
        print("No data to combine into a microbial abundance table.")
        return
    
    # Load metadata
    print("\n2. Loading metadata")
    if os.path.exists(args.metadata):
        metadata_df = load_metadata(args.metadata, args.sample_id_column)
        
        if metadata_df.empty:
            print("Error: Could not load metadata")
            return
            
        print(f"Loaded metadata for {len(metadata_df)} samples")
        
        # Check sample overlap
        common_samples = set(all_microbes_df.columns).intersection(set(metadata_df.index))
        print(f"Found {len(common_samples)} samples in both abundance data and metadata")
        
        if len(common_samples) < 5:
            print("Error: Not enough common samples for analysis")
            return
            
    else:
        print(f"Error: Metadata file not found at {args.metadata}")
        return
    
    # Filter to keep only common samples
    all_microbes_df = all_microbes_df[list(common_samples)]
    
    # Filter low-abundance and low-prevalence taxa
    print("\n3. Filtering low-abundance and low-prevalence taxa")
    filtered_df = filter_low_abundance(
        all_microbes_df, 
        min_prevalence=args.min_prevalence,
        min_abundance=args.min_abundance
    )
    
    # Save filtered abundance table
    filtered_file = abundance_dir / 'filtered_abundance.csv'
    filtered_df.to_csv(filtered_file)
    print(f"Filtered abundance table saved to {filtered_file}")
    
    # Analyze by timing (Prior, Acute, Post)
    print("\n4. Analyzing changes across timing (Prior, Acute, Post)")
    timing_var = 'Timing'  # Adjust if your metadata uses a different column name
    
    if timing_var not in metadata_df.columns:
        print(f"Error: Timing variable '{timing_var}' not found in metadata")
        return
        
    timing_results = analyze_by_timing(filtered_df, metadata_df, timing_var)
    
    if not timing_results.empty:
        # Save timing results
        timing_file = results_dir / 'timing_differential_abundance.csv'
        timing_results.to_csv(timing_file)
        print(f"Timing analysis results saved to {timing_file}")
        
        # Create plots for top taxa
        significant_taxa = timing_results[timing_results['Adjusted P-value'] < 0.05]['Taxon'].tolist()
        
        if significant_taxa:
            print(f"\nCreating boxplots for {len(significant_taxa)} taxa with significant timing changes")
            
            for i, taxon in enumerate(significant_taxa[:10]):  # Plot top 10
                try:
                    output_file = boxplots_dir / f"timing_{taxon.replace(' ', '_').replace(':', '_')}.pdf"
                    plot_timing_boxplot(filtered_df, metadata_df, taxon, timing_var, output_file)
                except Exception as e:
                    print(f"Error creating plot for {taxon}: {str(e)}")
    
    # Analyze by severity and symptoms at each time point
    print("\n5. Analyzing differences by severity and symptoms at each time point")
    group_vars = ['Severity', 'Symptoms']  # Adjust if your metadata uses different column names
    
    # Check if group variables exist in metadata
    missing_vars = [var for var in group_vars if var not in metadata_df.columns]
    if missing_vars:
        print(f"Warning: The following group variables are missing from metadata: {missing_vars}")
        group_vars = [var for var in group_vars if var not in missing_vars]
        
    if not group_vars:
        print("Error: No valid group variables for analysis")
        return
    
    # Perform analysis by time point for each group variable
    time_point_results = analyze_by_timepoint(filtered_df, metadata_df, timing_var, group_vars)
    
    # Save results for each time point and group variable
    for time_point, group_results in time_point_results.items():
        for group_var, results_df in group_results.items():
            if not results_df.empty:
                # Save results
                result_file = results_dir / f"{time_point}_{group_var}_differential_abundance.csv"
                results_df.to_csv(result_file)
                print(f"Results for {time_point} by {group_var} saved to {result_file}")
                
                # Create plots for significant taxa
                significant_taxa = results_df[results_df['Adjusted P-value'] < 0.05]['Taxon'].tolist()
                
                if significant_taxa:
                    print(f"\nCreating boxplots for {len(significant_taxa)} taxa with significant {group_var} differences at {time_point}")
                    
                    for i, taxon in enumerate(significant_taxa[:10]):  # Plot top 10
                        try:
                            output_file = boxplots_dir / f"{time_point}_{group_var}_{taxon.replace(' ', '_').replace(':', '_')}.pdf"
                            plot_taxa_boxplot(filtered_df, metadata_df, taxon, timing_var, group_var, output_file)
                        except Exception as e:
                            print(f"Error creating plot for {taxon}: {str(e)}")
    
    print("\nAnalysis complete!")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"Error in main execution: {str(e)}")
        traceback.print_exc()

