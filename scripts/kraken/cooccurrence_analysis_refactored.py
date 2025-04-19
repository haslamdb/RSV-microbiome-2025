#!/usr/bin/env python3
"""
Analyze co-occurrence patterns in RSV microbiome data using Kraken/Bracken abundance data.

This script:
1. Loads the previously generated abundance data from Kraken/Bracken
2. Extracts abundance data for specified bacterial species (default: S. pneumoniae and H. influenzae)
3. Analyzes their abundance patterns:
   - Across time points (Prior, Acute, Post)
   - Between severity groups at each time point
   - Between symptom groups at each time point
4. Examines co-occurrence patterns of both species
5. Generates visualizations regardless of statistical significance

Usage:
    python cooccurrence_analysis.py [--input-file FILTERED_ABUNDANCE_FILE] [--metadata METADATA_FILE]
"""

import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
import os
import sys
import re
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import mannwhitneyu, kruskal, spearmanr, ttest_1samp
from statsmodels.stats.multitest import multipletests
import traceback

# Add project root to Python path
project_root = Path(__file__).resolve().parents[2]
sys.path.append(str(project_root))

# Add kraken_tools to Python path
kraken_tools_dir = Path.home() / "Documents" / "Code" / "kraken_tools"
sys.path.append(str(kraken_tools_dir))

# Import functions from kraken_tools if available
try:
    from kraken_tools.utils.sample_utils import check_input_files_exist
    from kraken_tools.analysis.abundance import process_bracken_abundance, normalize_abundance
    from kraken_tools.analysis.metadata import read_and_process_metadata
    from kraken_tools.analysis.taxonomy import read_and_process_taxonomy
    from kraken_tools.analysis.statistical import run_statistical_tests
    kraken_tools_available = True
except ImportError:
    kraken_tools_available = False
    logging.info("Warning: kraken_tools module not available. Using built-in functions.")

# Set plotting style
plt.style.use('seaborn-v0_8-whitegrid')
sns.set(font_scale=1.1)

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Analyze co-occurrence patterns in RSV microbiome data using Kraken/Bracken')
    parser.add_argument('--input-file', type=str, default='results/kraken_analysis/filtered_kraken_s_abundance.tsv',
                      help='Path to filtered abundance table')
    parser.add_argument('--metadata', type=str, default='metadata.csv',
                      help='Path to metadata file')
    parser.add_argument('--output-dir', type=str, default='results/kraken_cooccurrence_analysis',
                      help='Directory to save analysis results')
    parser.add_argument('--sample-id-column', type=str, default='SampleID',
                      help='Column name for sample IDs in metadata')
    parser.add_argument('--target-species', type=str, nargs='+', 
                       default=['Streptococcus pneumoniae', 'Haemophilus influenzae', 'Moraxella catarrhalis'],
                       help='Target species to analyze for co-occurrence')
    parser.add_argument('--only-cooccurrence', action='store_true', 
                       help='Run only co-occurrence analysis')
    parser.add_argument('--skip-plots', action='store_true', 
                       help='Skip all plotting steps')
    parser.add_argument('--normalization', type=str, default='none',
                      choices=['none', 'clr', 'rarefaction', 'tss'],
                      help='Normalization method to use')
    parser.add_argument('--rarefaction-depth', type=int, default=None,
                      help='Sequencing depth for rarefaction (default: use minimum sample depth)')
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
            logging.info(f"Warning: Found {metadata_df.index.duplicated().sum()} duplicate sample IDs in metadata")
            metadata_df = metadata_df[~metadata_df.index.duplicated(keep='first')]
        
        # Convert categorical variables to string
        for col in metadata_df.columns:
            if metadata_df[col].dtype == 'object' or metadata_df[col].dtype.name == 'category':
                metadata_df[col] = metadata_df[col].astype(str)
        
        return metadata_df
    
    except Exception as e:
        logging.error(f"Error loading metadata file: {str(e)}")
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
        # Use word boundaries in regex for more precise matching
        pattern = r'\b' + re.escape(species.lower()) + r'\b'
        matches = [idx for idx in abundance_df.index if re.search(pattern, idx.lower())]

        if matches:
            # Use the first match as the index
            species_indices[species] = matches[0]
            logging.info(f"Found '{species}' as '{matches[0]}'")
            if len(matches) > 1:
                logging.info(f"  Note: Multiple matches found: {matches}")
        else:
            logging.warning(f"Could not find '{species}' in the data")
    
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

def rarefy_abundance_data(abundance_df, rarefaction_depth=None, seed=42):
    """
    Rarefy abundance data to a specified sequencing depth.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Taxa abundance DataFrame with taxa as index, samples as columns
        Values should be raw counts (not percentages or normalized values)
    rarefaction_depth : int, optional
        Sequencing depth to rarefy to. If None, uses the minimum sample depth.
    seed : int, optional
        Random seed for reproducibility
        
    Returns:
    --------
    pandas.DataFrame
        Rarefied abundance DataFrame
    """
    # Set random seed for reproducibility
    np.random.seed(seed)
    
    # Make a copy of the input data
    rarefied_df = abundance_df.copy()
    
    # Calculate the total counts for each sample
    sample_totals = rarefied_df.sum(axis=0)
    
    # If no rarefaction depth is given, use the minimum sample total
    if rarefaction_depth is None:
        rarefaction_depth = int(sample_totals.min())
        logging.info(f"Using minimum sample depth: {rarefaction_depth}")
    else:
        # Always exclude samples with fewer reads than the rarefaction depth
        samples_below_depth = sample_totals[sample_totals < rarefaction_depth].index.tolist()
        if samples_below_depth:
            logging.info(f"Excluding {len(samples_below_depth)} samples with fewer reads than the rarefaction depth ({rarefaction_depth})")
            logging.info(f"Excluded samples: {samples_below_depth}")
            rarefied_df = rarefied_df.drop(columns=samples_below_depth)
            # Recalculate sample totals after dropping samples
            sample_totals = rarefied_df.sum(axis=0)
            logging.info(f"Proceeding with {len(rarefied_df.columns)} samples")
    
    logging.info(f"Rarefying to depth of {rarefaction_depth} reads per sample")
    
    # Perform rarefaction for each sample
    for sample in rarefied_df.columns:
        # Get current counts for this sample
        counts = rarefied_df[sample].values
        
        # Create an array representing all individual reads
        reads_pool = np.repeat(np.arange(len(counts)), counts.astype(int))
        
        # Randomly select reads up to the rarefaction depth
        if len(reads_pool) > rarefaction_depth:
            # Randomly select indices without replacement
            selected_indices = np.random.choice(reads_pool, rarefaction_depth, replace=False)
            
            # Count occurrences of each taxon
            rarefied_counts = np.bincount(selected_indices, minlength=len(counts))
            
            # Update the dataframe with rarefied counts
            rarefied_df[sample] = rarefied_counts
    
    return rarefied_df

def apply_tss_normalization(abundance_df):
    """
    Apply Total Sum Scaling (TSS) normalization to abundance data (counts per million).
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Taxa abundance DataFrame with taxa as index, samples as columns
        
    Returns:
    --------
    pandas.DataFrame
        TSS normalized abundance DataFrame
    """
    # Make a copy of the input data
    normalized_df = abundance_df.copy()
    
    # Calculate the sum of counts for each sample
    sample_sums = normalized_df.sum(axis=0)
    
    # Divide each count by the sample sum and multiply by 1e6 (counts per million)
    normalized_df = normalized_df.div(sample_sums, axis=1) * 1e6
    
    return normalized_df

def apply_clr_transformation(abundance_df):
    """
    Apply Centered Log-Ratio (CLR) transformation to abundance data.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Taxa abundance DataFrame with taxa as index, samples as columns
        
    Returns:
    --------
    pandas.DataFrame
        CLR transformed abundance DataFrame
    """
    # Make a copy of the input data
    transformed_df = abundance_df.copy()
    
    # Add small pseudocount to zeros
    df_pseudo = transformed_df.replace(0, np.nextafter(0, 1))
    
    # Log transform
    df_log = np.log(df_pseudo)
    
    # Subtract column-wise mean (CLR transformation)
    transformed_df = df_log.subtract(df_log.mean(axis=0), axis=1)
    
    return transformed_df

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

    # Apply multiple testing correction
    if not results_df.empty and 'p-value' in results_df.columns and not results_df['p-value'].isna().all():
        try:
            # Filter out NaN values for correction
            valid_pvals = ~results_df['p-value'].isna()
            if valid_pvals.any():
                p_values = results_df.loc[valid_pvals, 'p-value'].values
                corrected = multipletests(p_values, method='fdr_bh')
                
                # Create a column for adjusted p-values (initialize with NaN)
                results_df['p-adjusted'] = np.nan
                
                # Update only the valid positions
                results_df.loc[valid_pvals, 'p-adjusted'] = corrected[1]
        except Exception as e:
            logging.warning(f"Could not apply multiple testing correction: {str(e)}")

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
        logging.info(f"Analyzing {time_point} by {group_var}")
        
        # Get samples for this time point
        time_samples = metadata_df[metadata_df[time_var] == time_point].index
        
        # Get common samples with abundance data
        common_samples = list(set(time_samples).intersection(set(abundance_df.columns)))
        
        if len(common_samples) < 3:
            logging.info(f"  Not enough samples for time point {time_point} (found {len(common_samples)}), skipping")
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

        # Apply multiple testing correction
        if not results_df.empty and 'p-value' in results_df.columns and not results_df['p-value'].isna().all():
            try:
                # Filter out NaN values for correction
                valid_pvals = ~results_df['p-value'].isna()
                if valid_pvals.any():
                    p_values = results_df.loc[valid_pvals, 'p-value'].values
                    corrected = multipletests(p_values, method='fdr_bh')
                    
                    # Create a column for adjusted p-values (initialize with NaN)
                    results_df['p-adjusted'] = np.nan
                    
                    # Update only the valid positions
                    results_df.loc[valid_pvals, 'p-adjusted'] = corrected[1]
            except Exception as e:
                logging.warning(f"Could not apply multiple testing correction: {str(e)}")
        
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
        logging.info("Need both species for co-occurrence analysis")
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
            logging.info(f"Error calculating p-value: {str(e)}")
        
        # Set title and labels
        ax.set_title(f'{species_name} Abundance Across Time Points', fontsize=14, fontweight='bold')
        ax.set_xlabel(time_var, fontsize=12)
        ax.set_ylabel('Abundance', fontsize=12)
        
        # Adjust layout
        plt.tight_layout()
        
        # Save figure
        output_file = os.path.join(output_dir, f"{species_name.replace(' ', '_')}_by_time.pdf")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logging.info(f"Saved plot to {output_file}")
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
    # Get common samples
    common_samples = list(set(abundance_df.columns).intersection(set(metadata_df.index)))
    
    # Get taxon abundance
    if taxon not in abundance_df.index:
        logging.info(f"Error: Taxon '{taxon}' not found in abundance data")
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
    
    # Create categorical variables with the correct order
    plot_data[time_var] = pd.Categorical(
        plot_data[time_var],
        categories=time_order,
        ordered=True
    )
    
    # Set order for Symptoms if that's the group_var
    if group_var == 'Symptoms':
        symptoms_order = ['Asymptomatic', 'Mild', 'Severe']
        # Convert to categorical with specific order
        plot_data[group_var] = pd.Categorical(
            plot_data[group_var],
            categories=symptoms_order,
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
            # Remove title and add time label on right
            axes[row_idx].set_title("")
            # Create right y-axis for time label
            ax_right = axes[row_idx].twinx()
            ax_right.set_yticks([])
            ax_right.set_ylabel(time_point, rotation=-90, labelpad=15, fontsize=12, fontweight='bold')
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
                    logging.info(f"Error calculating p-value: {str(e)}")
                    pass
        
        # Remove title from the main axis
        axes[row_idx].set_title("")
        
        # Create a twin axis on the right side
        ax_right = axes[row_idx].twinx()
        
        # Hide the tick marks on the right y-axis
        ax_right.set_yticks([])
        
        # Add time point as right-side ylabel rotated -90 degrees
        ax_right.set_ylabel(time_point, rotation=-90, labelpad=15, fontsize=12, fontweight='bold')
        
        # Set y-label for main axis (only for middle subplot)
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
        logging.info(f"Plot saved to {output_file}")
    
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
    # Check if we have at least two species
    if len(species_indices) < 2:
        logging.info("Need at least two species for co-occurrence heatmap")
        return
    
    # Get common samples
    common_samples = list(set(abundance_df.columns).intersection(set(metadata_df.index)))
    
    # Extract abundance for species
    species_names = list(species_indices.keys())
    
    # Create a dataframe to hold the data for the heatmap
    plot_data = pd.DataFrame()
    
    # Add data for each species that exists in the abundance dataframe
    for species_name in species_names:
        taxon = species_indices[species_name]
        if taxon in abundance_df.index:
            plot_data[species_name] = abundance_df.loc[taxon, common_samples]
    
    # If we don't have at least 2 species with data, we can't make a heatmap
    if plot_data.shape[1] < 2:
        logging.info("Not enough species with data to create a heatmap")
        return
    
    # Convert to numeric for heatmap (log transform if values are positive)
    plot_data_log = plot_data.apply(lambda x: np.log1p(x) if (x > 0).any() else x)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Create heatmap
    ax = sns.heatmap(
        plot_data_log.T,  # Transpose so species are rows
        cmap='YlGnBu',
        cbar_kws={'label': 'Log-transformed Abundance'},
        ax=ax
    )
    
    # The y-ticks are now automatically set to the correct number by seaborn
    # Just modify the labels to use friendly names
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    
    # Get time information if available
    if 'Timing' in metadata_df.columns:
        # Add time point labels above the heatmap
        time_colors = {'Prior': '#1f77b4', 'Acute': '#ff7f0e', 'Post': '#2ca02c'}
        
        for i, sample in enumerate(common_samples):
            if sample in metadata_df.index:
                time_point = metadata_df.loc[sample, 'Timing']
                if time_point in time_colors:
                    ax.add_patch(plt.Rectangle(
                        (i, plot_data.shape[1]),
                        1, 0.2,
                        color=time_colors[time_point],
                        ec='none'
                    ))
        
        # Add legend for time points
        from matplotlib.patches import Patch
        legend_elements = [Patch(facecolor=color, label=time) 
                          for time, color in time_colors.items()
                          if time in metadata_df['Timing'].values]
        
        if legend_elements:
            ax.legend(handles=legend_elements, loc='upper center', 
                     bbox_to_anchor=(0.5, 1.15), ncol=len(legend_elements))
    
    # Set title
    plt.title('Co-occurrence of Target Species', fontsize=14)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save figure
    output_file = os.path.join(output_dir, 'co_occurrence_heatmap.pdf')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    logging.info(f"Saved co-occurrence heatmap to {output_file}")
    plt.close(fig)

def analyze_prior_acute_delta(abundance_df, metadata_df, species_indices, output_dir):
    """
    Calculate and analyze the delta (change) in abundance from Prior to Acute timepoints.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Taxa abundance DataFrame with taxa as index, samples as columns
    metadata_df : pandas.DataFrame
        Metadata DataFrame with samples as index
    species_indices : dict
        Dictionary mapping species names to their indices in the abundance DataFrame
    output_dir : str
        Directory to save analysis results
        
    Returns:
    --------
    pandas.DataFrame
        DataFrame containing delta values for each subject and species
    """
    logging.info("\nAnalyzing Prior to Acute abundance changes")
    
    # Get subject ID column
    subject_col = 'SubjectID' if 'SubjectID' in metadata_df.columns else None
    if not subject_col:
        logging.error("SubjectID column not found in metadata")
        return None
    
    # Get timing column
    time_var = 'Timing'
    if time_var not in metadata_df.columns:
        logging.error(f"{time_var} column not found in metadata")
        return None
    
    # Filter to only include samples from Prior and Acute timepoints
    prior_acute_samples = metadata_df[metadata_df[time_var].isin(['Prior', 'Acute'])]
    
    # Create a dictionary to hold delta values for each subject and species
    deltas = {}
    
    # For each subject with both Prior and Acute samples
    subject_ids = prior_acute_samples[subject_col].unique()
    for subject in subject_ids:
        subject_samples = prior_acute_samples[prior_acute_samples[subject_col] == subject]
        
        # Check if subject has both Prior and Acute samples
        if len(subject_samples[subject_samples[time_var] == 'Prior']) >= 1 and len(subject_samples[subject_samples[time_var] == 'Acute']) >= 1:
            # Get sample IDs
            prior_sample = subject_samples[subject_samples[time_var] == 'Prior'].index[0]
            acute_sample = subject_samples[subject_samples[time_var] == 'Acute'].index[0]
            
            # Only process if both samples are in the abundance data
            if prior_sample in abundance_df.columns and acute_sample in abundance_df.columns:
                # For each target species, calculate delta
                for species_name, taxon in species_indices.items():
                    prior_value = abundance_df.loc[taxon, prior_sample]
                    acute_value = abundance_df.loc[taxon, acute_sample]
                    delta = acute_value - prior_value
                    
                    # Store in deltas dictionary
                    if species_name not in deltas:
                        deltas[species_name] = {}
                    deltas[species_name][subject] = {
                        'prior': prior_value,
                        'acute': acute_value,
                        'delta': delta,
                        'subject_id': subject
                    }
                    
                    # Get additional metadata if available
                    for col in ['Severity', 'Symptoms']:
                        if col in metadata_df.columns:
                            # Use the value from Acute sample
                            if acute_sample in metadata_df.index and pd.notna(metadata_df.loc[acute_sample, col]):
                                deltas[species_name][subject][col] = metadata_df.loc[acute_sample, col]
                            # Fallback to Prior sample if Acute doesn't have the value
                            elif prior_sample in metadata_df.index and pd.notna(metadata_df.loc[prior_sample, col]):
                                deltas[species_name][subject][col] = metadata_df.loc[prior_sample, col]
    
    # Convert nested dictionary to DataFrames for easier analysis
    delta_dfs = {}
    tables_dir = Path(output_dir) / "tables"
    figures_dir = Path(output_dir) / "figures"
    
    for species_name in deltas:
        # Convert to DataFrame
        df_data = [data for subject, data in deltas[species_name].items()]
        if df_data:
            delta_df = pd.DataFrame(df_data)
            delta_dfs[species_name] = delta_df
            
            # Save to file
            delta_file = tables_dir / f"{species_name.replace(' ', '_')}_prior_acute_delta.csv"
            delta_df.to_csv(delta_file, index=False)
            logging.info(f"Delta analysis for {species_name} saved to {delta_file}")
            
            # Create visualizations
            
            # 1. Simple boxplot of delta values
            plt.figure(figsize=(8, 6))
            plt.boxplot(delta_df['delta'].values, labels=[species_name])
            plt.axhline(y=0, color='r', linestyle='-', alpha=0.3)
            
            # Add individual points
            plt.scatter(np.ones(len(delta_df)) + np.random.normal(0, 0.05, size=len(delta_df)), 
                        delta_df['delta'].values, color='black', alpha=0.5)
            
            # Add p-value from one-sample t-test against 0
            try:
                t_stat, p_value = ttest_1samp(delta_df['delta'].values, 0)
                plt.title(f'Change in {species_name} from Prior to Acute\np={p_value:.3f}')
            except:
                plt.title(f'Change in {species_name} from Prior to Acute')
            
            plt.ylabel('Delta (Acute - Prior)')
            plt.tight_layout()
            plt.savefig(figures_dir / f"{species_name.replace(' ', '_')}_delta_boxplot.pdf", bbox_inches='tight')
            plt.close()
            
            # 2. Faceted boxplots by severity and symptoms if available
            for facet_var in ['Severity', 'Symptoms']:
                if facet_var in delta_df.columns:
                    plt.figure(figsize=(10, 6))
                    
                    # Define category order if possible
                    if facet_var == 'Symptoms':
                        category_order = ['Asymptomatic', 'Mild', 'Severe']
                        # Filter to valid categories
                        valid_cats = [cat for cat in category_order if cat in delta_df[facet_var].unique()]
                    else:
                        valid_cats = sorted(delta_df[facet_var].unique())
                    
                    # Create plot data
                    plot_data = []
                    labels = []
                    positions = []
                    
                    for i, cat in enumerate(valid_cats):
                        cat_data = delta_df[delta_df[facet_var] == cat]['delta'].values
                        if len(cat_data) > 0:
                            plot_data.append(cat_data)
                            labels.append(f"{cat} (n={len(cat_data)})")
                            positions.append(i+1)
                    
                    # Create boxplot
                    if plot_data:
                        plt.boxplot(plot_data, labels=labels, positions=positions)
                        plt.axhline(y=0, color='r', linestyle='-', alpha=0.3)
                        
                        # Add individual points
                        for i, data in enumerate(plot_data):
                            plt.scatter(np.ones(len(data)) * (i+1) + np.random.normal(0, 0.05, size=len(data)), 
                                        data, color='black', alpha=0.5)
                        
                        # Add p-values from Mann-Whitney U tests if we have at least 2 groups
                        if len(plot_data) >= 2:
                            try:
                                from scipy.stats import mannwhitneyu
                                for i in range(len(plot_data) - 1):
                                    for j in range(i + 1, len(plot_data)):
                                        stat, p_value = mannwhitneyu(plot_data[i], plot_data[j])
                                        plt.text((positions[i] + positions[j]) / 2, 
                                                max(np.max(plot_data[i]), np.max(plot_data[j])) * 1.1,
                                                f'p={p_value:.3f}',
                                                ha='center')
                            except Exception as e:
                                logging.info(f"Error calculating p-values: {str(e)}")
                        
                        plt.title(f'Change in {species_name} by {facet_var}')
                        plt.ylabel('Delta (Acute - Prior)')
                        plt.tight_layout()
                        plt.savefig(figures_dir / f"{species_name.replace(' ', '_')}_delta_by_{facet_var}.pdf", bbox_inches='tight')
                    plt.close()
    
    return delta_dfs

def plot_subject_trajectory(abundance_df, metadata_df, species_indices, output_dir):
    """
    Create trajectory plots showing individual subject changes from Prior to Acute to Post.
    
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
    logging.info("\nCreating subject trajectory plots")
    
    # Get subject ID column
    subject_col = 'SubjectID' if 'SubjectID' in metadata_df.columns else None
    if not subject_col:
        logging.error("SubjectID column not found in metadata")
        return
    
    # Get timing column
    time_var = 'Timing'
    if time_var not in metadata_df.columns:
        logging.error(f"{time_var} column not found in metadata")
        return
    
    # Define time order
    time_order = ['Prior', 'Acute', 'Post']
    
    # Output directory for figures
    figures_dir = Path(output_dir) / "figures"
    
    # Process each species
    for species_name, taxon in species_indices.items():
        logging.info(f"Creating trajectory plot for {species_name}")
        
        # Get all subjects that have at least 2 timepoints
        subject_timepoints = {}
        
        # First, collect all available data points
        for subject in metadata_df[subject_col].unique():
            subject_samples = metadata_df[metadata_df[subject_col] == subject]
            
            # Get samples for each timepoint that exist in abundance data
            timepoint_data = {}
            for time_point in time_order:
                if time_point in subject_samples[time_var].values:
                    time_samples = subject_samples[subject_samples[time_var] == time_point].index
                    valid_samples = [s for s in time_samples if s in abundance_df.columns]
                    
                    if valid_samples:
                        # Use the first valid sample
                        sample_id = valid_samples[0]
                        abundance_value = abundance_df.loc[taxon, sample_id]
                        timepoint_data[time_point] = {
                            'sample_id': sample_id,
                            'abundance': abundance_value
                        }
                        
                        # Add metadata if available
                        for col in ['Severity', 'Symptoms']:
                            if col in metadata_df.columns:
                                if sample_id in metadata_df.index and pd.notna(metadata_df.loc[sample_id, col]):
                                    timepoint_data[time_point][col] = metadata_df.loc[sample_id, col]
            
            # Only include subjects with at least 2 timepoints
            if len(timepoint_data) >= 2:
                subject_timepoints[subject] = timepoint_data
        
        # Create the trajectory plot
        plt.figure(figsize=(12, 8))
        
        # Map timepoints to x-axis positions
        time_to_x = {time: i for i, time in enumerate(time_order)}
        
        # Process each subject
        for subject, timepoints in subject_timepoints.items():
            # Extract x and y coordinates for this subject
            x_values = [time_to_x[time] for time in timepoints.keys()]
            y_values = [data['abundance'] for time, data in timepoints.items()]
            
            # Get color based on metadata if available (using first available timepoint)
            color = 'gray'  # default color
            
            if 'Severity' in next(iter(timepoints.values())):
                severity = next(iter(timepoints.values()))['Severity']
                if severity == 'Severe':
                    color = 'red'
                elif severity == 'Moderate':
                    color = 'orange'
                elif severity == 'Mild':
                    color = 'blue'
            
            # Plot subject trajectory
            plt.plot(x_values, y_values, 'o-', color=color, alpha=0.5, linewidth=1)
        
        # Add boxplots for each timepoint
        boxplot_data = []
        for time_point in time_order:
            # Collect all abundance values for this timepoint
            time_values = []
            for subject, timepoints in subject_timepoints.items():
                if time_point in timepoints:
                    time_values.append(timepoints[time_point]['abundance'])
            boxplot_data.append(time_values)
        
        # Add boxplots at each timepoint
        boxplot_positions = list(range(len(time_order)))
        boxplots = plt.boxplot(boxplot_data, positions=boxplot_positions, 
                              widths=0.5, patch_artist=True)
        
        # Set box colors to light gray for better visibility of the trajectories
        for box in boxplots['boxes']:
            box.set(facecolor='lightgray', alpha=0.3)
        
        # Set axis labels and title
        plt.xlabel('Timepoint')
        plt.ylabel('Abundance')
        plt.title(f'{species_name} Abundance Trajectory by Subject')
        plt.xticks(boxplot_positions, time_order)
        
        # Add legend if we used colors
        if 'Severity' in next(iter(next(iter(subject_timepoints.values())).values()), {}):
            from matplotlib.patches import Patch
            legend_elements = [
                Patch(facecolor='red', label='Severe'),
                Patch(facecolor='orange', label='Moderate'),
                Patch(facecolor='blue', label='Mild'),
                Patch(facecolor='gray', label='Unknown')
            ]
            plt.legend(handles=legend_elements, loc='upper right')
        
        # Save figure
        plt.tight_layout()
        output_file = figures_dir / f"{species_name.replace(' ', '_')}_subject_trajectory.pdf"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logging.info(f"Saved trajectory plot to {output_file}")
        plt.close()
        
        # Create faceted trajectory plots if we have metadata
        for facet_var in ['Severity', 'Symptoms']:
            # Check if at least one subject has this metadata
            has_metadata = False
            for subject, timepoints in subject_timepoints.items():
                for time_point, data in timepoints.items():
                    if facet_var in data:
                        has_metadata = True
                        break
                if has_metadata:
                    break
            
            if has_metadata:
                # Get unique values for this variable
                unique_values = set()
                for subject, timepoints in subject_timepoints.items():
                    for time_point, data in timepoints.items():
                        if facet_var in data:
                            unique_values.add(data[facet_var])
                
                # Define order for Symptoms if needed
                if facet_var == 'Symptoms':
                    symptoms_order = ['Asymptomatic', 'Mild', 'Severe']
                    # Filter to valid categories
                    ordered_values = [val for val in symptoms_order if val in unique_values]
                else:
                    ordered_values = sorted(unique_values)
                
                # Create a subplot for each value
                fig, axes = plt.subplots(1, len(ordered_values), figsize=(15, 6), sharey=True)
                if len(ordered_values) == 1:
                    axes = [axes]  # Make sure axes is always a list
                
                # Process each group
                for i, value in enumerate(ordered_values):
                    ax = axes[i]
                    
                    # Collect subjects for this group
                    group_subjects = []
                    for subject, timepoints in subject_timepoints.items():
                        # Check if any timepoint has this value
                        for time_point, data in timepoints.items():
                            if facet_var in data and data[facet_var] == value:
                                group_subjects.append(subject)
                                break
                    
                    # Plot trajectories for this group
                    for subject in group_subjects:
                        timepoints = subject_timepoints[subject]
                        # Extract x and y coordinates for this subject
                        x_values = [time_to_x[time] for time in timepoints.keys()]
                        y_values = [data['abundance'] for time, data in timepoints.items()]
                        ax.plot(x_values, y_values, 'o-', color='blue', alpha=0.5, linewidth=1)
                    
                    # Add boxplots for each timepoint
                    boxplot_data = []
                    for time_point in time_order:
                        # Collect all abundance values for this timepoint and group
                        time_values = []
                        for subject in group_subjects:
                            if subject in subject_timepoints and time_point in subject_timepoints[subject]:
                                time_values.append(subject_timepoints[subject][time_point]['abundance'])
                        boxplot_data.append(time_values)
                    
                    # Add boxplots at each timepoint
                    boxplot_positions = list(range(len(time_order)))
                    if any(len(values) > 0 for values in boxplot_data):
                        boxplots = ax.boxplot(boxplot_data, positions=boxplot_positions, 
                                             widths=0.5, patch_artist=True)
                        # Set box colors
                        for box in boxplots['boxes']:
                            box.set(facecolor='lightgray', alpha=0.3)
                    
                    # Set title and x-axis labels
                    ax.set_title(f'{value} (n={len(group_subjects)})')
                    ax.set_xticks(boxplot_positions)
                    ax.set_xticklabels(time_order)
                    
                    # Add ylabel only to the first subplot
                    if i == 0:
                        ax.set_ylabel('Abundance')
                
                # Set overall title
                fig.suptitle(f'{species_name} Abundance Trajectory by {facet_var}', fontsize=14)
                plt.tight_layout()
                fig.subplots_adjust(top=0.85)
                
                # Save figure
                output_file = figures_dir / f"{species_name.replace(' ', '_')}_trajectory_by_{facet_var}.pdf"
                plt.savefig(output_file, dpi=300, bbox_inches='tight')
                logging.info(f"Saved faceted trajectory plot by {facet_var} to {output_file}")
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
        logging.info("Need both species for correlation plot")
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
    logging.info(f"Saved correlation plot to {output_file}")
    plt.close()


def run_analysis_pipeline(args):
    """Main function for analyzing species co-occurrence in microbiome data."""
    # Set up output directories
    output_dir = Path(args.output_dir)
    tables_dir = output_dir / "tables"
    figures_dir = output_dir / "figures"
    
    # Create output directories
    for dir_path in [output_dir, tables_dir, figures_dir]:
        os.makedirs(dir_path, exist_ok=True)
    
    logging.info(f"Output will be saved to {output_dir}")
    
    # Load abundance data
    logging.info("\n1. Loading abundance data")
    if os.path.exists(args.input_file):
        # First read data without setting the index to determine its format
        raw_data = pd.read_csv(args.input_file, sep='\t')
        
        # Check if the data contains metadata columns
        metadata_cols = ['SampleID', 'SubjectID', 'CollectionDate', 'Timing', 'Severity', 'Symptoms']
        has_metadata_cols = any(col in raw_data.columns for col in metadata_cols)
        
        if has_metadata_cols:
            logging.info("Detected metadata columns in abundance file. Data appears to be in samples-as-rows format.")
            logging.info("Transposing data to get taxa-as-rows format required for analysis...")
            
            # Identify which metadata columns are present
            present_metadata_cols = [col for col in metadata_cols if col in raw_data.columns]
            
            # Extract metadata
            metadata = raw_data[present_metadata_cols].copy()
            
            # Extract abundance data
            abundance_cols = [col for col in raw_data.columns if col not in present_metadata_cols]
            abundance = raw_data[abundance_cols]
            
            # Transpose abundance data (taxa as rows, samples as columns)
            transposed = abundance.T
            transposed.columns = raw_data['SampleID'] if 'SampleID' in raw_data.columns else raw_data.index
            transposed.index.name = 'Species'
            
            # Use the transposed data
            abundance_df = transposed
            logging.info(f"Successfully transposed data: now has {len(abundance_df.index)} taxa as rows and {len(abundance_df.columns)} samples as columns")
        else:
            # Data is likely already in the correct format or has a non-standard structure
            abundance_df = pd.read_csv(args.input_file, sep='\t', index_col=0)
            
            # Check if data appears to be in the expected format
            likely_taxa_patterns = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom', 'domain',
                          'streptococcus', 'haemophilus', 'staphylococcus', 'moraxella', 'corynebacterium',
                          'bacteria', 'virus', 'fungi', 'bacteroides', 'bifidobacterium', 'escherichia']
            
            # Check if index entries contain species names with dots (e.g., "Streptococcus.pneumoniae")
            first_indices = abundance_df.index[:5].astype(str)
            lower_indices = [idx.lower() for idx in first_indices]
            
            # Look for taxonomy patterns or Species-like entries
            has_taxonomy_pattern = any(any(taxon in idx for taxon in likely_taxa_patterns) for idx in lower_indices)
            has_species_dot_format = any('.' in idx for idx in first_indices)
            
            if has_taxonomy_pattern or has_species_dot_format:
                logging.info(f"Data appears to be in the correct format with {len(abundance_df.index)} taxa and {len(abundance_df.columns)} samples")
            else:
                logging.warning("Data format unclear. Attempting to proceed, but results may not be as expected.")
                logging.info(f"Ideal format: Taxa as rows, samples as columns. Current shape: {abundance_df.shape}")
    else:
        logging.error(f"Input file not found at {args.input_file}")
        return None
    
    # Apply normalization if requested
    if args.normalization == 'rarefaction':
        logging.info("\nApplying rarefaction normalization")
        original_shape = abundance_df.shape
        abundance_df = rarefy_abundance_data(abundance_df, args.rarefaction_depth)
        logging.info(f"Data shape before rarefaction: {original_shape}")
        logging.info(f"Data shape after rarefaction: {abundance_df.shape}")
    elif args.normalization == 'clr':
        logging.info("\nApplying CLR transformation")
        abundance_df = apply_clr_transformation(abundance_df)
        logging.info(f"Applied CLR transformation to data")
    elif args.normalization == 'tss':
        logging.info("\nApplying Total Sum Scaling (TSS) normalization")
        abundance_df = apply_tss_normalization(abundance_df)
        logging.info(f"TSS normalization applied to {abundance_df.shape[1]} samples")
    else:
        logging.info("\nUsing raw counts (no normalization)")
    
    # Load metadata
    logging.info("\n2. Loading metadata")
    if os.path.exists(args.metadata):
        metadata_df = load_metadata(args.metadata, args.sample_id_column)
        
        if metadata_df.empty:
            logging.error("Could not load metadata")
            return None
            
        logging.info(f"Loaded metadata for {len(metadata_df)} samples")
        
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
            logging.warning(f"The following variables are missing from metadata: {missing_vars}")
    else:
        logging.error(f"Metadata file not found at {args.metadata}")
        return None
    
    # Find the target species in the data
    logging.info("\n3. Finding target species in the data")
    target_species = args.target_species
    species_indices = find_species_in_data(abundance_df, target_species)
    
    if not species_indices:
        logging.error("Could not find any target species in the data")
        return None
    
    # If only co-occurrence analysis is requested, skip other analyses
    if args.only_cooccurrence:
        logging.info("\nSkipping individual species analyses and performing only co-occurrence analysis")
        if len(species_indices) >= 2:
            logging.info("Analyzing co-occurrence patterns")
            co_occurrence_results = analyze_co_occurrence(abundance_df, metadata_df, species_indices)
            
            if co_occurrence_results:
                # Save results
                result_file = tables_dir / 'co_occurrence_results.csv'
                pd.DataFrame(co_occurrence_results['co_occurrence_rates']).T.to_csv(result_file)
                logging.info(f"Co-occurrence results saved to {result_file}")
                
                # Create visualizations if not skipped
                if not args.skip_plots:
                    logging.info("Creating co-occurrence visualizations")
                    plot_co_occurrence_heatmap(abundance_df, metadata_df, species_indices, figures_dir)
                    plot_scatter_correlation(abundance_df, metadata_df, species_indices, figures_dir)
        else:
            logging.error("At least two species are required for co-occurrence analysis")
            return None
    else:
        # Perform full analysis pipeline
        
        # Analysis by timing
        logging.info("\n4.1 Analyzing species abundance across timing")
        timing_results = analyze_species_by_variable(abundance_df, metadata_df, species_indices, time_var)
        
        if not timing_results.empty:
            # Save results
            timing_file = tables_dir / 'species_by_timing.csv'
            timing_results.to_csv(timing_file)
            logging.info(f"Timing analysis results saved to {timing_file}")
            
            # Create plots if not skipped
            if not args.skip_plots:
                logging.info("Creating boxplots for species abundance by timing")
                plot_species_by_time(abundance_df, metadata_df, species_indices, time_var, figures_dir)
        
        # Analysis by severity
        if severity_var in metadata_df.columns:
            logging.info(f"\n4.2 Analyzing species abundance by {severity_var}")
            severity_results = analyze_species_by_variable(abundance_df, metadata_df, species_indices, severity_var)
            
            if not severity_results.empty:
                # Save results
                severity_file = tables_dir / 'species_by_severity.csv'
                severity_results.to_csv(severity_file)
                logging.info(f"Severity analysis results saved to {severity_file}")
                
                # Create plots for each species by severity if not skipped
                if not args.skip_plots:
                    for species_name, taxon in species_indices.items():
                        try:
                            output_file = figures_dir / f"{species_name.replace(' ', '_')}_by_{severity_var}.pdf"
                            plot_taxa_facet(abundance_df, metadata_df, taxon, time_var, severity_var, output_file)
                        except Exception as e:
                            logging.error(f"Error creating plot for {species_name} by {severity_var}: {str(e)}")
        
        # Analysis by symptoms
        if symptoms_var in metadata_df.columns:
            logging.info(f"\n4.3 Analyzing species abundance by {symptoms_var}")
            symptoms_results = analyze_species_by_variable(abundance_df, metadata_df, species_indices, symptoms_var)
            
            if not symptoms_results.empty:
                # Save results
                symptoms_file = tables_dir / 'species_by_symptoms.csv'
                symptoms_results.to_csv(symptoms_file)
                logging.info(f"Symptoms analysis results saved to {symptoms_file}")
                
                # Create plots for each species by symptoms if not skipped
                if not args.skip_plots:
                    for species_name, taxon in species_indices.items():
                        try:
                            output_file = figures_dir / f"{species_name.replace(' ', '_')}_by_{symptoms_var}.pdf"
                            plot_taxa_facet(abundance_df, metadata_df, taxon, time_var, symptoms_var, output_file)
                        except Exception as e:
                            logging.error(f"Error creating plot for {species_name} by {symptoms_var}: {str(e)}")
        
        # Co-occurrence analysis
        if len(species_indices) >= 2:
            logging.info("\n4.4 Analyzing co-occurrence patterns")
            co_occurrence_results = analyze_co_occurrence(abundance_df, metadata_df, species_indices)
            
            if co_occurrence_results:
                # Save results
                result_file = tables_dir / 'co_occurrence_results.csv'
                pd.DataFrame(co_occurrence_results['co_occurrence_rates']).T.to_csv(result_file)
                logging.info(f"Co-occurrence results saved to {result_file}")
                
                # Create visualizations if not skipped
                if not args.skip_plots:
                    logging.info("Creating co-occurrence visualizations")
                    plot_co_occurrence_heatmap(abundance_df, metadata_df, species_indices, figures_dir)
                    plot_scatter_correlation(abundance_df, metadata_df, species_indices, figures_dir)
        
        # Analysis by timepoint and severity/symptoms
        logging.info("\n4.5 Analyzing species abundance by severity and symptoms at each time point")
        
        # By severity
        if severity_var in metadata_df.columns:
            severity_time_results = analyze_species_by_timepoint_and_variable(
                abundance_df, metadata_df, species_indices, time_var, severity_var
            )
            
            for time_point, results_df in severity_time_results.items():
                # Save results
                result_file = tables_dir / f"{time_point}_{severity_var}_results.csv"
                results_df.to_csv(result_file)
                logging.info(f"Results for {time_point} by {severity_var} saved to {result_file}")
        
        # By symptoms
        if symptoms_var in metadata_df.columns:
            symptoms_time_results = analyze_species_by_timepoint_and_variable(
                abundance_df, metadata_df, species_indices, time_var, symptoms_var
            )
            
            for time_point, results_df in symptoms_time_results.items():
                # Save results
                result_file = tables_dir / f"{time_point}_{symptoms_var}_results.csv"
                results_df.to_csv(result_file)
                logging.info(f"Results for {time_point} by {symptoms_var} saved to {result_file}")
        
        # Analysis of Prior to Acute delta
        logging.info("\n4.6 Analyzing changes from Prior to Acute timepoints")
        if not args.skip_plots:
            delta_dfs = analyze_prior_acute_delta(abundance_df, metadata_df, species_indices, output_dir)
        
        # Create subject trajectory plots if not skipped
        if not args.skip_plots:
            logging.info("\n4.7 Creating subject trajectory plots")
            plot_subject_trajectory(abundance_df, metadata_df, species_indices, output_dir)
    
    # Include info about the normalization method used in the results
    normalization_info = {
        'method': args.normalization
    }
    if args.normalization == 'rarefaction' and args.rarefaction_depth:
        normalization_info['rarefaction_depth'] = args.rarefaction_depth
    
    # Save normalization info
    with open(output_dir / 'normalization_info.txt', 'w') as f:
        for key, value in normalization_info.items():
            f.write(f"{key}: {value}\n")
    
    logging.info(f"\nNormalization method '{args.normalization}' info saved to {output_dir / 'normalization_info.txt'}")
    
    logging.info("\nAnalysis complete!")
    
    return {
        'abundance_df': abundance_df,
        'metadata_df': metadata_df,
        'species_indices': species_indices
    }

if __name__ == "__main__":
    try:
        # Parse command line arguments
        args = parse_args()
        
        # Run the full analysis
        run_analysis_pipeline(args)
        
    except Exception as e:
        logging.error(f"Error in main execution: {str(e)}")
        traceback.print_exc()
        sys.exit(1)