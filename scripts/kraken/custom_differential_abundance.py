#!/usr/bin/env python
# scripts/kraken/custom_differential_abundance.py

import argparse
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multitest import multipletests
import logging

def setup_logger(log_file=None, log_level=logging.INFO):
    """Setup logging."""
    logger = logging.getLogger('differential_abundance')
    logger.setLevel(log_level)
    
    # Create formatter
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    
    # Create console handler
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    # Create file handler if log_file is provided
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    
    return logger

def log_print(message, level="info"):
    """Print message to console and log to logger."""
    print(message)
    logger = logging.getLogger('differential_abundance')
    getattr(logger, level)(message)

def parse_arguments():
    parser = argparse.ArgumentParser(description="Run differential abundance analysis on microbiome data")
    
    parser.add_argument(
        "--abundance-file", 
        required=True,
        help="Path to abundance file with metadata (TSV format)"
    )
    
    parser.add_argument(
        "--output-dir",
        default="results/differential_abundance",
        help="Directory for output files (default: results/differential_abundance)"
    )
    
    parser.add_argument(
        "--group-col", 
        default="Timing",
        help="Column name for grouping samples (default: Timing)"
    )
    
    parser.add_argument(
        "--log-file",
        default=None,
        help="Path to log file (default: log to console only)"
    )
    
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Logging level (default: INFO)"
    )
    
    parser.add_argument(
        "--method",
        default="kruskal",
        choices=["kruskal", "anova", "t-test"],
        help="Statistical method to use (default: kruskal)"
    )
    
    parser.add_argument(
        "--filter-groups",
        help="Comma-separated list of groups to include (default: use all groups)"
    )
    
    parser.add_argument(
        "--p-threshold", 
        type=float, 
        default=0.05, 
        help="P-value threshold for significance (default: 0.05)"
    )
    
    return parser.parse_args()

def kruskal_wallis_test(data, taxa_col, group_col, taxa_list):
    """
    Perform Kruskal-Wallis H test for each taxon.
    
    Args:
        data: DataFrame with taxa and group information
        taxa_col: Column name with taxon names
        group_col: Column name with group information
        taxa_list: List of taxa to test
        
    Returns:
        DataFrame with test results
    """
    results = []
    
    for taxon in taxa_list:
        # Group data by the group column
        groups = []
        for group_name, group_data in data.groupby(group_col):
            if taxon in group_data.columns:
                groups.append(group_data[taxon].values)
        
        if len(groups) < 2:
            continue
            
        # Perform Kruskal-Wallis H test
        h_stat, p_value = stats.kruskal(*groups)
        
        results.append({
            'Taxon': taxon,
            'H_statistic': h_stat,
            'p_value': p_value,
        })
    
    # Convert results to DataFrame
    results_df = pd.DataFrame(results)
    
    # Apply multiple testing correction
    if len(results_df) > 0:
        results_df['adjusted_p_value'] = multipletests(results_df['p_value'], method='fdr_bh')[1]
    
    return results_df

def anova_test(data, taxa_col, group_col, taxa_list):
    """
    Perform ANOVA test for each taxon.
    
    Args:
        data: DataFrame with taxa and group information
        taxa_col: Column name with taxon names
        group_col: Column name with group information
        taxa_list: List of taxa to test
        
    Returns:
        DataFrame with test results
    """
    results = []
    
    for taxon in taxa_list:
        # Create a dataframe for the current taxon
        taxon_data = pd.DataFrame({
            'abundance': data[taxon],
            'group': data[group_col]
        })
        
        # Fit the ANOVA model
        try:
            model = ols(f'abundance ~ C(group)', data=taxon_data).fit()
            anova_table = sm.stats.anova_lm(model)
            
            f_stat = anova_table.loc['C(group)', 'F']
            p_value = anova_table.loc['C(group)', 'PR(>F)']
            
            results.append({
                'Taxon': taxon,
                'F_statistic': f_stat,
                'p_value': p_value,
            })
        except Exception as e:
            log_print(f"Error in ANOVA for taxon {taxon}: {str(e)}", level="warning")
    
    # Convert results to DataFrame
    results_df = pd.DataFrame(results)
    
    # Apply multiple testing correction
    if len(results_df) > 0:
        results_df['adjusted_p_value'] = multipletests(results_df['p_value'], method='fdr_bh')[1]
    
    return results_df

def t_test(data, taxa_col, group_col, taxa_list):
    """
    Perform t-test for each taxon between two groups.
    
    Args:
        data: DataFrame with taxa and group information
        taxa_col: Column name with taxon names
        group_col: Column name with group information
        taxa_list: List of taxa to test
        
    Returns:
        DataFrame with test results
    """
    results = []
    
    # Check if there are exactly two groups
    unique_groups = data[group_col].unique()
    if len(unique_groups) != 2:
        log_print(f"Error: t-test requires exactly 2 groups, but found {len(unique_groups)}", level="error")
        return pd.DataFrame()
    
    # Get the two groups
    group1 = unique_groups[0]
    group2 = unique_groups[1]
    
    for taxon in taxa_list:
        # Get abundance values for each group
        group1_values = data[data[group_col] == group1][taxon].values
        group2_values = data[data[group_col] == group2][taxon].values
        
        # Perform t-test
        t_stat, p_value = stats.ttest_ind(group1_values, group2_values, equal_var=False)
        
        results.append({
            'Taxon': taxon,
            'Group1': group1,
            'Group2': group2,
            't_statistic': t_stat,
            'p_value': p_value,
        })
    
    # Convert results to DataFrame
    results_df = pd.DataFrame(results)
    
    # Apply multiple testing correction
    if len(results_df) > 0:
        results_df['adjusted_p_value'] = multipletests(results_df['p_value'], method='fdr_bh')[1]
    
    return results_df

def generate_boxplots(data, significant_taxa, group_col, output_dir):
    """
    Generate boxplots for significant taxa.
    
    Args:
        data: DataFrame with abundance and group information
        significant_taxa: List of significant taxa
        group_col: Column name with group information
        output_dir: Directory to save the plots
    """
    # Create directory for boxplots
    boxplot_dir = os.path.join(output_dir, 'boxplots')
    os.makedirs(boxplot_dir, exist_ok=True)
    
    for taxon in significant_taxa:
        plt.figure(figsize=(10, 6))
        sns.boxplot(x=group_col, y=taxon, data=data)
        plt.title(f'Abundance of {taxon} by {group_col}')
        plt.ylabel('Normalized Abundance')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(os.path.join(boxplot_dir, f'{taxon.replace(".", "_")}_boxplot.png'), dpi=300)
        plt.close()

def main():
    # Parse command line arguments
    args = parse_arguments()
    
    # Setup logging
    logger = setup_logger(log_file=args.log_file, log_level=getattr(logging, args.log_level))
    log_print("Starting differential abundance analysis", level="info")
    
    # Check if abundance file exists
    if not os.path.exists(args.abundance_file):
        log_print(f"Error: Abundance file not found: {args.abundance_file}", level="error")
        sys.exit(1)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Parse filter groups if provided
    filter_groups = None
    if args.filter_groups:
        filter_groups = args.filter_groups.split(',')
        log_print(f"Filtering to groups: {', '.join(filter_groups)}", level="info")
    
    # Load data
    try:
        # Load abundance data with metadata
        data = pd.read_csv(args.abundance_file, sep='\t')
        
        # Check if the group column exists
        if args.group_col not in data.columns:
            log_print(f"Error: Group column '{args.group_col}' not found in the data", level="error")
            sys.exit(1)
            
        # Filter groups if specified
        if filter_groups:
            data = data[data[args.group_col].isin(filter_groups)]
            log_print(f"Filtered to {len(data)} samples after group filtering", level="info")
        
        # Identify taxa columns (all columns that are numeric and not the group column)
        metadata_cols = ['SampleID', 'SubjectID', 'CollectionDate', 'Timing', 'Severity', 'Symptoms']
        available_metadata_cols = [col for col in metadata_cols if col in data.columns]
        
        # All columns that are not metadata are taxonomy columns
        taxonomy_cols = [col for col in data.columns if col not in available_metadata_cols]
        
        log_print(f"Loaded data with {len(data)} samples and {len(taxonomy_cols)} taxa", level="info")
        
        # Run the appropriate statistical test
        if args.method == 'kruskal':
            log_print("Running Kruskal-Wallis H test", level="info")
            results = kruskal_wallis_test(data, None, args.group_col, taxonomy_cols)
        elif args.method == 'anova':
            log_print("Running ANOVA test", level="info")
            results = anova_test(data, None, args.group_col, taxonomy_cols)
        elif args.method == 't-test':
            log_print("Running t-test", level="info")
            results = t_test(data, None, args.group_col, taxonomy_cols)
        
        # Sort results by p-value
        if len(results) > 0:
            results = results.sort_values('p_value')
            
            # Save results to file
            results_file = os.path.join(args.output_dir, f"{args.group_col}_{args.method}_results.tsv")
            results.to_csv(results_file, sep='\t', index=False)
            log_print(f"Saved results to {results_file}", level="info")
            
            # Count significant taxa
            significant_taxa = results[results['adjusted_p_value'] < args.p_threshold]['Taxon'].tolist()
            log_print(f"Found {len(significant_taxa)} significant taxa (adjusted p < {args.p_threshold})", level="info")
            
            # Generate boxplots for significant taxa
            if significant_taxa:
                log_print("Generating boxplots for significant taxa", level="info")
                generate_boxplots(data, significant_taxa, args.group_col, args.output_dir)
        else:
            log_print("No statistical test results produced", level="warning")
    
    except Exception as e:
        log_print(f"Error during differential abundance analysis: {str(e)}", level="error")
        import traceback
        log_print(traceback.format_exc(), level="debug")
        sys.exit(1)
    
    log_print(f"Results saved to {args.output_dir}", level="info")

if __name__ == "__main__":
    main()