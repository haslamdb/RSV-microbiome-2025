#!/usr/bin/env python3
"""
Calculate and analyze alpha and beta diversity of the microbiome data.

This script:
1. Loads the combined abundance table and metadata
2. Calculates alpha diversity metrics (Shannon, Simpson, Observed species)
3. Compares alpha diversity between clinical groups
4. Calculates beta diversity distance matrix
5. Performs PERMANOVA to test for community composition differences
6. Generates visualizations of diversity metrics

Usage:
    python scripts/02_calculate_diversity.py [--config CONFIG_FILE]
"""

import os
import sys
import argparse
import yaml
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# Add project root to Python path
project_root = Path(__file__).resolve().parents[1]
sys.path.append(str(project_root))

# Add tools directory to Python path
tools_dir = project_root / 'tools'
sys.path.append(str(tools_dir))

# Now import metaphlan_tools functions
from metaphlan_tools import parse_metaphlan_file, combine_samples, load_metadata

# Import our patched function
sys.path.append(str(project_root / 'scripts' / 'utils'))
from parser_patch import patched_combine_samples


# Import functions from metaphlan_tools
from metaphlan_tools import (
    load_metadata,
    calculate_alpha_diversity,
    compare_alpha_diversity,
    calculate_beta_diversity,
    perform_permanova,
    plot_alpha_diversity_boxplot
)

# Import additional visualization functions
from metaphlan_tools.stats import plot_beta_diversity_ordination


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Calculate and analyze microbiome diversity')
    parser.add_argument('--config', type=str, default='config/analysis_parameters.yml',
                       help='Path to configuration file')
    return parser.parse_args()


def main():
    """Main function to calculate diversity metrics."""
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
    
    # Create directories if they don't exist
    figures_dir.mkdir(exist_ok=True, parents=True)
    tables_dir.mkdir(exist_ok=True, parents=True)
    
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
    
    # Calculate alpha diversity
    print("Calculating alpha diversity metrics...")
    alpha_metrics = config['diversity']['alpha_metrics']
    alpha_df = calculate_alpha_diversity(abundance_df, metrics=alpha_metrics)
    
    # Save alpha diversity results
    alpha_file = tables_dir / 'alpha_diversity.csv'
    alpha_df.to_csv(alpha_file)
    print(f"Alpha diversity saved to {alpha_file}")
    
    # Compare alpha diversity between clinical groups
    group_vars = config['metadata']['group_variables']
    alpha_stats_results = {}
    
    for var in group_vars:
        if var in metadata_df.columns:
            print(f"\nAnalyzing differences in alpha diversity by {var}")
            results = compare_alpha_diversity(alpha_df, metadata_df, var)
            alpha_stats_results[var] = results
            
            # Print results
            for metric, stats in results.items():
                if stats['p-value'] is not None:
                    print(f"  {metric}: {stats['test']} p-value = {stats['p-value']:.4f}")
                else:
                    note = stats.get('note', 'Statistical test could not be performed')
                    print(f"  {metric}: {stats['test']} - {note}")
            
            # Create and save boxplot
            fig = plot_alpha_diversity_boxplot(alpha_df, metadata_df, var)
            boxplot_file = figures_dir / f"alpha_diversity_{var}.png"
            fig.savefig(boxplot_file, dpi=config['visualization']['figure_dpi'], bbox_inches='tight')
            plt.close(fig)
            print(f"  Boxplot saved to {boxplot_file}")
        else:
            print(f"Warning: Variable '{var}' not found in metadata")
    
    # Save alpha diversity comparison results
    for var, results in alpha_stats_results.items():
        # Convert nested dict to DataFrame
        results_df = pd.DataFrame({metric: {k: v for k, v in stats.items()} 
                                for metric, stats in results.items()}).T
        results_file = tables_dir / f"alpha_diversity_{var}_stats.csv"
        results_df.to_csv(results_file)
        print(f"Statistical results saved to {results_file}")
    
    # Calculate beta diversity
    print("\nCalculating beta diversity...")
    beta_metric = config['diversity']['beta_metric']
    beta_dm = calculate_beta_diversity(abundance_df, metric=beta_metric)
    
    # Perform PERMANOVA tests
    permanova_results = {}
    
    for var in group_vars:
        if var in metadata_df.columns:
            print(f"\nPerforming PERMANOVA for {var}")
            result = perform_permanova(beta_dm, metadata_df, var)
            permanova_results[var] = result
            
            # Print results
            print(f"  Test statistic: {result['test-statistic']:.4f}")
            print(f"  p-value: {result['p-value']:.4f}")
            print(f"  Sample size: {result['sample size']}")
        else:
            print(f"Warning: Variable '{var}' not found in metadata")
    
    # Save PERMANOVA results
    permanova_df = pd.DataFrame.from_dict(permanova_results, orient='index')
    permanova_file = tables_dir / 'permanova_results.csv'
    permanova_df.to_csv(permanova_file)
    print(f"\nPERMANOVA results saved to {permanova_file}")
    
    # Create ordination plots
    print("\nCreating ordination plots...")
    for var in group_vars:
        if var in metadata_df.columns:
            print(f"  Creating ordination plot for {var}")
            try:
                fig = plot_beta_diversity_ordination(beta_dm, metadata_df, var, method='PCoA')
                ordination_file = figures_dir / f"beta_diversity_pcoa_{var}.png"
                fig.savefig(ordination_file, dpi=config['visualization']['figure_dpi'], bbox_inches='tight')
                plt.close(fig)
                print(f"  Ordination plot saved to {ordination_file}")
            except Exception as e:
                print(f"  Error creating ordination plot for {var}: {str(e)}")
    
    # Time-based analysis if time variable exists
    time_var = config['metadata']['time_variable']
    if time_var in metadata_df.columns:
        # Create alpha diversity time series plots
        print(f"\nAnalyzing diversity changes over {time_var}")
        
        # Join alpha diversity with metadata for time plot
        alpha_time_df = alpha_df.join(metadata_df[[time_var]], how='inner')
        
        # Plot for each metric
        # Create boxplots for time analysis
        for metric in alpha_metrics:
            fig, ax = plt.subplots(figsize=(10, 6))
            
            # Count samples in each time point
            time_counts = alpha_time_df[time_var].value_counts()
            
            # Only include time points with enough samples
            valid_times = time_counts[time_counts >= 2].index
            
            if len(valid_times) >= 2:
                # Filter to valid time points
                plot_data = alpha_time_df[alpha_time_df[time_var].isin(valid_times)]
                
                try:
                    # Create boxplot
                    sns_plot = sns.boxplot(x=time_var, y=metric, data=plot_data, ax=ax)
                    ax.set_title(f'{metric} Diversity Over {time_var}')
                    ax.set_xlabel(time_var)
                    ax.set_ylabel(f'{metric} Diversity')
                    
                    # Rotate x-axis labels if needed
                    plt.xticks(rotation=45)
                    plt.tight_layout()
                    
                    # Save the figure
                    fig.savefig(os.path.join(output_dir, f'alpha_{metric}_{time_var}.png'), dpi=300)
                    print(f"  Time series plot for {metric} saved")
                    
                except Exception as e:
                    # Handle plotting error
                    ax.text(0.5, 0.5, f"Error creating plot: {str(e)}",
                        horizontalalignment='center', verticalalignment='center',
                        transform=ax.transAxes, fontsize=12, color='red')
                    ax.set_title(f'{metric} Diversity Over {time_var}')
                    ax.axis('off')
                    
                    # Save the error message figure
                    fig.savefig(os.path.join(output_dir, f'alpha_{metric}_{time_var}_error.png'), dpi=300)
                    print(f"  Error creating time series plot for {metric}: {str(e)}")
            else:
                # Not enough valid time points
                ax.text(0.5, 0.5, "Insufficient data points\nfor time series analysis",
                    horizontalalignment='center', verticalalignment='center',
                    transform=ax.transAxes, fontsize=12)
                ax.set_title(f'{metric} Diversity Over {time_var}')
                ax.axis('off')
                
                # Save the message figure
                fig.savefig(os.path.join(output_dir, f'alpha_{metric}_{time_var}_insufficient.png'), dpi=300)
                print(f"  Insufficient data for time series plot for {metric}")
            
            plt.close(fig)
    
    print("\nDiversity analysis complete!")


if __name__ == "__main__":
    # Import seaborn here to avoid import error in function
    import seaborn as sns
    sns.set(style="whitegrid")
    
    main()
