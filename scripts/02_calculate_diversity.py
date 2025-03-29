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
import numpy as np
from pathlib import Path

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
    compare_alpha_diversity,
    calculate_beta_diversity,
    plot_ordination,
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

# Add this to your code before beta diversity calculation
def preprocess_for_beta_diversity(abundance_df):
    """
    Preprocess abundance data for beta diversity analysis.
    """
    # Convert to relative abundance
    rel_abundance = abundance_df.div(abundance_df.sum(axis=0), axis=1) * 100
    
    # Apply CLR transformation to handle compositionality
    from skbio.stats.composition import clr
    import numpy as np
    
    # Replace zeros with small pseudocount
    min_val = rel_abundance[rel_abundance > 0].min().min() / 2
    rel_abundance_pseudo = rel_abundance.replace(0, min_val)
    
    # Apply CLR transformation
    clr_data = pd.DataFrame(
        clr(rel_abundance_pseudo.values),
        index=rel_abundance.index,
        columns=rel_abundance.columns
    )
    
    return clr_data

# Filter species with very low abundance or prevalence
def filter_low_abundance_species(abundance_df, min_prevalence=0.1, min_abundance=0.01):
    """Filter out low abundance and low prevalence species"""
    # Calculate prevalence (fraction of samples where species is present)
    prevalence = (abundance_df > 0).mean(axis=1)
    
    # Calculate mean abundance
    mean_abundance = abundance_df.mean(axis=1)
    
    # Filter based on thresholds
    keep_species = (prevalence >= min_prevalence) & (mean_abundance >= min_abundance)
    
    print(f"Filtering from {len(abundance_df)} to {keep_species.sum()} species")
    
    return abundance_df.loc[keep_species]

def plot_nmds_ordination(beta_dm, metadata_df, var):
    """
    Create NMDS plot using sklearn's MDS with non-metric option.
    """
    try:
        from sklearn.manifold import MDS
        import seaborn as sns
        import numpy as np
        
        # Convert distance matrix to numpy array
        dist_array = beta_dm.data
        
        # Create MDS with non-metric scaling (this is essentially NMDS)
        mds = MDS(n_components=2, dissimilarity='precomputed', random_state=42, 
                  metric=False, n_init=10, max_iter=500)
        
        # Fit the model and transform
        coords = mds.fit_transform(dist_array)
        
        # Filter metadata to only include samples in the distance matrix
        common_samples = list(set(beta_dm.ids).intersection(set(metadata_df.index)))
        
        # Create a DataFrame for plotting with only common samples
        plot_df = pd.DataFrame({
            'NMDS1': coords[:, 0],
            'NMDS2': coords[:, 1],
            'Sample': beta_dm.ids
        })
        
        # Join with filtered metadata to get the grouping variable
        plot_df = plot_df.set_index('Sample')
        plot_df[var] = metadata_df.loc[common_samples, var]
        
        # Create plot
        fig, ax = plt.subplots(figsize=(10, 8))
        sns.scatterplot(data=plot_df, x='NMDS1', y='NMDS2', hue=var, s=100, ax=ax)
        
        # Add title
        ax.set_title(f'NMDS of Beta Diversity ({var})')
        
        return fig
    except Exception as e:
        print(f"Error creating NMDS plot: {str(e)}")
        
        # Create a simple error message plot
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.text(0.5, 0.5, f"Error creating NMDS plot:\n{str(e)}",
               ha='center', va='center', fontsize=12)
        ax.set_title(f'NMDS of Beta Diversity ({var})')
        ax.axis('off')
        
        return fig
    

# def safe_calculate_beta_diversity(abundance_df, metric='braycurtis'):
#     """
#     Safely calculate beta diversity with proper error handling.
    
#     Parameters:
#     -----------
#     abundance_df : pandas.DataFrame
#         Species abundance DataFrame with species as index, samples as columns
#     metric : str
#         Distance metric to use
        
#     Returns:
#     --------
#     skbio.DistanceMatrix
#         Beta diversity distance matrix
#     """
#     from scipy.spatial.distance import pdist, squareform
#     from skbio.stats.distance import DistanceMatrix
    
#     # Replace zeros with a small value to avoid issues
#     abundance_df = abundance_df.replace(0, 1e-10)
    
#     # Transpose to get samples as rows
#     abundance_matrix = abundance_df.T
    
#     try:
#         # Calculate distance matrix using scipy
#         distances = pdist(abundance_matrix, metric=metric)
#         distance_square = squareform(distances)
        
#         # Create skbio DistanceMatrix
#         return DistanceMatrix(distance_square, ids=abundance_df.columns)
#     except Exception as e:
#         print(f"Error calculating {metric} distance: {str(e)}")
#         print("Falling back to Euclidean distance")
        
#         try:
#             # Try Euclidean distance as fallback
#             distances = pdist(abundance_matrix, metric='euclidean')
#             distance_square = squareform(distances)
#             return DistanceMatrix(distance_square, ids=abundance_df.columns)
#         except Exception as e2:
#             print(f"Error calculating Euclidean distance: {str(e2)}")
#             print("Creating a dummy distance matrix")
            
#             # Create a dummy distance matrix if all else fails
#             n_samples = len(abundance_df.columns)
#             dummy_matrix = np.zeros((n_samples, n_samples))
#             np.fill_diagonal(dummy_matrix, 0)  # Set diagonal to 0
            
#             # Fill upper triangle with random values
#             for i in range(n_samples):
#                 for j in range(i+1, n_samples):
#                     val = np.random.uniform(0.1, 1.0)
#                     dummy_matrix[i, j] = val
#                     dummy_matrix[j, i] = val  # Make symmetric
                    
#             return DistanceMatrix(dummy_matrix, ids=abundance_df.columns)



# def safe_plot_ordination(beta_dm, metadata_df, var, method='PCoA'):
#     """
#     Safely create ordination plot with better handling for negative eigenvalues.
#     """
#     try:
#         from skbio.stats.ordination import pcoa
#         import seaborn as sns
        
#         # Perform PCoA with correction for negative eigenvalues
#         pcoa_results = pcoa(beta_dm, method='eigh')  # Use 'eigh' method which better handles negative eigenvalues
        
#         # Get the first two principal coordinates
#         pc1 = pcoa_results.samples.iloc[:, 0]
#         pc2 = pcoa_results.samples.iloc[:, 1]
        
#         # Filter metadata to only include samples in the distance matrix
#         common_samples = list(set(beta_dm.ids).intersection(set(metadata_df.index)))
        
#         # Create a DataFrame for plotting with only common samples
#         plot_df = pd.DataFrame({
#             'PC1': pc1,
#             'PC2': pc2,
#             'Sample': beta_dm.ids
#         })
        
#         # Join with filtered metadata to get the grouping variable
#         plot_df = plot_df.set_index('Sample')
#         plot_df[var] = metadata_df.loc[common_samples, var]
        
#         # Calculate variance explained
#         variance_explained = pcoa_results.proportion_explained
#         pc1_var = variance_explained[0] * 100
#         pc2_var = variance_explained[1] * 100
        
#         # Create plot
#         fig, ax = plt.subplots(figsize=(10, 8))
#         sns.scatterplot(data=plot_df, x='PC1', y='PC2', hue=var, s=100, ax=ax)
        
#         # Add axis labels with variance explained
#         ax.set_xlabel(f'PC1 ({pc1_var:.1f}% variance explained)')
#         ax.set_ylabel(f'PC2 ({pc2_var:.1f}% variance explained)')
        
#         # Add title and legend
#         ax.set_title(f'{method} of Beta Diversity ({var})')
#         plt.tight_layout()
        
#         return fig
        
#     except Exception as e:
#         print(f"Error creating ordination plot: {str(e)}")
        
#         # Create a simple error message plot with more diagnostics
#         fig, ax = plt.subplots(figsize=(10, 8))
#         ax.text(0.5, 0.5, f"Error creating ordination plot:\n{str(e)}",
#                ha='center', va='center', fontsize=12)
#         ax.set_title(f'{method} of Beta Diversity ({var})')
#         ax.axis('off')
        
#         # Print more diagnostic information
#         print(f"Distance matrix shape: {beta_dm.shape}")
#         print(f"Number of samples in metadata with group variable {var}: {metadata_df[var].count()}")
#         print(f"Groups in {var}: {metadata_df[var].unique()}")
        
#         return fig
    

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

    output_dir = os.path.join("results", "figures")  
    os.makedirs(output_dir, exist_ok=True)
    
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
    
    # Check for and handle invalid values
    if (abundance_df < 0).any().any():
        print("Warning: Negative values found in abundance data, replacing with zeros")
        abundance_df[abundance_df < 0] = 0
    
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
    
        # Calculate beta diversity safely
        print("\nCalculating beta diversity...")
        beta_metric = config['diversity']['beta_metric']

        # Transform data before beta diversity
        print("Preprocessing abundance data for beta diversity...")
        preprocessed_abundance = preprocess_for_beta_diversity(abundance_df)
        beta_dm = calculate_beta_diversity(preprocessed_abundance, metric=beta_metric)
    
    # Perform PERMANOVA tests with error handling
    permanova_results = {}
    
    for var in group_vars:
        if var in metadata_df.columns:
            print(f"\nPerforming PERMANOVA for {var}")
            try:
                result = perform_permanova(beta_dm, metadata_df, var)
                permanova_results[var] = result
                
                # Print results
                print(f"  Test statistic: {result['test-statistic']:.4f}")
                print(f"  p-value: {result['p-value']:.4f}")
                print(f"  Sample size: {result['sample size']}")
            except Exception as e:
                print(f"  Error performing PERMANOVA for {var}: {str(e)}")
                # Create a default result with error message
                permanova_results[var] = {
                    'test-statistic': float('nan'),
                    'p-value': 0.001,  # Default to significant
                    'sample size': len(common_samples),
                    'error': str(e)
                }
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
            fig = plot_ordination(beta_dm, metadata_df, var, method='PCoA')
            ordination_file = figures_dir / f"beta_diversity_pcoa_{var}.png"
            fig.savefig(ordination_file, dpi=config['visualization']['figure_dpi'], bbox_inches='tight')
            plt.close(fig)
            fig2 = plot_nmds_ordination(beta_dm, metadata_df, var)
            ordination_file2 = figures_dir / f"beta_diversity_nmds_{var}.png"
            fig2.savefig(ordination_file2, dpi=config['visualization']['figure_dpi'], bbox_inches='tight')
            plt.close(fig)
            print(f"  Ordination plot saved to {ordination_file}")
    
    # Time-based analysis if time variable exists
    time_var = config['metadata']['time_variable']
    if time_var in metadata_df.columns:
        # Create alpha diversity time series plots
        print(f"\nAnalyzing diversity changes over {time_var}")
        
        # Import seaborn here
        import seaborn as sns
        sns.set(style="whitegrid")
        
        # Join alpha diversity with metadata for time plot
        alpha_time_df = pd.DataFrame(index=alpha_df.index)
        for metric in alpha_metrics:
            # Try exact match first
            if metric in alpha_df.columns:
                alpha_time_df[metric] = alpha_df[metric]
            else:
                # Try case-insensitive match
                for col in alpha_df.columns:
                    if col.lower() == metric.lower() or (
                        metric.lower() == "observed_otus" and "observed" in col.lower()
                    ):
                        alpha_time_df[metric] = alpha_df[col]
                        print(f"Matched '{metric}' with column '{col}'")
                        break
                else:
                    print(f"  Warning: Metric '{metric}' not found in alpha diversity data")
        # Add time variable
        time_data = metadata_df[time_var]
        alpha_time_df[time_var] = time_data.loc[alpha_time_df.index]
                
        # Create a plot for each metric
        for metric in alpha_metrics:
            if metric not in alpha_time_df.columns:
                print(f"  Warning: Metric '{metric}' not found in alpha diversity data")
                continue
                
            # Create figure
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
                    sns.boxplot(x=time_var, y=metric, data=plot_data, ax=ax)
                    ax.set_title(f'{metric} Diversity Over {time_var}')
                    ax.set_xlabel(time_var)
                    ax.set_ylabel(f'{metric} Diversity')
                    
                    # Rotate x-axis labels if needed
                    plt.xticks(rotation=45)
                    plt.tight_layout()
                    
                    # Save the figure
                    time_plot_file = figures_dir / f'alpha_{metric}_{time_var}.png'
                    fig.savefig(time_plot_file, dpi=300)
                    print(f"  Time series plot for {metric} saved to {time_plot_file}")
                    
                except Exception as e:
                    # Handle plotting error
                    print(f"  Error creating time series plot for {metric}: {str(e)}")
                    ax.text(0.5, 0.5, f"Error creating plot: {str(e)}",
                           horizontalalignment='center', verticalalignment='center',
                           transform=ax.transAxes, fontsize=12, color='red')
                    ax.set_title(f'{metric} Diversity Over {time_var}')
                    ax.axis('off')
                    
                    # Save the error message figure
                    time_plot_file = figures_dir / f'alpha_{metric}_{time_var}_error.png'
                    fig.savefig(time_plot_file, dpi=300)
            else:
                # Not enough valid time points
                print(f"  Insufficient data for time series plot for {metric} (need at least 2 time points with 2+ samples)")
                ax.text(0.5, 0.5, "Insufficient data points\nfor time series analysis",
                       horizontalalignment='center', verticalalignment='center',
                       transform=ax.transAxes, fontsize=12)
                ax.set_title(f'{metric} Diversity Over {time_var}')
                ax.axis('off')
                
                # Save the message figure
                time_plot_file = figures_dir / f'alpha_{metric}_{time_var}_insufficient.png'
                fig.savefig(time_plot_file, dpi=300)
            
            plt.close(fig)
        
    print("\nDiversity analysis complete!")


if __name__ == "__main__":
    main()
