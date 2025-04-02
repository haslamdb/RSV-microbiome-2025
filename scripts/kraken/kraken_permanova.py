#!/usr/bin/env python
# scripts/kraken/kraken_permanova.py

import argparse
import os
import sys
import logging
import pandas as pd
import numpy as np
import tempfile
import scipy.spatial.distance as ssd
from skbio.stats.distance import DistanceMatrix
import skbio.stats.distance

# Add kraken_tools to path
kraken_tools_path = os.path.expanduser("~/Documents/Code/kraken_tools")
if os.path.exists(kraken_tools_path):
    sys.path.append(kraken_tools_path)
else:
    print(f"Error: kraken_tools directory not found at {kraken_tools_path}")
    sys.exit(1)

try:
    from kraken_tools.analysis.permanova import run_permanova_analysis
    from kraken_tools.logger import setup_logger, log_print
except ImportError as e:
    print(f"Error importing from kraken_tools: {e}")
    print("Make sure kraken_tools is installed and accessible")
    sys.exit(1)
    
def hellinger_transform(df):
    """Apply Hellinger transformation to abundance data."""
    # Calculate sample sums (column sums)
    sample_sums = df.sum(axis=0)
    # Divide each value by its sample sum
    df_rel = df.div(sample_sums, axis=1)
    # Apply square root
    return np.sqrt(df_rel)

def run_permanova_direct(abundance_file, metadata_file, output_dir, categorical_vars=None,
                         distance_metric="bray", transform="hellinger",
                         permutations=999, min_group_size=3):
    """
    Run PERMANOVA analysis directly on abundance data.
    
    This implementation computes distance between samples correctly.
    
    Parameters:
    -----------
    abundance_file : str
        Path to abundance file (taxa as rows, samples as columns)
    metadata_file : str
        Path to metadata file (samples as rows)
    output_dir : str
        Directory for output files
    categorical_vars : str or list
        Categorical variables to test
    distance_metric : str
        Distance metric to use (bray, jaccard, euclidean)
    transform : str
        Transformation to apply (clr, hellinger, log, none)
    permutations : int
        Number of permutations for PERMANOVA
    min_group_size : int
        Minimum group size to include
    
    Returns:
    --------
    dict
        Dictionary of PERMANOVA results
    """
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Load abundance data
    abundance_df = pd.read_csv(abundance_file, sep='\t', index_col=0)
    log_print(f"Loaded abundance data: {abundance_df.shape[0]} taxa, {abundance_df.shape[1]} samples", level="info")
    
    # Load metadata
    metadata_df = pd.read_csv(metadata_file)
    if 'SampleID' in metadata_df.columns:
        metadata_df.set_index('SampleID', inplace=True)
    else:
        log_print("Warning: No SampleID column found in metadata, using first column as index", level="warning")
        metadata_df.set_index(metadata_df.columns[0], inplace=True)
    
    log_print(f"Loaded metadata: {metadata_df.shape[0]} samples, {metadata_df.shape[1]} variables", level="info")
    
    # Make sure indices are strings for consistent matching
    abundance_df.columns = abundance_df.columns.astype(str)
    metadata_df.index = metadata_df.index.astype(str)
    
    # Find common samples
    common_samples = list(set(abundance_df.columns).intersection(set(metadata_df.index)))
    log_print(f"Found {len(common_samples)} common samples between abundance and metadata", level="info")
    
    if len(common_samples) < 10:
        log_print("Error: Not enough common samples found", level="error")
        return None
    
    # Filter to common samples
    abundance_df = abundance_df[common_samples]
    metadata_df = metadata_df.loc[common_samples]
    
    # Apply transformation
    if transform.lower() == 'hellinger':
        log_print("Applying Hellinger transformation to abundance data", level="info")
        abundance_df = hellinger_transform(abundance_df)
    elif transform.lower() == 'clr':
        log_print("Applying CLR transformation to abundance data", level="info")
        # Add small pseudocount to zeros
        df_pseudo = abundance_df.replace(0, np.nextafter(0, 1))
        # Log transform
        df_log = np.log(df_pseudo)
        # Subtract row-wise mean (CLR transformation)
        abundance_df = df_log.subtract(df_log.mean(axis=0), axis=1)
    elif transform.lower() == 'log':
        log_print("Applying log transformation to abundance data", level="info")
        abundance_df = np.log1p(abundance_df)
    elif transform.lower() == 'none':
        log_print("No transformation applied to abundance data", level="info")
    else:
        log_print(f"Unknown transformation: {transform}. Using no transformation.", level="warning")
    
    # Transpose to get samples as rows for distance calculation
    abundance_df_t = abundance_df.T
    
    # Calculate distance matrix between samples (not taxa!)
    if distance_metric.lower() == 'bray':
        log_print("Calculating Bray-Curtis distance between samples", level="info")
        distances = ssd.pdist(abundance_df_t.values, metric='braycurtis')
    elif distance_metric.lower() == 'jaccard':
        log_print("Calculating Jaccard distance between samples", level="info")
        distances = ssd.pdist(abundance_df_t.values, metric='jaccard')
    elif distance_metric.lower() == 'euclidean':
        log_print("Calculating Euclidean distance between samples", level="info")
        distances = ssd.pdist(abundance_df_t.values, metric='euclidean')
    else:
        log_print(f"Unknown distance metric: {distance_metric}. Using Euclidean.", level="warning")
        distances = ssd.pdist(abundance_df_t.values, metric='euclidean')
    
    # Create DistanceMatrix object
    dm = DistanceMatrix(distances, ids=abundance_df_t.index)
    
    # Split categorical variables
    if isinstance(categorical_vars, str):
        cat_vars = [v.strip() for v in categorical_vars.split(',')]
    else:
        cat_vars = categorical_vars
    
    log_print(f"Running PERMANOVA for: {', '.join(cat_vars)}", level="info")
    
    # Dictionary to store results
    results = {}
    
    # Run PERMANOVA for each categorical variable
    for var in cat_vars:
        if var not in metadata_df.columns:
            log_print(f"Variable '{var}' not found in metadata", level="warning")
            continue
        
        # Get grouping variable
        grouping = metadata_df[var].astype(str)
        
        # Skip if too few unique groups
        if len(grouping.unique()) < 2:
            log_print(f"Skipping {var} - fewer than 2 unique groups", level="warning")
            continue
        
        # Run PERMANOVA
        log_print(f"Running PERMANOVA for {var}", level="info")
        permanova_result = skbio.stats.distance.permanova(dm, grouping, permutations=permutations)
        
        # Debug the permanova result
        log_print(f"PERMANOVA result keys: {list(permanova_result.keys())}", level="debug")
        
        # Get p-value and test statistic
        p_value = permanova_result['p-value']
        test_stat = permanova_result['test statistic']
        
        # Calculate R²
        # different versions of scikit-bio might have different keys
        if 'R2' in permanova_result:
            r_squared = permanova_result['R2']
        else:
            # Approximate R² based on the F-statistic formula: F = (R²/(k-1))/((1-R²)/(n-k))
            # Using n = 'sample size', and k = 'number of groups'
            # Solving for R²: R² = F * (k-1) / (F * (k-1) + (n-k))
            n = permanova_result['sample size']
            k = permanova_result['number of groups']
            f_stat = permanova_result['test statistic']
            
            if k > 1 and n > k:
                numerator = f_stat * (k - 1)
                denominator = numerator + (n - k)
                r_squared = numerator / denominator if denominator > 0 else 0.0
            else:
                r_squared = 0.0
        
        # Store results
        results[var] = {
            'p_value': p_value,
            'test_statistic': test_stat,
            'R2': r_squared,
            'n_samples': len(grouping)
        }
        
        # Print results
        significance = "significant" if p_value < 0.05 else "not significant"
        log_print(f"PERMANOVA result for {var}: R2={r_squared:.4f}, p-value={p_value:.4f} ({significance})", level="info")
    
    # Save results to file
    results_df = pd.DataFrame({
        'Variable': [var for var in results],
        'R2': [results[var]['R2'] for var in results],
        'p_value': [results[var]['p_value'] for var in results],
        'Significant': ['Yes' if results[var]['p_value'] < 0.05 else 'No' for var in results],
        'n_samples': [results[var]['n_samples'] for var in results]
    })
    
    # Sort by p-value
    results_df = results_df.sort_values('p_value')
    
    # Save to CSV file
    results_file = os.path.join(output_dir, 'permanova_results.csv')
    results_df.to_csv(results_file, index=False)
    log_print(f"Results saved to {results_file}", level="info")
    
    # Create a summary file in markdown format
    summary_file = os.path.join(output_dir, 'permanova_summary.md')
    
    with open(summary_file, 'w') as f:
        f.write(f"# PERMANOVA Analysis Results Summary\n\n")
        f.write(f"## Analysis Parameters\n")
        f.write(f"- **Transformation:** {transform}\n")
        f.write(f"- **Distance Metric:** {distance_metric}\n")
        f.write(f"- **Number of Permutations:** {permutations}\n")
        f.write(f"- **Number of Samples:** {len(abundance_df.columns)}\n")
        f.write(f"- **Number of Taxa:** {len(abundance_df.index)}\n\n")
        
        f.write(f"## Results Table\n\n")
        f.write("| Variable | R² | p-value | Significant | n_samples |\n")
        f.write("|----------|-----|---------|-------------|----------|\n")
        
        for _, row in results_df.iterrows():
            var = row['Variable']
            r2 = f"{row['R2']:.3f}"
            p = f"{row['p_value']:.3f}"
            sig = row['Significant']
            n = str(row['n_samples'])
            f.write(f"| {var} | {r2} | {p} | {sig} | {n} |\n")
        
        f.write("\n## Interpretation\n\n")
        
        # Count significant results
        n_sig = sum(1 for var in results if results[var]['p_value'] < 0.05)
        
        if n_sig == 0:
            f.write("No significant associations were found between the microbiome composition and the tested variables.\n")
        else:
            sig_vars = [var for var in results if results[var]['p_value'] < 0.05]
            f.write(f"Significant associations were found for {n_sig} variables: {', '.join(sig_vars)}.\n")
            f.write("This indicates that these factors are associated with differences in the microbiome composition.\n\n")
            
            # Add effect size interpretation
            f.write("### Effect Sizes (R²)\n")
            f.write("The R² value represents the proportion of variation in community structure explained by each variable:\n")
            
            for var in sig_vars:
                r2 = results[var]['R2']
                effect = "small" if r2 < 0.05 else "medium" if r2 < 0.15 else "large"
                f.write(f"- **{var}**: R² = {r2:.3f} ({effect} effect size, explains {r2*100:.1f}% of variation)\n")
            
    log_print(f"Summary report saved to {summary_file}", level="info")
    
    return results

def check_and_transpose_data(input_file, logger=None):
    """
    Check the format of the abundance data and transpose if necessary.
    
    Parameters:
    -----------
    input_file : str
        Path to the input abundance file
    logger : function, optional
        Logging function to use
        
    Returns:
    --------
    str
        Path to a file containing properly formatted data (taxa as rows, samples as columns)
    """
    def log(message, level="info"):
        if logger:
            logger(message, level=level)
        else:
            print(message)
    
    # First read data without setting the index to determine its format
    log(f"Checking format of input file: {input_file}", level="info")
    raw_data = pd.read_csv(input_file, sep='\t')
    
    # Check if the data contains metadata columns
    metadata_cols = ['SampleID', 'SubjectID', 'CollectionDate', 'Timing', 'Severity', 'Symptoms']
    has_metadata_cols = any(col in raw_data.columns for col in metadata_cols)
    
    if has_metadata_cols:
        log("Detected metadata columns in abundance file. Data appears to be in samples-as-rows format.", level="info")
        log("Transposing data to get taxa-as-rows format required for analysis...", level="info")
        
        # Identify which metadata columns are present
        present_metadata_cols = [col for col in metadata_cols if col in raw_data.columns]
        
        # Create a new dataframe to preserve original data
        abundance_metadata_df = raw_data.copy()
        
        # Make sure we have a SampleID column, which we need for matching with metadata
        if 'SampleID' not in abundance_metadata_df.columns:
            log("WARNING: No SampleID column found. Using index as SampleID.", level="warning")
            abundance_metadata_df['SampleID'] = abundance_metadata_df.index
        
        # Extract abundance data (all columns not in metadata)
        abundance_cols = [col for col in abundance_metadata_df.columns if col not in present_metadata_cols]
        
        # Set SampleID as index before transposing
        abundance_metadata_df.set_index('SampleID', inplace=True)
        
        # Transpose only the abundance data (taxa as rows, samples as columns)
        transposed = abundance_metadata_df[abundance_cols].T
        
        # The index of the transposed dataframe becomes the columns (taxa names)
        # and the columns become the index (sample IDs)
        transposed.index.name = 'Species'
        
        # Generate a small sample of the transposed data for debugging
        sample_data = transposed.iloc[:min(5, transposed.shape[0]), :min(5, transposed.shape[1])]
        log(f"Sample of transposed data:\n{sample_data}", level="debug")
        
        # Create a temporary file for the transposed data
        temp_file = tempfile.NamedTemporaryFile(delete=False, suffix='.tsv')
        temp_filename = temp_file.name
        temp_file.close()
        
        # Save the transposed data
        transposed.to_csv(temp_filename, sep='\t')
        log(f"Successfully transposed data: now has {len(transposed.index)} taxa as rows and {len(transposed.columns)} samples as columns", level="info")
        log(f"Transposed data saved to temporary file: {temp_filename}", level="debug")
        
        return temp_filename
    else:
        # Data is likely already in the correct format or has a non-standard structure
        abundance_df = pd.read_csv(input_file, sep='\t', index_col=0)
        
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
            log(f"Data appears to be in the correct format with {len(abundance_df.index)} taxa and {len(abundance_df.columns)} samples", level="info")
        else:
            log("WARNING: Data format unclear. Will attempt to proceed with the data as is, but results may not be as expected.", level="warning")
            log("Ideal format: Taxa as rows, samples as columns. Current shape:" + str(abundance_df.shape), level="warning")
        
        return input_file

def parse_arguments():
    parser = argparse.ArgumentParser(description="Run PERMANOVA analysis on Kraken2/Bracken data")
    
    parser.add_argument(
        "--abundance-file", 
        required=True,
        help="Path to abundance file (from process_kraken_data.py)"
    )
    
    parser.add_argument(
        "--metadata", 
        required=True,
        help="Path to metadata CSV file"
    )
    
    parser.add_argument(
        "--output-dir",
        default="results/permanova",
        help="Directory for output files (default: results/permanova)"
    )
    
    parser.add_argument(
        "--categorical-vars", 
        default="Timing,Severity,Symptoms",
        help="Comma-separated list of categorical variables to test (default: Timing,Severity,Symptoms)"
    )
    
    parser.add_argument(
        "--distance-metric", 
        default="bray", 
        choices=["bray", "jaccard", "euclidean"],
        help="Distance metric to use (default: bray)"
    )
    
    parser.add_argument(
        "--transform", 
        default="clr", 
        choices=["clr", "hellinger", "log", "none"],
        help="Transformation to apply to abundance data (default: clr)"
    )
    
    parser.add_argument(
        "--permutations", 
        type=int, 
        default=999, 
        help="Number of permutations for significance testing (default: 999)"
    )
    
    parser.add_argument(
        "--min-group-size", 
        type=int, 
        default=3, 
        help="Minimum number of samples per group (default: 3)"
    )
    
    parser.add_argument(
        "--make-pcoa", 
        action="store_true", 
        default=True, 
        help="Generate PCoA plots (default: True)"
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
    
    return parser.parse_args()

def main():
    # Parse command line arguments
    args = parse_arguments()
    
    # Setup logging
    logger = setup_logger(log_file=args.log_file, log_level=getattr(logging, args.log_level))
    log_print("Starting PERMANOVA analysis", level="info")
    
    # Check if abundance file exists
    if not os.path.exists(args.abundance_file):
        log_print(f"Error: Abundance file not found: {args.abundance_file}", level="error")
        sys.exit(1)
    
    # Check if metadata file exists
    if not os.path.exists(args.metadata):
        log_print(f"Error: Metadata file not found: {args.metadata}", level="error")
        sys.exit(1)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Parse categorical variables
    categorical_vars = args.categorical_vars.split(',')
    log_print(f"Testing categorical variables: {', '.join(categorical_vars)}", level="info")
    
    # Check metadata first to see sample IDs
    try:
        metadata_df = pd.read_csv(args.metadata)
        if 'SampleID' in metadata_df.columns:
            log_print(f"Metadata file contains {len(metadata_df)} samples", level="info")
            log_print(f"First few sample IDs in metadata: {metadata_df['SampleID'][:5].tolist()}", level="debug")
        else:
            log_print("WARNING: Metadata file does not contain a SampleID column", level="warning")
    except Exception as e:
        log_print(f"Error reading metadata file: {str(e)}", level="error")
    
    # Check and potentially transpose the abundance data if needed
    abundance_file = check_and_transpose_data(args.abundance_file, logger=log_print)
    
    # Run PERMANOVA analysis using our direct implementation
    try:
        log_print("Using direct PERMANOVA implementation to ensure correct sample distances", level="info")
        results = run_permanova_direct(
            abundance_file=abundance_file,
            metadata_file=args.metadata,
            output_dir=args.output_dir,
            categorical_vars=args.categorical_vars,
            distance_metric=args.distance_metric,
            transform=args.transform,
            permutations=args.permutations,
            min_group_size=args.min_group_size
        )
        
        if results:
            log_print("PERMANOVA analysis completed successfully", level="info")
            
            # Print summary of results
            for var, result in results.items():
                if 'p_value' in result:
                    sig = "significant" if result['p_value'] < 0.05 else "not significant"
                    r2_val = result.get('r2', result.get('R2', 'NA'))
                    if isinstance(r2_val, str):
                        r2_str = f"R2={r2_val}"
                    else:
                        r2_str = f"R2={r2_val:.3f}"
                    log_print(f"Variable {var}: {r2_str}, p-value={result['p_value']:.4f} ({sig})", level="info")
        else:
            log_print("Warning: No results produced from PERMANOVA analysis", level="warning")
    
    except Exception as e:
        log_print(f"Error during PERMANOVA analysis: {str(e)}", level="error")
        import traceback
        log_print(traceback.format_exc(), level="debug")
        sys.exit(1)
    
    log_print(f"Results saved to {args.output_dir}", level="info")
    log_print(f"Summary report available at {os.path.join(args.output_dir, 'permanova_summary.md')}", level="info")
    
    # Clean up temporary files if they were created
    if abundance_file != args.abundance_file and os.path.exists(abundance_file):
        try:
            os.remove(abundance_file)
            log_print(f"Cleaned up temporary abundance file", level="debug")
        except Exception as e:
            log_print(f"Warning: Could not remove temporary file {abundance_file}: {str(e)}", level="warning")

if __name__ == "__main__":
    main()