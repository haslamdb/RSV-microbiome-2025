#!/usr/bin/env python
# scripts/kraken/process_sylph_data.py

import argparse
import os
import sys
import yaml
import pandas as pd
import glob
from pathlib import Path
import logging

def setup_logger(log_file=None, log_level=logging.INFO):
    """Set up logger for the script."""
    logger = logging.getLogger('sylph_processor')
    logger.setLevel(log_level)
    
    # Create formatter
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    
    # Create console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(log_level)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    # Create file handler if log_file specified
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(log_level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    
    return logger

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process Sylph taxonomic profiling data")
    
    parser.add_argument(
        "--config", 
        required=True,
        help="Path to configuration file (YAML)"
    )
    
    parser.add_argument(
        "--input-dir", 
        default=None,
        help="Directory containing Sylph profile files (default: data/SylphProfiles from project root)"
    )
    
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Directory for output files (default: results/kraken_analysis from project root)"
    )
    
    parser.add_argument(
        "--kingdom",
        choices=["bacteria", "fungi", "viruses", "all"],
        default="bacteria",
        help="Kingdom to process (default: bacteria)"
    )
    
    parser.add_argument(
        "--min-abundance",
        type=float,
        default=None,
        help="Minimum relative abundance threshold (default: from config)"
    )
    
    parser.add_argument(
        "--min-prevalence",
        type=float,
        default=None,
        help="Minimum prevalence threshold (default: from config)"
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

def load_config(config_path):
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
    return config

def get_sylph_files(input_dir, kingdom="bacteria"):
    """Get Sylph profile files for the specified kingdom."""
    if kingdom == "all":
        file_pattern = f"{input_dir}/*_profiled_*.tsv"
    else:
        file_pattern = f"{input_dir}/*_profiled_{kingdom}.tsv"
    
    files = glob.glob(file_pattern)
    return sorted(files)

def extract_sample_id(filepath, delimiter="-DNA"):
    """Extract sample ID from filename."""
    # Get the filename without directory
    filename = os.path.basename(filepath)
    # Split by delimiter and get the first part
    sample_id = filename.split(delimiter)[0] + delimiter
    return sample_id

def process_sylph_file(filepath, logger):
    """Process a single Sylph profile file."""
    try:
        # Read the file
        df = pd.read_csv(filepath, sep='\t')
        
        # Check if it's empty
        if df.empty:
            logger.warning(f"Empty file: {filepath}")
            return None
        
        # Extract sample ID from filename
        sample_id = extract_sample_id(filepath)
        
        # Process the data
        # Rename columns to match expected format
        if 'taxid' in df.columns and 'relabund' in df.columns:
            # Create a new DataFrame with taxid and abundance
            result_df = df[['taxid', 'name', 'relabund']].copy()
            result_df.rename(columns={'relabund': sample_id}, inplace=True)
            return result_df
        else:
            logger.warning(f"Required columns not found in {filepath}")
            return None
    
    except Exception as e:
        logger.error(f"Error processing {filepath}: {str(e)}")
        return None

def merge_abundance_data(file_dfs, logger):
    """Merge abundance data from multiple files."""
    if not file_dfs:
        logger.error("No valid data frames to merge")
        return None
    
    # Start with the first dataframe
    merged_df = file_dfs[0]
    
    # Merge with the rest based on taxid and name
    for df in file_dfs[1:]:
        if df is not None:
            merged_df = pd.merge(
                merged_df, 
                df, 
                on=['taxid', 'name'], 
                how='outer'
            )
    
    # Fill NaN values with 0
    numeric_cols = [col for col in merged_df.columns if col not in ['taxid', 'name']]
    merged_df[numeric_cols] = merged_df[numeric_cols].fillna(0)
    
    # Set taxid and name as index
    merged_df.set_index(['taxid', 'name'], inplace=True)
    
    return merged_df

def filter_abundance_data(abundance_df, min_abundance, min_prevalence, logger):
    """Filter abundance data based on minimum abundance and prevalence."""
    logger.info(f"Filtering with min_abundance={min_abundance}, min_prevalence={min_prevalence}")
    
    # Get only the abundance columns (excluding index)
    abundance_cols = abundance_df.columns
    
    # Calculate metrics for filtering
    # Maximum relative abundance across samples
    max_abundance = abundance_df[abundance_cols].max(axis=1)
    
    # Prevalence (proportion of samples where taxon is present)
    prevalence = (abundance_df[abundance_cols] > 0).mean(axis=1)
    
    # Apply filters
    keep_taxa = (max_abundance >= min_abundance) & (prevalence >= min_prevalence)
    filtered_df = abundance_df.loc[keep_taxa]
    
    logger.info(f"Filtered from {len(abundance_df)} to {len(filtered_df)} taxa")
    
    return filtered_df

def main():
    # Parse command line arguments
    args = parse_arguments()
    
    # Setup logging
    logger = setup_logger(
        log_file=args.log_file, 
        log_level=getattr(logging, args.log_level)
    )
    
    logger.info("Starting Sylph data processing")
    
    # Load configuration
    try:
        config = load_config(args.config)
        logger.info(f"Loaded configuration from {args.config}")
    except Exception as e:
        logger.error(f"Error loading configuration: {str(e)}")
        sys.exit(1)
    
    # Set input and output directories
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.dirname(os.path.dirname(script_dir))
    
    input_dir = args.input_dir
    if input_dir is None:
        input_dir = os.path.join(project_dir, "data", "SylphProfiles")
    
    output_dir = args.output_dir
    if output_dir is None:
        output_dir = os.path.join(project_dir, "results", "kraken_analysis")
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Get filtering parameters
    min_abundance = args.min_abundance
    if min_abundance is None:
        min_abundance = config["differential_abundance"]["min_abundance"]
    
    min_prevalence = args.min_prevalence
    if min_prevalence is None:
        min_prevalence = config["differential_abundance"]["min_prevalence"]
    
    # Get Sylph files
    sylph_files = get_sylph_files(input_dir, args.kingdom)
    logger.info(f"Found {len(sylph_files)} Sylph files for kingdom: {args.kingdom}")
    
    if not sylph_files:
        logger.error(f"No Sylph files found in {input_dir} for kingdom: {args.kingdom}")
        sys.exit(1)
    
    # Process each file
    file_dfs = []
    for sylph_file in sylph_files:
        logger.debug(f"Processing {sylph_file}")
        df = process_sylph_file(sylph_file, logger)
        if df is not None:
            file_dfs.append(df)
    
    # Merge abundance data
    logger.info("Merging abundance data from all files")
    merged_df = merge_abundance_data(file_dfs, logger)
    
    if merged_df is None:
        logger.error("Failed to merge abundance data")
        sys.exit(1)
    
    # Save raw abundance data
    raw_abundance_file = os.path.join(output_dir, f"raw_{args.kingdom}_abundance.tsv")
    merged_df.to_csv(raw_abundance_file, sep='\t')
    logger.info(f"Raw abundance data saved to {raw_abundance_file}")
    
    # Filter abundance data
    logger.info("Filtering abundance data")
    filtered_df = filter_abundance_data(merged_df, min_abundance, min_prevalence, logger)
    
    # Save filtered abundance data
    filtered_abundance_file = os.path.join(output_dir, f"filtered_{args.kingdom}_abundance.tsv")
    filtered_df.to_csv(filtered_abundance_file, sep='\t')
    logger.info(f"Filtered abundance data saved to {filtered_abundance_file}")
    
    # Save a general filtered_abundance.tsv file if we're processing bacteria (or all)
    if args.kingdom == "bacteria" or args.kingdom == "all":
        general_filtered_file = os.path.join(output_dir, "filtered_abundance.tsv")
        filtered_df.to_csv(general_filtered_file, sep='\t')
        logger.info(f"Filtered abundance data saved to {general_filtered_file} (for compatibility)")
    
    logger.info("Sylph data processing completed successfully")

if __name__ == "__main__":
    main()