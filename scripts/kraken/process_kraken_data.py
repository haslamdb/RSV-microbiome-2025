#!/usr/bin/env python
# scripts/kraken/process_kraken_data.py

import argparse
import os
import sys
import yaml
import pandas as pd
import numpy as np
import glob
import logging
from pathlib import Path

# Add kraken_tools to path
kraken_tools_path = os.path.expanduser("~/Documents/Code/kraken_tools")
if os.path.exists(kraken_tools_path):
    sys.path.append(kraken_tools_path)
else:
    print(f"Error: kraken_tools directory not found at {kraken_tools_path}")
    sys.exit(1)

try:
    from kraken_tools.analysis.data_processing import standardize_species_name
    from kraken_tools.analysis.abundance import normalize_abundance
    from kraken_tools.logger import setup_logger, log_print
except ImportError as e:
    print(f"Error importing from kraken_tools: {e}")
    print("Make sure kraken_tools is installed and accessible")
    sys.exit(1)

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process Kraken2/Bracken output files")
    
    parser.add_argument(
        "--config", 
        help="Path to configuration file (YAML)"
    )
    
    parser.add_argument(
        "--kreport-dir", 
        help="Directory containing Kraken2 report files"
    )
    
    parser.add_argument(
        "--bracken-dir", 
        help="Directory containing Bracken abundance files"
    )
    
    parser.add_argument(
        "--output-dir",
        default="results/kraken_analysis",
        help="Directory for output files (default: results/kraken_analysis)"
    )
    
    parser.add_argument(
        "--taxonomic-level", 
        default="S", 
        choices=["D", "P", "C", "O", "F", "G", "S"],
        help="Taxonomic level for analysis - D=domain, P=phylum, C=class, O=order, F=family, G=genus, S=species (default: S)"
    )
    
    parser.add_argument(
        "--min-abundance", 
        type=float, 
        default=0.01, 
        help="Minimum relative abundance threshold (default: 0.01)"
    )
    
    parser.add_argument(
        "--min-prevalence", 
        type=float, 
        default=0.1, 
        help="Minimum prevalence threshold (default: 0.1)"
    )
    
    parser.add_argument(
        "--normalize", 
        action="store_true",
        default=True,
        help="Enable data normalization (default: True)"
    )
    
    parser.add_argument(
        "--normalization-method", 
        default="clr", 
        choices=["relabundance", "cpm", "log10", "clr"],
        help="Normalization method (default: clr)"
    )
    
    parser.add_argument(
        "--metadata",
        help="Path to metadata file"
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

def load_kraken_file(filepath, logger):
    """Load a single Kraken report file."""
    try:
        # Determine columns based on file type
        if filepath.endswith('.kreport'):
            # Kraken2 report format
            df = pd.read_csv(filepath, sep='\t', header=None)
            if len(df.columns) >= 6:
                df.columns = ['percentage', 'clade_reads', 'taxon_reads', 'rank', 'taxid', 'name']
                # Filter by taxonomic rank (S for species)
                df = df[df['rank'] == 'S'].copy()
                # Extract sample ID from filename
                sample_id = os.path.basename(filepath).split('.')[0]
                # Standardize names and create result dataframe
                df['name'] = df['name'].apply(lambda x: standardize_species_name(x.strip()))
                result_df = df[['name', 'taxon_reads']].copy()
                result_df.columns = ['Species', sample_id]
                return result_df
            else:
                logger.warning(f"Not enough columns in file {filepath}")
                return None
                
        elif filepath.endswith('.bracken'):
            # Bracken format
            df = pd.read_csv(filepath, sep='\t')
            # Extract sample ID from filename
            sample_id = os.path.basename(filepath).split('.')[0]
            # Standardize species names
            df['name'] = df['name'].apply(standardize_species_name)
            # Create result dataframe
            result_df = df[['name', 'new_est_reads']].copy()
            result_df.columns = ['Species', sample_id]
            return result_df
            
        else:
            # Try generic format
            df = pd.read_csv(filepath, sep='\t')
            if 'name' in df.columns and len(df.columns) >= 6:
                # Extract sample ID from filename
                sample_id = os.path.basename(filepath).split('.')[0]
                # Standardize species names
                df['name'] = df['name'].apply(standardize_species_name)
                # Create result dataframe with the 6th column (assumed to be read counts)
                read_col = df.columns[5]
                result_df = df[['name', read_col]].copy()
                result_df.columns = ['Species', sample_id]
                return result_df
            else:
                logger.warning(f"Unrecognized format in file {filepath}")
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
    
    # Merge with the rest based on Species
    for df in file_dfs[1:]:
        if df is not None:
            merged_df = pd.merge(
                merged_df, 
                df, 
                on='Species', 
                how='outer'
            )
    
    # Fill NaN values with 0
    numeric_cols = [col for col in merged_df.columns if col != 'Species']
    merged_df[numeric_cols] = merged_df[numeric_cols].fillna(0)
    
    return merged_df

def filter_abundance_data(abundance_df, min_abundance, min_prevalence, logger):
    """Filter abundance data based on minimum abundance and prevalence."""
    logger.info(f"Filtering with min_abundance={min_abundance}, min_prevalence={min_prevalence}")
    
    # Get only the abundance columns (excluding Species)
    abundance_cols = [col for col in abundance_df.columns if col != 'Species']
    
    # Calculate total reads per sample for relative abundance
    sample_totals = abundance_df[abundance_cols].sum()
    
    # Create a copy with relative abundances
    rel_abundance_df = abundance_df.copy()
    for col in abundance_cols:
        if sample_totals[col] > 0:
            rel_abundance_df[col] = abundance_df[col] / sample_totals[col]
        else:
            rel_abundance_df[col] = 0
    
    # Calculate metrics for filtering
    # Maximum relative abundance across samples
    max_abundance = rel_abundance_df[abundance_cols].max(axis=1)
    
    # Prevalence (proportion of samples where taxon is present)
    prevalence = (abundance_df[abundance_cols] > 0).mean(axis=1)
    
    # Apply filters
    keep_taxa = (max_abundance >= min_abundance) & (prevalence >= min_prevalence)
    filtered_df = abundance_df.loc[keep_taxa].copy()
    
    logger.info(f"Filtered from {len(abundance_df)} to {len(filtered_df)} taxa")
    
    return filtered_df

def merge_with_metadata(abundance_df, metadata_file, logger):
    """Merge abundance data with sample metadata."""
    try:
        # Read metadata
        metadata_df = pd.read_csv(metadata_file)
        logger.info(f"Read metadata with {len(metadata_df)} samples")
        
        # Transpose abundance data for merging (samples as rows)
        abundance_t = abundance_df.set_index('Species').T.reset_index()
        abundance_t.rename(columns={'index': 'SampleID'}, inplace=True)
        
        # Merge with metadata
        merged_df = pd.merge(metadata_df, abundance_t, on='SampleID', how='inner')
        logger.info(f"Merged data with {len(merged_df)} samples")
        
        return merged_df
    except Exception as e:
        logger.error(f"Error merging with metadata: {str(e)}")
        return None

def perform_normalization(abundance_df, method, output_file, logger):
    """Apply normalization to the abundance data"""
    # First, create a copy of the dataframe with Species as index
    abundance_for_norm = abundance_df.copy()
    abundance_for_norm.set_index('Species', inplace=True)
    
    # Convert to format expected by normalize_abundance function
    temp_file = output_file.replace('.tsv', '_temp.tsv')
    abundance_for_norm.to_csv(temp_file, sep='\t')
    
    # Apply normalization
    try:
        normalized_file = normalize_abundance(
            abundance_file=temp_file,
            output_file=output_file,
            method=method,
            logger=logger
        )
        
        # Read back the normalized data and add the Species column back as a regular column
        normalized_df = pd.read_csv(normalized_file, sep='\t', index_col=0)
        normalized_df.reset_index(inplace=True)
        
        # Clean up temporary file
        if os.path.exists(temp_file):
            os.remove(temp_file)
        
        return normalized_df
    except Exception as e:
        logger.error(f"Error during normalization: {str(e)}")
        # Clean up temporary file
        if os.path.exists(temp_file):
            os.remove(temp_file)
        return abundance_df  # Return original on error

def main():
    # Parse command line arguments
    args = parse_arguments()
    
    # Setup logging
    logger = setup_logger(log_file=args.log_file, log_level=getattr(logging, args.log_level))
    log_print("Starting Kraken/Bracken data processing", level="info")
    
    # Load configuration if provided
    config = None
    if args.config:
        try:
            config = load_config(args.config)
            log_print(f"Loaded configuration from {args.config}", level="info")
            
            # Override command line arguments with config values if not explicitly provided
            if not args.min_abundance and 'differential_abundance' in config and 'min_abundance' in config['differential_abundance']:
                args.min_abundance = config['differential_abundance']['min_abundance']
            
            if not args.min_prevalence and 'differential_abundance' in config and 'min_prevalence' in config['differential_abundance']:
                args.min_prevalence = config['differential_abundance']['min_prevalence']
            
            if not args.metadata and 'metadata' in config and 'filename' in config['metadata']:
                args.metadata = os.path.join(os.path.dirname(args.config), '..', config['metadata']['filename'])
        except Exception as e:
            log_print(f"Error loading configuration: {str(e)}", level="error")
            sys.exit(1)
    
    # Check that at least one of kreport_dir or bracken_dir is provided
    if not args.kreport_dir and not args.bracken_dir:
        log_print("Error: At least one of --kreport-dir or --bracken-dir must be provided", level="error")
        sys.exit(1)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Process Kraken2 reports if provided
    kraken_dfs = []
    if args.kreport_dir:
        log_print(f"Processing Kraken2 reports from {args.kreport_dir}", level="info")
        kraken_files = glob.glob(os.path.join(args.kreport_dir, f"*.{args.taxonomic_level}.kreport"))
        
        if not kraken_files:
            # Try other common extensions
            kraken_files = glob.glob(os.path.join(args.kreport_dir, "*.kreport"))
            kraken_files.extend(glob.glob(os.path.join(args.kreport_dir, "*.kreport.txt")))
        
        log_print(f"Found {len(kraken_files)} Kraken2 report files", level="info")
        
        for kraken_file in kraken_files:
            df = load_kraken_file(kraken_file, logger)
            if df is not None:
                kraken_dfs.append(df)
        
        log_print(f"Successfully processed {len(kraken_dfs)} Kraken2 report files", level="info")
    
    # Process Bracken files if provided
    bracken_dfs = []
    if args.bracken_dir:
        log_print(f"Processing Bracken files from {args.bracken_dir}", level="info")
        bracken_files = glob.glob(os.path.join(args.bracken_dir, f"*.{args.taxonomic_level}.bracken"))
        
        if not bracken_files:
            # Try other common extensions
            bracken_files = glob.glob(os.path.join(args.bracken_dir, "*.bracken"))
            bracken_files.extend(glob.glob(os.path.join(args.bracken_dir, "*.bracken.txt")))
        
        log_print(f"Found {len(bracken_files)} Bracken files", level="info")
        
        for bracken_file in bracken_files:
            df = load_kraken_file(bracken_file, logger)
            if df is not None:
                bracken_dfs.append(df)
        
        log_print(f"Successfully processed {len(bracken_dfs)} Bracken files", level="info")
    
    # Combine all dataframes
    all_dfs = kraken_dfs + bracken_dfs
    
    if not all_dfs:
        log_print("Error: No valid Kraken2 or Bracken files were processed", level="error")
        sys.exit(1)
    
    # Merge abundance data
    log_print("Merging abundance data", level="info")
    merged_df = merge_abundance_data(all_dfs, logger)
    
    if merged_df is None:
        log_print("Error: Failed to merge abundance data", level="error")
        sys.exit(1)
    
    # Save raw abundance data
    raw_abundance_file = os.path.join(args.output_dir, "raw_abundance.tsv")
    merged_df.to_csv(raw_abundance_file, sep='\t', index=False)
    log_print(f"Raw abundance data saved to {raw_abundance_file}", level="info")
    
    # Filter abundance data
    log_print("Filtering abundance data", level="info")
    filtered_df = filter_abundance_data(merged_df, args.min_abundance, args.min_prevalence, logger)
    
    # Save filtered abundance data
    filtered_abundance_file = os.path.join(args.output_dir, f"filtered_{args.taxonomic_level}_abundance.tsv")
    filtered_df.to_csv(filtered_abundance_file, sep='\t', index=False)
    log_print(f"Filtered abundance data saved to {filtered_abundance_file}", level="info")
    
    # Apply normalization if requested
    normalized_df = filtered_df
    if args.normalize:
        log_print(f"Applying {args.normalization_method} normalization", level="info")
        normalized_file = os.path.join(args.output_dir, f"normalized_{args.normalization_method}_{args.taxonomic_level}_abundance.tsv")
        normalized_df = perform_normalization(filtered_df, args.normalization_method, normalized_file, logger)
        log_print(f"Normalized abundance data saved to {normalized_file}", level="info")
        
        # For compatibility, also save as normalized_abundance.tsv
        general_normalized_file = os.path.join(args.output_dir, "normalized_abundance.tsv")
        normalized_df.to_csv(general_normalized_file, sep='\t', index=False)
        log_print(f"Normalized abundance data saved to {general_normalized_file} (for compatibility)", level="info")
    else:
        # For compatibility, still save as filtered_abundance.tsv
        general_filtered_file = os.path.join(args.output_dir, "filtered_abundance.tsv")
        filtered_df.to_csv(general_filtered_file, sep='\t', index=False)
        log_print(f"Filtered abundance data saved to {general_filtered_file} (for compatibility)", level="info")
    
    # Merge with metadata if provided
    if args.metadata:
        log_print(f"Merging with metadata from {args.metadata}", level="info")
        # Use normalized data if available, otherwise use filtered data
        data_to_merge = normalized_df if args.normalize else filtered_df
        merged_with_metadata = merge_with_metadata(data_to_merge, args.metadata, logger)
        
        if merged_with_metadata is not None:
            # Save merged data
            merged_file = os.path.join(args.output_dir, "abundance_with_metadata.tsv")
            merged_with_metadata.to_csv(merged_file, sep='\t', index=False)
            log_print(f"Data with metadata saved to {merged_file}", level="info")
    
    log_print("Kraken/Bracken data processing completed successfully", level="info")

if __name__ == "__main__":
    main()