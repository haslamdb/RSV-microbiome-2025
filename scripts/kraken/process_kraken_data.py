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
from typing import Dict, List, Optional, Union, Set

# Add kraken_tools to path
kraken_tools_path = os.path.expanduser("~/Documents/Code/kraken_tools")
if os.path.exists(kraken_tools_path):
    sys.path.append(kraken_tools_path)
else:
    print(f"Error: kraken_tools directory not found at {kraken_tools_path}")
    sys.exit(1)

try:
    from kraken_tools.analysis.data_processing import standardize_species_name
    from kraken_tools.analysis.abundance import (
        read_kraken_report,
        read_bracken_abundance_file,
        normalize_abundance,
        merge_kraken_reports,
        merge_bracken_files
    )
    from kraken_tools.logger import setup_logger, log_print
    from kraken_tools.utils.metadata_utils import (
        collect_samples_from_metadata,
        read_samples_file
    )
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
    
    # Added arguments for sample selection
    parser.add_argument(
        "--sample-list",
        help="Path to file containing list of sample IDs (one per line)"
    )
    
    parser.add_argument(
        "--sample-col",
        help="Column name for sample identifiers in metadata file"
    )
    
    parser.add_argument(
        "--select-by-metadata",
        help="Filter samples by metadata criteria (e.g., 'StudyGroup==RSV' or 'Age>2')"
    )
    
    parser.add_argument(
        "--r1-col",
        help="Column name for R1 file paths in metadata file"
    )
    
    parser.add_argument(
        "--r2-col",
        help="Column name for R2 file paths in metadata file"
    )
    
    parser.add_argument(
        "--file-pattern",
        help="Pattern for finding sequence files (e.g., '{sample}_S*_R*.fastq.gz')"
    )
    
    parser.add_argument(
        "--r1-suffix",
        help="Suffix for R1 files (e.g., '_R1.fastq.gz')"
    )
    
    parser.add_argument(
        "--r2-suffix",
        help="Suffix for R2 files (e.g., '_R2.fastq.gz')"
    )
    
    parser.add_argument(
        "--paired",
        action="store_true",
        default=False,
        help="Whether samples are paired-end reads (default: False)"
    )
    
    parser.add_argument(
        "--seq-dir",
        help="Directory containing sequence files"
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
        "--use-direct-merge",
        action="store_true",
        default=False,
        help="Use the direct merge functionality from kraken_tools (default: False)"
    )
    
    return parser.parse_args()

def load_config(config_path):
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
    return config

def load_sample_list(sample_list_file: str, logger) -> Set[str]:
    """
    Load a list of sample IDs from a file.
    
    Args:
        sample_list_file: Path to a file containing one sample ID per line
        logger: Logger instance
        
    Returns:
        Set of sample IDs
    """
    try:
        with open(sample_list_file, 'r') as f:
            samples = {line.strip() for line in f if line.strip() and not line.startswith('#')}
        logger.info(f"Loaded {len(samples)} sample IDs from {sample_list_file}")
        return samples
    except Exception as e:
        logger.error(f"Error loading sample list from {sample_list_file}: {str(e)}")
        return set()

def filter_samples_by_metadata_criteria(metadata_df: pd.DataFrame, criteria: str, logger) -> Set[str]:
    """
    Filter sample IDs based on metadata criteria.
    
    Args:
        metadata_df: DataFrame containing sample metadata
        criteria: String expression for filtering (e.g., 'StudyGroup=="RSV"' or 'Age>2')
        logger: Logger instance
        
    Returns:
        Set of sample IDs that match the criteria
    """
    try:
        # Try to determine the sample ID column
        sample_cols = [col for col in metadata_df.columns if 'sample' in col.lower() or 'id' in col.lower()]
        if not sample_cols:
            sample_col = metadata_df.columns[0]
            logger.warning(f"Could not identify sample ID column, using first column: {sample_col}")
        else:
            sample_col = sample_cols[0]
            logger.info(f"Using {sample_col} as sample ID column")
        
        # Apply query to filter metadata
        filtered_df = metadata_df.query(criteria)
        samples = set(filtered_df[sample_col].astype(str).values)
        
        logger.info(f"Filtered {len(samples)} samples matching criteria: {criteria}")
        return samples
    except Exception as e:
        logger.error(f"Error filtering samples by criteria '{criteria}': {str(e)}")
        return set()

def get_samples_from_metadata(metadata_file: str, criteria: str = None, sample_col: str = None, logger=None) -> Set[str]:
    """
    Get sample IDs from a metadata file, optionally filtered by criteria.
    
    Args:
        metadata_file: Path to metadata CSV file
        criteria: Optional criteria for filtering samples
        sample_col: Column name for sample identifiers
        logger: Logger instance
        
    Returns:
        Set of sample IDs
    """
    try:
        # Read metadata file
        metadata_df = pd.read_csv(metadata_file)
        logger.info(f"Read metadata file with {len(metadata_df)} rows and columns: {metadata_df.columns.tolist()}")
        
        # Identify sample ID column if not specified
        if not sample_col:
            # Try to auto-detect
            common_sample_cols = ['SampleID', 'Sample_ID', 'SampleName', 'sample_id', 'sample_name', 'Sample', 'ID']
            for col in common_sample_cols:
                if col in metadata_df.columns:
                    sample_col = col
                    logger.info(f"Auto-detected sample ID column: {sample_col}")
                    break
            
            if not sample_col:
                sample_col = metadata_df.columns[0]
                logger.warning(f"Could not auto-detect sample ID column, using the first column: {sample_col}")
        
        # Apply criteria if provided
        if criteria:
            samples = filter_samples_by_metadata_criteria(metadata_df, criteria, logger)
        else:
            # If no criteria, include all samples from the metadata
            samples = set(metadata_df[sample_col].astype(str).values)
            logger.info(f"Including all {len(samples)} samples from metadata")
        
        # Filter out empty or NaN values
        samples = {s for s in samples if s and s != 'nan'}
        
        return samples
    except Exception as e:
        logger.error(f"Error processing metadata file {metadata_file}: {str(e)}")
        return set()

def find_matching_files(sample_ids: Set[str], search_dir: str, file_pattern: str = None,
                       r1_suffix: str = None, r2_suffix: str = None, paired: bool = False,
                       logger=None) -> Dict[str, List[str]]:
    """
    Find files matching the given sample IDs based on patterns or suffixes.
    
    Args:
        sample_ids: Set of sample IDs to search for
        search_dir: Directory to search in
        file_pattern: Pattern to use for file search
        r1_suffix: Suffix for R1 files
        r2_suffix: Suffix for R2 files
        paired: Whether to look for paired-end files
        logger: Logger instance
        
    Returns:
        Dictionary mapping sample IDs to file paths
    """
    samples_dict = {}
    
    from kraken_tools.utils.metadata_utils import find_sample_files
    
    for sample_id in sample_ids:
        # Use find_sample_files directly from kraken_tools
        files = find_sample_files(
            sample_id=sample_id,
            search_dir=search_dir,
            file_pattern=file_pattern,
            r1_suffix=r1_suffix,
            r2_suffix=r2_suffix,
            paired=paired
        )
        
        # If files were found, add to our dictionary
        if files:
            samples_dict[sample_id] = files
    
    logger.info(f"Found files for {len(samples_dict)} out of {len(sample_ids)} samples")
    return samples_dict

def filter_files_by_samples(files_dict: Dict[str, str], sample_ids: Set[str], logger) -> Dict[str, str]:
    """
    Filter a dictionary of files to only include specified sample IDs.
    
    Args:
        files_dict: Dictionary mapping sample IDs to file paths
        sample_ids: Set of sample IDs to include
        logger: Logger instance
        
    Returns:
        Filtered dictionary
    """
    filtered_dict = {sample_id: files_dict[sample_id] for sample_id in sample_ids if sample_id in files_dict}
    logger.info(f"Filtered to {len(filtered_dict)} files for selected samples (from {len(files_dict)} total files)")
    return filtered_dict

def load_kraken_file(filepath, logger, taxonomic_level='S'):
    """
    Load a single Kraken report or Bracken file.
    
    Args:
        filepath: Path to the Kraken or Bracken file
        logger: Logger instance
        taxonomic_level: Taxonomic level to extract (default: 'S' for species)
        
    Returns:
        DataFrame with species and abundance information for the sample
    """
    try:
        # Get filename and extract sample_id
        file_base = os.path.basename(filepath)
        if "_abundance.txt" in file_base:
            # Format like: 3000829042-DNA_S_abundance.txt
            sample_id = file_base.split('_')[0]
        else:
            # Traditional format
            sample_id = file_base.split('.')[0]
        
        if filepath.endswith('.kreport'):
            # Use the new read_kraken_report function from kraken_tools
            df = read_kraken_report(filepath, taxonomic_level, logger)
            
            if df is not None:
                # Standardize names
                df['taxon_name'] = df['taxon_name'].apply(standardize_species_name)
                
                # Create result dataframe with the taxon reads (use actual read counts)
                result_df = df[['taxon_name', 'taxon_reads']].copy()
                result_df.columns = ['Species', sample_id]
                return result_df
            else:
                logger.warning(f"Unable to read Kraken report file {filepath}")
                return None
                
        elif filepath.endswith('.bracken'):
            # Use the new read_bracken_abundance_file function from kraken_tools
            df = read_bracken_abundance_file(filepath, logger)
            
            if df is not None:
                # Standardize names
                df['taxon_name'] = df['taxon_name'].apply(standardize_species_name)
                
                # Create result dataframe with the breads (use breads count instead of percentage)
                result_df = df[['taxon_name', 'breads']].copy()
                result_df.columns = ['Species', sample_id]
                return result_df
            else:
                logger.warning(f"Unable to read Bracken file {filepath}")
                return None
        
        elif filepath.endswith('_abundance.txt'):
            # Handle custom format like 3000829042-DNA_S_abundance.txt
            try:
                df = pd.read_csv(filepath, sep='\t')
                if 'name' in df.columns and 'new_est_reads' in df.columns:
                    # Standardize species names
                    df['name'] = df['name'].apply(standardize_species_name)
                    # Use new_est_reads column for read counts
                    result_df = df[['name', 'new_est_reads']].copy()
                    result_df.columns = ['Species', sample_id]
                    return result_df
                else:
                    logger.warning(f"Missing expected columns in abundance file {filepath}")
                    return None
            except Exception as e:
                logger.error(f"Error reading abundance file {filepath}: {str(e)}")
                return None
            
        else:
            # Try generic format
            df = pd.read_csv(filepath, sep='\t')
            if 'name' in df.columns and len(df.columns) >= 6:
                # Extract sample ID from filename
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

def filter_human_reads(abundance_df, logger, taxonomic_level='S'):
    """
    Filter out human genome reads from abundance data.
    
    Args:
        abundance_df: DataFrame with abundance data (raw read counts)
        logger: Logger instance
        taxonomic_level: Taxonomic level of the data (default: 'S' for species)
        
    Returns:
        Filtered abundance DataFrame with human reads removed
    """
    logger.info("Filtering out human genome reads")
    
    # Get only the abundance columns (excluding Species)
    abundance_cols = [col for col in abundance_df.columns if col != 'Species']
    
    # Create copy of the dataframe
    filtered_df = abundance_df.copy()
    
    # Set human patterns based on taxonomic level
    if taxonomic_level == 'S':
        human_patterns = ['Homo sapiens']
        logger.info("Using species-level filter: 'Homo sapiens'")
    elif taxonomic_level == 'G':
        human_patterns = ['Homo']
        logger.info("Using genus-level filter: 'Homo'")
    else:
        human_patterns = ['Homo sapiens', 'Homo']
        logger.info("Using default filters: 'Homo sapiens', 'Homo'")
    
    human_taxa = []
    
    for pattern in human_patterns:
        matches = filtered_df[filtered_df['Species'].str.contains(pattern, case=False)]
        if not matches.empty:
            human_taxa.extend(matches['Species'].tolist())
    
    if human_taxa:
        # Drop human taxa rows
        filtered_df = filtered_df[~filtered_df['Species'].isin(human_taxa)]
        logger.info(f"Removed {len(human_taxa)} human taxa: {', '.join(human_taxa)}")
    else:
        logger.info("No human taxa found to remove")
    
    # Calculate read counts before and after removal
    before_total = abundance_df[abundance_cols].sum().sum()
    after_total = filtered_df[abundance_cols].sum().sum()
    removed_percentage = 100 * (before_total - after_total) / before_total if before_total > 0 else 0
    
    logger.info(f"Removed {before_total - after_total:,} human reads ({removed_percentage:.2f}% of total)")
    
    return filtered_df

def filter_abundance_data(abundance_df, min_abundance, min_prevalence, logger):
    """
    Filter abundance data based on minimum abundance and prevalence.
    Works with raw read count data by converting to relative abundance for filtering.
    
    Args:
        abundance_df: DataFrame with abundance data (raw read counts)
        min_abundance: Minimum relative abundance threshold
        min_prevalence: Minimum prevalence threshold
        logger: Logger instance
        
    Returns:
        Filtered abundance DataFrame (still raw read counts)
    """
    logger.info(f"Filtering with min_abundance={min_abundance}, min_prevalence={min_prevalence}")
    
    # Get only the abundance columns (excluding Species)
    abundance_cols = [col for col in abundance_df.columns if col != 'Species']
    
    # Calculate total reads per sample
    sample_totals = abundance_df[abundance_cols].sum()
    logger.info(f"Sample read count stats: min={sample_totals.min()}, max={sample_totals.max()}, mean={sample_totals.mean():.1f}")
    
    # Create a copy with relative abundances (for filtering purposes)
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
    """
    Apply normalization to the abundance data using the updated normalize_abundance function.
    Handles different normalization methods for read count data.
    
    Args:
        abundance_df: DataFrame with abundance data (assumes read counts, not percentages)
        method: Normalization method ('relabundance', 'cpm', 'log10', 'clr')
        output_file: Path to save normalized data
        logger: Logger instance
        
    Returns:
        Normalized DataFrame
    """
    # First, create a copy of the dataframe with Species as index
    abundance_for_norm = abundance_df.copy()
    abundance_for_norm.set_index('Species', inplace=True)
    
    # For CLR transformation, we need to make sure we have non-zero values
    # Add pseudocount of 1 to avoid issues with zeros in the data
    if method == 'clr':
        logger.info("Adding pseudocount of 1 to abundance data before CLR transformation")
        abundance_for_norm = abundance_for_norm + 1
    
    # Convert to format expected by normalize_abundance function
    temp_file = output_file.replace('.tsv', '_temp.tsv')
    abundance_for_norm.to_csv(temp_file, sep='\t')
    
    # Apply normalization using the updated function
    try:
        # Note: Our data is now raw read counts, not percentages
        # This requires adjustment in how we approach normalization
        
        if method == 'relabundance':
            logger.info("Converting read counts to relative abundance")
            # Calculate column sums for normalization
            col_sums = abundance_for_norm.sum(axis=0)
            # Normalize by column sum to get relative abundance (0-1)
            norm_df = abundance_for_norm.div(col_sums, axis=1)
            norm_df.to_csv(output_file, sep='\t')
            normalized_file = output_file
        elif method == 'cpm':
            logger.info("Normalizing to counts per million")
            # Calculate column sums
            col_sums = abundance_for_norm.sum(axis=0)
            # Normalize to counts per million
            norm_df = abundance_for_norm.div(col_sums, axis=1) * 1000000
            norm_df.to_csv(output_file, sep='\t')
            normalized_file = output_file
        else:
            # For log10 and clr transformations, use the library function
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
            
            # Get sample selection parameters from config if provided
            if 'sample_selection' in config:
                sample_config = config['sample_selection']
                if 'sample_list' in sample_config and not args.sample_list:
                    args.sample_list = os.path.join(os.path.dirname(args.config), '..', sample_config['sample_list'])
                if 'sample_col' in sample_config and not args.sample_col:
                    args.sample_col = sample_config['sample_col']
                if 'select_by_metadata' in sample_config and not args.select_by_metadata:
                    args.select_by_metadata = sample_config['select_by_metadata']
                if 'seq_dir' in sample_config and not args.seq_dir:
                    args.seq_dir = sample_config['seq_dir']
                if 'file_pattern' in sample_config and not args.file_pattern:
                    args.file_pattern = sample_config['file_pattern']
                if 'r1_suffix' in sample_config and not args.r1_suffix:
                    args.r1_suffix = sample_config['r1_suffix']
                if 'r2_suffix' in sample_config and not args.r2_suffix:
                    args.r2_suffix = sample_config['r2_suffix']
                if 'paired' in sample_config and not args.paired:
                    args.paired = sample_config['paired']
        except Exception as e:
            log_print(f"Error loading configuration: {str(e)}", level="error")
            sys.exit(1)
    
    # Check that at least one of kreport_dir or bracken_dir is provided
    if not args.kreport_dir and not args.bracken_dir:
        log_print("Error: At least one of --kreport-dir or --bracken-dir must be provided", level="error")
        sys.exit(1)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Determine which samples to include in the analysis
    selected_samples = set()
    
    # Method 1: Get samples from a sample list file
    if args.sample_list:
        log_print(f"Loading sample list from {args.sample_list}", level="info")
        sample_list_samples = load_sample_list(args.sample_list, logger)
        selected_samples.update(sample_list_samples)
        log_print(f"Added {len(sample_list_samples)} samples from sample list", level="info")
    
    # Method 2: Get samples from metadata (optionally filtered by criteria)
    if args.metadata:
        log_print(f"Loading samples from metadata {args.metadata}", level="info")
        metadata_samples = get_samples_from_metadata(
            metadata_file=args.metadata,
            criteria=args.select_by_metadata,
            sample_col=args.sample_col,
            logger=logger
        )
        # If no samples selected yet, use all from metadata, otherwise use the intersection
        if not selected_samples:
            selected_samples.update(metadata_samples)
            log_print(f"Added {len(metadata_samples)} samples from metadata", level="info")
        else:
            # Intersection - only keep samples that are in both sets
            initial_count = len(selected_samples)
            selected_samples = selected_samples.intersection(metadata_samples)
            log_print(f"Filtered selected samples from {initial_count} to {len(selected_samples)} using metadata", level="info")
    
    # Process Kraken2 reports if provided
    kraken_dfs = []
    kraken_files_dict = {}
    
    if args.kreport_dir:
        log_print(f"Processing Kraken2 reports from {args.kreport_dir}", level="info")
        kraken_files = glob.glob(os.path.join(args.kreport_dir, f"*.{args.taxonomic_level}.kreport"))
        
        if not kraken_files:
            # Try other common extensions
            kraken_files = glob.glob(os.path.join(args.kreport_dir, "*.kreport"))
            kraken_files.extend(glob.glob(os.path.join(args.kreport_dir, "*.kreport.txt")))
        
        log_print(f"Found {len(kraken_files)} Kraken2 report files", level="info")
        
        # Create a dictionary of sample IDs to file paths
        for kraken_file in kraken_files:
            sample_id = os.path.basename(kraken_file).split('.')[0]
            kraken_files_dict[sample_id] = kraken_file
        
        # Filter kraken files to only include selected samples
        if selected_samples:
            kraken_files_dict = filter_files_by_samples(kraken_files_dict, selected_samples, logger)
            log_print(f"Filtered to {len(kraken_files_dict)} Kraken files based on sample selection", level="info")
            
            # Update kraken_files list to use only selected files
            kraken_files = list(kraken_files_dict.values())
        
        # Option 1: Use our updated load_kraken_file function 
        for kraken_file in kraken_files:
            df = load_kraken_file(kraken_file, logger, args.taxonomic_level)
            if df is not None:
                kraken_dfs.append(df)
        
        # Option 2: We could alternatively use merge_kraken_reports from the kraken_tools
        # to directly merge all files, but we'll use Option 1 for compatibility
        # with the rest of the code
        
        log_print(f"Successfully processed {len(kraken_dfs)} Kraken2 report files", level="info")
    
    # Process Bracken files if provided
    bracken_dfs = []
    bracken_files_dict = {}
    
    if args.bracken_dir:
        log_print(f"Processing Bracken files from {args.bracken_dir}", level="info")
        bracken_files = glob.glob(os.path.join(args.bracken_dir, f"*_{args.taxonomic_level}.bracken"))
        
        if not bracken_files:
            # Try other common extensions
            bracken_files = glob.glob(os.path.join(args.bracken_dir, "*.bracken"))
            bracken_files.extend(glob.glob(os.path.join(args.bracken_dir, "*.bracken.txt")))
            # Check for custom format like *_S_abundance.txt
            bracken_files.extend(glob.glob(os.path.join(args.bracken_dir, f"*_{args.taxonomic_level}_abundance.txt")))
        
        log_print(f"Found {len(bracken_files)} Bracken files", level="info")
        
        # Create a dictionary of sample IDs to file paths
        for bracken_file in bracken_files:
            file_base = os.path.basename(bracken_file)
            # Handle different filename formats
            if "_abundance.txt" in file_base:
                # Format like: 3000829042-DNA_S_abundance.txt
                sample_id = file_base.split('_')[0]
            else:
                # Traditional format
                sample_id = file_base.split('.')[0]
            bracken_files_dict[sample_id] = bracken_file
        
        # Filter bracken files to only include selected samples
        if selected_samples:
            bracken_files_dict = filter_files_by_samples(bracken_files_dict, selected_samples, logger)
            log_print(f"Filtered to {len(bracken_files_dict)} Bracken files based on sample selection", level="info")
            
            # Update bracken_files list to use only selected files
            bracken_files = list(bracken_files_dict.values())
            
        # Option 1: Use our updated load_kraken_file function
        for bracken_file in bracken_files:
            df = load_kraken_file(bracken_file, logger, args.taxonomic_level)
            if df is not None:
                bracken_dfs.append(df)
                
        # Option 2: We could alternatively use merge_bracken_files from kraken_tools
        # to directly merge all files, but we'll use Option 1 for compatibility
        # with the rest of the code
        
        log_print(f"Successfully processed {len(bracken_dfs)} Bracken files", level="info")
    
    # Flag to determine which approach to use
    use_direct_merge = args.use_direct_merge
    
    if use_direct_merge and (kraken_files_dict or bracken_files_dict):
        # Use the new direct merge functionality if any dictionaries are populated
        log_print("Using direct merge functionality from kraken_tools", level="info")
        
        merged_file = None
        
        # Process Kraken files directly if available
        # Note: The default merge_kraken_reports uses percentage, but we want read counts
        # We'll need to implement our own direct merge for Kraken that uses read counts
        if kraken_files_dict:
            log_print("Merging Kraken reports directly with read counts", level="info")
            kraken_merged_file = os.path.join(args.output_dir, f"kraken_{args.taxonomic_level}_merged.tsv")
            
            # Create a dict to store abundance data by taxon
            abundance_data = {}
            
            # Process each sample
            for sample_id, file_path in kraken_files_dict.items():
                log_print(f"Processing Kraken report for sample {sample_id}", level="info")
                
                # Read the Kraken report
                kraken_df = read_kraken_report(file_path, args.taxonomic_level, logger)
                if kraken_df is None:
                    log_print(f"Skipping sample {sample_id} due to file reading error", level="warning")
                    continue
                
                # Extract relevant columns - using taxon_reads instead of percentage
                taxa = kraken_df['taxon_name'].apply(standardize_species_name).values
                reads = kraken_df['taxon_reads'].values
                
                # Add to abundance dict
                for taxon, read_count in zip(taxa, reads):
                    if taxon not in abundance_data:
                        abundance_data[taxon] = {}
                    abundance_data[taxon][sample_id] = read_count
            
            if not abundance_data:
                log_print("No valid data found in any Kraken report", level="error")
                use_direct_merge = False
                merged_file = None
            else:
                # Convert to DataFrame
                merged_df_temp = pd.DataFrame(abundance_data).T
                
                # Fill missing values with 0
                merged_df_temp = merged_df_temp.fillna(0)
                
                # Save to file
                try:
                    os.makedirs(os.path.dirname(kraken_merged_file), exist_ok=True)
                    merged_df_temp.to_csv(kraken_merged_file, sep='\t')
                    log_print(f"Merged Kraken report data saved to {kraken_merged_file}", level="info")
                    merged_file = kraken_merged_file
                except Exception as e:
                    log_print(f"Error saving merged Kraken data: {str(e)}", level="error")
                    use_direct_merge = False
                    merged_file = None
        
        # Process Bracken files directly if available
        # Similar to Kraken, we need to implement our own direct merge for Bracken to use read counts
        elif bracken_files_dict:
            log_print("Merging Bracken files directly with read counts", level="info")
            bracken_merged_file = os.path.join(args.output_dir, f"bracken_{args.taxonomic_level}_merged.tsv")
            
            # Create a dict to store abundance data by taxon
            abundance_data = {}
            
            # Process each sample
            for sample_id, file_path in bracken_files_dict.items():
                log_print(f"Processing Bracken file for sample {sample_id}", level="info")
                
                # Read the Bracken file
                bracken_df = read_bracken_abundance_file(file_path, logger)
                if bracken_df is None:
                    log_print(f"Skipping sample {sample_id} due to file reading error", level="warning")
                    continue
                
                # Extract relevant columns - using breads (read counts) instead of breads_percentage
                taxa = bracken_df['taxon_name'].apply(standardize_species_name).values
                reads = bracken_df['breads'].values
                
                # Add to abundance dict
                for taxon, read_count in zip(taxa, reads):
                    if taxon not in abundance_data:
                        abundance_data[taxon] = {}
                    abundance_data[taxon][sample_id] = read_count
            
            if not abundance_data:
                log_print("No valid data found in any Bracken file", level="error")
                use_direct_merge = False
                merged_file = None
            else:
                # Convert to DataFrame
                merged_df_temp = pd.DataFrame(abundance_data).T
                
                # Fill missing values with 0
                merged_df_temp = merged_df_temp.fillna(0)
                
                # Save to file
                try:
                    os.makedirs(os.path.dirname(bracken_merged_file), exist_ok=True)
                    merged_df_temp.to_csv(bracken_merged_file, sep='\t')
                    log_print(f"Merged Bracken data saved to {bracken_merged_file}", level="info")
                    merged_file = bracken_merged_file
                except Exception as e:
                    log_print(f"Error saving merged Bracken data: {str(e)}", level="error")
                    use_direct_merge = False
                    merged_file = None
        
        # If we successfully merged files directly
        if use_direct_merge and merged_file:
            # Read the merged file
            merged_df = pd.read_csv(merged_file, sep='\t', index_col=0)
            # Add 'Species' column by creating a new DataFrame with index as column
            merged_df = merged_df.reset_index().rename(columns={'index': 'Species'})
            
            # Save raw abundance data (for consistency with the traditional approach)
            raw_abundance_file = os.path.join(args.output_dir, "raw_abundance.tsv")
            merged_df.to_csv(raw_abundance_file, sep='\t', index=False)
            log_print(f"Raw abundance data saved to {raw_abundance_file}", level="info")
            
            # Filter out human reads
            log_print("Removing human genome reads", level="info")
            human_filtered_df = filter_human_reads(merged_df, logger)
            
            # Save human-filtered abundance data
            human_filtered_file = os.path.join(args.output_dir, "human_filtered_abundance.tsv")
            human_filtered_df.to_csv(human_filtered_file, sep='\t', index=False)
            log_print(f"Human-filtered abundance data saved to {human_filtered_file}", level="info")
            
            # Replace the merged_df with the human_filtered_df for further processing
            merged_df = human_filtered_df
        else:
            # If direct merge failed, fall back to the traditional approach
            use_direct_merge = False
    
    # If not using direct merge or it failed, use the traditional approach
    if not use_direct_merge:
        # Combine all dataframes
        all_dfs = kraken_dfs + bracken_dfs
        
        if not all_dfs:
            log_print("Error: No valid Kraken2 or Bracken files were processed", level="error")
            sys.exit(1)
        
        # Merge abundance data
        log_print("Merging abundance data using traditional approach", level="info")
        merged_df = merge_abundance_data(all_dfs, logger)
        
        if merged_df is None:
            log_print("Error: Failed to merge abundance data", level="error")
            sys.exit(1)
        
        # Save raw abundance data (before human filtering)
        raw_abundance_file = os.path.join(args.output_dir, "raw_abundance.tsv")
        merged_df.to_csv(raw_abundance_file, sep='\t', index=False)
        log_print(f"Raw abundance data saved to {raw_abundance_file}", level="info")
    
    # Filter out human reads
    log_print("Removing human genome reads", level="info")
    human_filtered_df = filter_human_reads(merged_df, logger)
    
    # Save human-filtered abundance data
    human_filtered_file = os.path.join(args.output_dir, "human_filtered_abundance.tsv")
    human_filtered_df.to_csv(human_filtered_file, sep='\t', index=False)
    log_print(f"Human-filtered abundance data saved to {human_filtered_file}", level="info")
    
    # Replace the merged_df with the human_filtered_df for further processing
    # This ensures the rest of the pipeline continues with human-filtered data
    merged_df = human_filtered_df
    
    # Filter abundance data based on min abundance and prevalence
    log_print("Filtering abundance data", level="info")
    filtered_df = filter_abundance_data(human_filtered_df, args.min_abundance, args.min_prevalence, logger)
    
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