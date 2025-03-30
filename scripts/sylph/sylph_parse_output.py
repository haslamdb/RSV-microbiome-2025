#!/usr/bin/env python3
"""
Parse and process Sylph output files for microbiome analysis.

This script:
1. Parses Sylph profiling output files (bacterial, viral, fungal)
2. Extracts taxonomic information and abundance data
3. Combines results from multiple samples into unified abundance tables
4. Creates separate tables for bacteria, viruses, and fungi
5. Optionally joins with metadata for downstream analysis

Usage:
    python process_sylph_outputs.py [--config CONFIG_FILE]
"""

import os
import sys
import glob
import argparse
import pandas as pd
import numpy as np
import yaml
from pathlib import Path
import re

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Process Sylph output files')
    parser.add_argument('--input-dir', type=str, default='data/SylphProfiles',
                       help='Directory containing Sylph output files')
    parser.add_argument('--output-dir', type=str, default='data/processed',
                       help='Directory to save processed files')
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
        # Debug: Print first few lines of file
        print(f"\nDebug - First 5 lines of {os.path.basename(file_path)}:")
        with open(file_path, 'r') as f:
            for i, line in enumerate(f):
                if i < 5:
                    print(f"  Line {i+1}: {line.strip()}")
                else:
                    break
        
        # The file doesn't have a clean header - need to handle specially
        with open(file_path, 'r') as f:
            header_line = f.readline().strip()
        
        # Clean up header line and split into column names
        header_columns = re.split(r'\s+', header_line.strip())
        print(f"  Detected columns: {header_columns}")
        
        # Read the file with the cleaned column names
        try:
            df = pd.read_csv(file_path, sep='\t', skiprows=1, names=header_columns)
            print(f"  Successfully read file with {len(df)} rows")
            
            # Print column names to debug
            print(f"  Actual columns in DataFrame: {df.columns.tolist()}")
            
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
            match = re.search(r'(?:profiled_|)(\d+-DNA)_', os.path.basename(f))
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
        # Instead of raising an error, return an empty DataFrame with a warning
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
        metadata_df.set_index(sample_id_column, inplace=True)
        return metadata_df
    except Exception as e:
        print(f"Error loading metadata file: {str(e)}")
        return pd.DataFrame()

def process_sylph_profile_type(input_dir, output_dir, profile_type, pattern, args):
    """
    Process one type of Sylph profile (bacterial, viral, or fungal).
    
    Parameters:
    -----------
    input_dir : str
        Directory containing Sylph output files
    output_dir : str
        Directory to save processed files
    profile_type : str
        Type of profile ('bacteria', 'viruses', or 'fungi')
    pattern : str
        Glob pattern to match files
    args : argparse.Namespace
        Command line arguments
        
    Returns:
    --------
    pandas.DataFrame or None
        Combined abundance table for this profile type
    """
    # Find all matching files
    files = glob.glob(os.path.join(input_dir, pattern))
    
    if not files:
        print(f"No {profile_type} profile files found matching pattern: {pattern}")
        return None
    
    print(f"\nProcessing {len(files)} {profile_type} profile files...")
    
    # Sample IDs will be determined inside the combine function
    combined_df = combine_sylph_samples(
        files, 
        abundance_type=args.abundance_type,
        min_ani=args.min_ani,
        min_coverage=args.min_coverage
    )
    
    if combined_df.empty:
        print(f"No usable data found in {profile_type} profile files.")
        return None
    
    # Save combined table
    output_file = os.path.join(output_dir, f'sylph_{profile_type}_abundance.csv')
    combined_df.to_csv(output_file)
    print(f"Combined {profile_type} abundance table saved to {output_file}")
    print(f"Table contains {len(combined_df.index)} taxa across {len(combined_df.columns)} samples")
    
    return combined_df


def main():
    """Main function to process Sylph output files."""
    # Parse arguments
    args = parse_args()
    
    # Set up directories
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Process each type of profile with the correct file naming pattern
    bacteria_df = process_sylph_profile_type(
        args.input_dir,
        args.output_dir,
        'bacteria',
        '*_profiled_bacteria.tsv',  # Fixed pattern to match actual files
        args
    )
    
    viruses_df = process_sylph_profile_type(
        args.input_dir,
        args.output_dir,
        'viruses',
        '*_profiled_viruses.tsv',  # Fixed pattern to match actual files
        args
    )
    
    fungi_df = process_sylph_profile_type(
        args.input_dir,
        args.output_dir,
        'fungi',
        '*_profiled_fungi.tsv',  # Fixed pattern to match actual files
        args
    )    
    # Load metadata if available
    if os.path.exists(args.metadata):
        print(f"\nLoading metadata from {args.metadata}")
        metadata_df = load_metadata(args.metadata, args.sample_id_column)
        
        if not metadata_df.empty:
            print("Checking sample overlap with metadata...")
            
            # Check overlap with each abundance table
            for name, df in [
                ('bacteria', bacteria_df),
                ('viruses', viruses_df),
                ('fungi', fungi_df)
            ]:
                if df is not None and not df.empty:
                    abundance_samples = set(df.columns)
                    metadata_samples = set(metadata_df.index)
                    common_samples = abundance_samples.intersection(metadata_samples)
                    
                    print(f"{name.capitalize()}: Found {len(common_samples)} samples with both abundance data and metadata")
                    
                    # Filter to common samples and save
                    if common_samples:
                        filtered_df = df[list(common_samples)]
                        output_file = os.path.join(args.output_dir, f'sylph_{name}_with_metadata_samples.csv')
                        filtered_df.to_csv(output_file)
                        print(f"Filtered {name} abundance table saved to {output_file}")
    
    # Create a combined table with all microbial types
    print("\nCreating combined microbial abundance table...")
    combined_tables = []
    
    if bacteria_df is not None and not bacteria_df.empty:
        bacteria_df.index = ['Bacteria: ' + idx for idx in bacteria_df.index]
        combined_tables.append(bacteria_df)
        
    if viruses_df is not None and not viruses_df.empty:
        viruses_df.index = ['Virus: ' + idx for idx in viruses_df.index]
        combined_tables.append(viruses_df)
        
    if fungi_df is not None and not fungi_df.empty:
        fungi_df.index = ['Fungus: ' + idx for idx in fungi_df.index]
        combined_tables.append(fungi_df)
    
    if combined_tables:
        all_microbes_df = pd.concat(combined_tables)
        output_file = os.path.join(args.output_dir, 'sylph_all_microbes_abundance.csv')
        all_microbes_df.to_csv(output_file)
        print(f"Combined microbial abundance table saved to {output_file}")
        print(f"Table contains {len(all_microbes_df.index)} taxa across {len(all_microbes_df.columns)} samples")
    else:
        print("No data to combine into a microbial abundance table.")
    
    print("\nProcessing complete!")

if __name__ == "__main__":
    main()
