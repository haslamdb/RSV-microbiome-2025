#!/usr/bin/env python
# scripts/01_process_metaphlan_files.py

"""
Process raw MetaPhlAn output files and generate a combined abundance table.

This script:
1. Finds all MetaPhlAn output files in the data/raw directory
2. Parses them to extract species-level abundances
3. Combines them into a unified abundance table
4. Saves the result to data/processed/combined_abundance.csv

Usage:
    python scripts/01_process_metaphlan_files.py

"""

import os
import glob
import pandas as pd
import yaml
import sys
from pathlib import Path

# Add project root to Python path
project_root = Path(__file__).resolve().parents[1]
sys.path.append(str(project_root))

# Add tools directory to Python path
tools_dir = project_root / 'tools'
sys.path.append(str(tools_dir))

# Now import metaphlan_tools functions
from metaphlan_tools import parse_metaphlan_file, combine_samples, load_metadata


# Import functions from metaphlan_tools
from metaphlan_tools import (
    load_metadata,
    combine_samples
)


def main():
    # Load configuration
    config_path = project_root / 'config' / 'analysis_parameters.yml'
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    # Set up paths
    raw_data_dir = project_root / 'data' / 'raw'
    processed_data_dir = project_root / 'data' / 'processed'
    processed_data_dir.mkdir(exist_ok=True, parents=True)
    
    # Find all MetaPhlAn output files
    file_pattern = os.path.join(raw_data_dir, config['metaphlan']['file_pattern'])
    files = glob.glob(file_pattern)
    
    if not files:
        print(f"No files found matching pattern: {file_pattern}")
        return
    
    print(f"Found {len(files)} MetaPhlAn output files.")
    
    # Extract sample IDs from filenames
    sample_ids = []
    for file_path in files:
        # Extract sample ID from filename based on pattern in config
        sample_id = os.path.basename(file_path).split(config['metaphlan']['sample_id_delimiter'])[0]
        sample_ids.append(sample_id)
    
    # Combine files into a single abundance table
    print("Combining files into abundance table...")
    abundance_df = combine_samples(files, sample_ids)
    
    # Save the combined table
    output_file = processed_data_dir / 'combined_abundance.csv'
    abundance_df.to_csv(output_file)
    print(f"Combined abundance table saved to {output_file}")
    print(f"Table contains {len(abundance_df.index)} species across {len(abundance_df.columns)} samples")
    
    # Optionally join with metadata
    if config['metadata']['join_with_abundance']:
        metadata_path = project_root / config['metadata']['filename']
        if metadata_path.exists():
            print("Joining with metadata...")
            # Convert Path object to string
            metadata_path_str = str(metadata_path)
            metadata_df = load_metadata(metadata_path_str, config['metadata']['sample_id_column'])
            
            
            # Check for sample overlap
            abundance_samples = set(abundance_df.columns)
            metadata_samples = set(metadata_df.index)
            common_samples = abundance_samples.intersection(metadata_samples)
            
            print(f"Found {len(common_samples)} samples with both abundance data and metadata")
            
            # Filter to common samples and save
            if common_samples:
                abundance_df_filtered = abundance_df[list(common_samples)]
                output_file_filtered = processed_data_dir / 'abundance_with_metadata_samples.csv'
                abundance_df_filtered.to_csv(output_file_filtered)
                print(f"Filtered abundance table saved to {output_file_filtered}")
        else:
            print(f"Metadata file not found: {metadata_path}")


    # Add diagnostic section to identify sample discrepancies
    print("\n=== Sample ID Diagnostics ===")
    metadata_path = project_root / config['metadata']['filename']
    if metadata_path.exists():
        # Load metadata
        metadata_path_str = str(metadata_path)
        metadata_df = load_metadata(metadata_path_str, config['metadata']['sample_id_column'])
        
        # Get sample IDs from both datasets
        abundance_samples = set(abundance_df.columns)
        metadata_samples = set(metadata_df.index)
        
        # Calculate differences
        only_in_abundance = abundance_samples - metadata_samples
        only_in_metadata = metadata_samples - abundance_samples
        common_samples = abundance_samples.intersection(metadata_samples)
        
        # Print summary
        print(f"Total samples in abundance data: {len(abundance_samples)}")
        print(f"Total samples in metadata: {len(metadata_samples)}")
        print(f"Samples in both datasets: {len(common_samples)}")
        print(f"Samples only in abundance data: {len(only_in_abundance)}")
        print(f"Samples only in metadata: {len(only_in_metadata)}")
        
        # Print some examples of mismatched samples
        if only_in_abundance:
            print("\nExample samples in abundance data but missing from metadata:")
            for sample in list(only_in_abundance)[:5]:  # Show first 5 examples
                print(f"  - {sample}")
            
        if only_in_metadata:
            print("\nExample samples in metadata but missing from abundance data:")
            for sample in list(only_in_metadata)[:5]:  # Show first 5 examples
                print(f"  - {sample}")
        
        # Check for case sensitivity issues
        if only_in_abundance and only_in_metadata:
            print("\nChecking for case sensitivity issues...")
            abundance_samples_lower = {s.lower() for s in abundance_samples}
            metadata_samples_lower = {s.lower() for s in metadata_df.index}
            
            case_insensitive_matches = 0
            for a_sample in only_in_abundance:
                if a_sample.lower() in metadata_samples_lower:
                    case_insensitive_matches += 1
                    
            if case_insensitive_matches > 0:
                print(f"Found {case_insensitive_matches} potential case sensitivity mismatches!")
                print("Consider standardizing sample IDs to the same case.")
        
        # Check for minor formatting differences (spaces, dashes, etc.)
        if only_in_abundance and only_in_metadata:
            print("\nChecking for formatting differences...")
            # Create normalized versions (remove special chars)
            import re
            
            def normalize_id(sample_id):
                return re.sub(r'[^a-zA-Z0-9]', '', str(sample_id))
            
            norm_abundance = {normalize_id(s): s for s in abundance_samples}
            norm_metadata = {normalize_id(s): s for s in metadata_samples}
            
            format_matches = 0
            for norm_a, orig_a in norm_abundance.items():
                if norm_a in norm_metadata and orig_a not in metadata_samples:
                    format_matches += 1
                    if format_matches <= 3:  # Show first 3 examples
                        print(f"  Abundance: '{orig_a}' --> Metadata: '{norm_metadata[norm_a]}'")
            
            if format_matches > 0:
                print(f"Found {format_matches} potential formatting differences in sample IDs!")
                print("Consider standardizing special characters and spacing in sample IDs.")
        
        # Write mismatch lists to files for further investigation
        if only_in_abundance:
            mismatch_file = processed_data_dir / 'samples_only_in_abundance.txt'
            with open(mismatch_file, 'w') as f:
                for sample in sorted(only_in_abundance):
                    f.write(f"{sample}\n")
            print(f"\nSamples only in abundance data saved to: {mismatch_file}")
        
        if only_in_metadata:
            # Save as TXT file
            mismatch_file_txt = processed_data_dir / 'samples_only_in_metadata.txt'
            with open(mismatch_file_txt, 'w') as f:
                for sample in sorted(only_in_metadata):
                    f.write(f"{sample}\n")
            print(f"Samples only in metadata saved to: {mismatch_file_txt}")
            
            # Also save as CSV file with more details
            mismatch_file_csv = processed_data_dir / 'samples_only_in_metadata.csv'
            # Get the metadata for these samples
            metadata_only_df = metadata_df.loc[list(only_in_metadata)]
            metadata_only_df.to_csv(mismatch_file_csv)
            print(f"Detailed metadata for missing samples saved to: {mismatch_file_csv}")
            
            # Create a subset of metadata file with only matched samples
            matched_metadata_file = processed_data_dir / 'metadata_matched_samples.csv'
            metadata_df.loc[list(common_samples)].to_csv(matched_metadata_file)
            print(f"Metadata for matched samples saved to: {matched_metadata_file}")

    else:
        print(f"Metadata file not found: {metadata_path}")    
    print("Processing complete!")

if __name__ == "__main__":
    main()
