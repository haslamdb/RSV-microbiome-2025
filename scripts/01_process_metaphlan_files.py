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

# Import our patched function
sys.path.append(str(project_root / 'scripts' / 'utils'))
from parser_patch import patched_combine_samples

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
    abundance_df = patched_combine_samples(files, sample_ids)
    
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
            metadata_df = load_metadata(metadata_path, config['metadata']['sample_id_column'])
            
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
    
    print("Processing complete!")

if __name__ == "__main__":
    main()
