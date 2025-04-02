#!/usr/bin/env python3
"""
Transpose abundance data for co-occurrence analysis.
This script converts data from samples-as-rows to taxa-as-rows format.
"""

import pandas as pd
import sys

# Input and output files
input_file = "results/kraken_analysis/normalized_clr_abundance_with_metadata.tsv"
output_file = "results/data_prep/transposed_abundance.tsv"

try:
    # Read the input file
    print(f"Reading input file: {input_file}")
    data = pd.read_csv(input_file, sep='\t')
    
    # Extract metadata columns
    metadata_cols = ['SampleID', 'SubjectID', 'CollectionDate', 'Timing', 'Severity', 'Symptoms']
    metadata = data[metadata_cols].copy()
    
    # Extract abundance data
    abundance_cols = [col for col in data.columns if col not in metadata_cols]
    abundance = data[abundance_cols]
    
    # Transpose abundance data (taxa as rows, samples as columns)
    print("Transposing data...")
    transposed = abundance.T
    transposed.columns = data['SampleID']
    transposed.index.name = 'Species'
    
    # Save to output file
    print(f"Saving to: {output_file}")
    transposed.reset_index().to_csv(output_file, sep='\t', index=False)
    print(f"Transposed data shape: {transposed.shape}")
    print("Done!")
    
except Exception as e:
    print(f"Error: {e}")
    sys.exit(1)