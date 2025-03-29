# scripts/utils/parser_patch.py

import pandas as pd
import os
import csv
import re

def patched_parse_metaphlan_file(file_path, taxonomic_level='species'):
    """
    Parse a MetaPhlAn output file and extract abundances at a specific taxonomic level.
    This is a robust parser that handles various MetaPhlAn file formats and edge cases.
    
    Parameters:
    -----------
    file_path : str
        Path to the MetaPhlAn output file
    taxonomic_level : str, optional
        Taxonomic level to extract (default: 'species')
        
    Returns:
    --------
    pandas.DataFrame
        DataFrame with species as index and abundance as values
    """
    # Define taxonomic prefixes for filtering
    level_prefixes = {
        'kingdom': 'k__',
        'phylum': 'p__',
        'class': 'c__',
        'order': 'o__',
        'family': 'f__',
        'genus': 'g__',
        'species': 's__'
    }
    
    if taxonomic_level not in level_prefixes:
        raise ValueError(f"Invalid taxonomic level: {taxonomic_level}. "
                        f"Must be one of {list(level_prefixes.keys())}")
    
    target_prefix = level_prefixes[taxonomic_level]
    abundance_data = {}
    
    # Manually parse the file line by line for maximum robustness
    with open(file_path, 'r') as f:
        for line_num, line in enumerate(f, 1):
            # Skip comment lines except header comment
            if line.startswith('#') and not line.startswith('#clade_name'):
                continue
                
            # Skip empty lines
            if not line.strip():
                continue
            
            # Split by tabs
            parts = line.strip().split('\t')
            
            # Need at least 2 parts for parsing (taxonomy and abundance)
            if len(parts) < 2:
                print(f"Warning: Line {line_num} in {file_path} has fewer than 2 columns, skipping")
                continue
                
            # First column is taxonomy, last column is abundance (for safety)
            taxonomy = parts[0]
            try:
                # Try parsing with different column positions
                if len(parts) >= 3 and re.match(r'^[\d.]+$', parts[2]):
                    # Standard format: taxonomy, NCBI IDs, abundance
                    abundance = float(parts[2])
                elif len(parts) >= 2 and re.match(r'^[\d.]+$', parts[1]):
                    # Simplified format: taxonomy, abundance
                    abundance = float(parts[1])
                else:
                    # Last chance: try last column
                    abundance = float(parts[-1])
            except (ValueError, IndexError):
                print(f"Warning: Could not parse abundance in line {line_num} of {file_path}, skipping")
                continue
                
            # Check if this line contains the target taxonomic level
            if target_prefix in taxonomy:
                # Extract the name of the taxon at the target level
                taxon_parts = taxonomy.split('|')
                for part in taxon_parts:
                    if part.startswith(target_prefix):
                        taxon_name = part.replace(target_prefix, '')
                        abundance_data[taxon_name] = abundance
                        break
    
    if not abundance_data:
        raise ValueError(f"No {taxonomic_level}-level taxa found in {file_path}")
    
    # Convert to DataFrame
    df = pd.DataFrame(list(abundance_data.items()), columns=['Taxon', 'abundance'])
    df.set_index('Taxon', inplace=True)
    
    return df


def patched_combine_samples(files, sample_ids=None, taxonomic_level='species'):
    """
    Combine multiple MetaPhlAn output files into a single abundance table.
    This patched version handles duplicate species indices and various file formats.
    
    Parameters:
    -----------
    files : list
        List of file paths to MetaPhlAn output files
    sample_ids : list, optional
        List of sample IDs corresponding to each file
    taxonomic_level : str, optional
        Taxonomic level to extract (default: 'species')
        
    Returns:
    --------
    pandas.DataFrame
        Combined abundance table with species as rows and samples as columns
    """
    dfs = []
    successful_files = []
    successful_sample_ids = []
    
    if sample_ids is None:
        # Use file names as sample IDs
        sample_ids = [os.path.basename(f).split('.')[0] for f in files]
    
    if len(files) != len(sample_ids):
        raise ValueError("Number of files must match number of sample IDs")
    
    for i, file_path in enumerate(files):
        try:
            # Parse the MetaPhlAn file using our robust parser
            print(f"Processing file {i+1}/{len(files)}: {os.path.basename(file_path)}")
            df = patched_parse_metaphlan_file(file_path, taxonomic_level)
            
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
        # Print more helpful diagnostic information
        print("\nDiagnostic information:")
        print(f"Total files attempted: {len(files)}")
        print(f"Files that failed: {len(files)}")
        
        # Try to read the first few lines of the first file
        if files:
            try:
                print(f"\nFirst few lines of {files[0]}:")
                with open(files[0], 'r') as f:
                    for i, line in enumerate(f):
                        if i < 5:  # Print first 5 lines
                            print(f"Line {i+1}: {line.strip()}")
                        else:
                            break
            except Exception as e:
                print(f"Could not read file for diagnostics: {str(e)}")
                
        raise ValueError("No valid data frames to combine. Check file formats and error messages above.")
    
    print(f"\nSuccessfully processed {len(dfs)}/{len(files)} files")
    
    # Combine along columns (each sample is a column)
    combined_df = pd.concat(dfs, axis=1)
    
    # Fill missing values with zeros
    combined_df = combined_df.fillna(0)
    
    return combined_df


def diagnostic_file_check(file_path):
    """
    Perform a diagnostic check on a MetaPhlAn file to understand its structure.
    """
    try:
        print(f"\nDiagnostic check of {file_path}:")
        
        # Try to count lines and columns in the file
        with open(file_path, 'r') as f:
            lines = f.readlines()
            
        print(f"Total lines in file: {len(lines)}")
        
        # Check first 5 lines
        for i, line in enumerate(lines[:5]):
            stripped = line.strip()
            parts = stripped.split('\t')
            print(f"Line {i+1}: {len(parts)} columns, Content: {stripped[:80]}{'...' if len(stripped) > 80 else ''}")
            
        # Look for species lines
        species_count = 0
        for line in lines:
            if 's__' in line:
                species_count += 1
                
        print(f"Lines containing species entries (s__): {species_count}")
        
        return True
    except Exception as e:
        print(f"Diagnostic check failed: {str(e)}")
        return False
