# Create a new file in your project directory:
# ~/Documents/Code/RSV-microbiome-2025/scripts/utils/parser_patch.py

import pandas as pd
from metaphlan_tools.parser import parse_metaphlan_file

def patched_combine_samples(files, sample_ids=None, taxonomic_level='species'):
    """
    Combine multiple MetaPhlAn output files into a single abundance table.
    This patched version handles duplicate species indices.
    
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
    
    if sample_ids is None:
        # Use file names as sample IDs
        sample_ids = [os.path.basename(f).split('.')[0] for f in files]
    
    if len(files) != len(sample_ids):
        raise ValueError("Number of files must match number of sample IDs")
    
    for i, file_path in enumerate(files):
        try:
            # Parse the MetaPhlAn file
            df = parse_metaphlan_file(file_path, taxonomic_level)
            
            # Set column name to sample ID
            df.columns = [sample_ids[i]]
            
            # Check for duplicate indices and handle them
            if df.index.duplicated().any():
                print(f"Warning: Found duplicate species in {sample_ids[i]}, keeping first occurrence")
                df = df[~df.index.duplicated(keep='first')]
            
            dfs.append(df)
        except Exception as e:
            print(f"Error processing {file_path}: {str(e)}")
            continue
    
    if not dfs:
        raise ValueError("No valid data frames to combine")
    
    # Combine along columns (each sample is a column)
    combined_df = pd.concat(dfs, axis=1)
    
    # Fill missing values with zeros
    combined_df = combined_df.fillna(0)
    
    return combined_df