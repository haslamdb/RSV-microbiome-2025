"""
Utility functions for sylph data processing and analysis.
"""

import pandas as pd
import numpy as np
import os


def load_metadata(filepath, sample_id_column='SampleID'):
    """
    Load metadata from a CSV file.
    
    Parameters:
    -----------
    filepath : str
        Path to the metadata file
    sample_id_column : str
        Column name for sample IDs
        
    Returns:
    --------
    pandas.DataFrame
        Metadata DataFrame with sample IDs as index
    """
    try:
        metadata_df = pd.read_csv(filepath)
        
        # Check if the sample ID column exists
        if sample_id_column not in metadata_df.columns:
            raise ValueError(f"Sample ID column '{sample_id_column}' not found in metadata")
            
        # Set index and remove any duplicate sample IDs
        metadata_df = metadata_df.set_index(sample_id_column)
        if metadata_df.index.duplicated().any():
            print(f"Warning: Found {metadata_df.index.duplicated().sum()} duplicate sample IDs in metadata")
            metadata_df = metadata_df[~metadata_df.index.duplicated(keep='first')]
        
        # Convert categorical variables to string
        for col in metadata_df.columns:
            if metadata_df[col].dtype == 'object' or metadata_df[col].dtype.name == 'category':
                metadata_df[col] = metadata_df[col].astype(str)
        
        return metadata_df
    
    except Exception as e:
        print(f"Error loading metadata: {str(e)}")
        return pd.DataFrame()


def preprocess_abundance_data(abundance_df, normalize=True, log_transform=False, clr_transform=False):
    """
    Preprocess abundance data.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Taxa abundance DataFrame with taxa as index, samples as columns
    normalize : bool
        Whether to normalize to relative abundance
    log_transform : bool
        Whether to apply log transformation
    clr_transform : bool
        Whether to apply centered log-ratio transformation
        
    Returns:
    --------
    pandas.DataFrame
        Preprocessed abundance DataFrame
    """
    processed_df = abundance_df.copy()
    
    # Replace NaNs with zeros
    processed_df = processed_df.fillna(0)
    
    # Normalize to relative abundance
    if normalize:
        for col in processed_df.columns:
            col_sum = processed_df[col].sum()
            if col_sum > 0:
                processed_df[col] = processed_df[col] / col_sum * 100
    
    # Apply log transformation
    if log_transform:
        if clr_transform:
            # CLR transformation requires compositional data
            from skbio.stats.composition import clr
            
            # Add small pseudocount to zeros
            min_val = processed_df[processed_df > 0].min().min() / 2
            processed_df = processed_df.replace(0, min_val)
            
            # Apply CLR transformation (samples as rows)
            processed_df = pd.DataFrame(
                clr(processed_df.T.values),
                index=processed_df.columns,
                columns=processed_df.index
            ).T
        else:
            # Simple log transformation with pseudocount
            processed_df = np.log1p(processed_df)
    
    return processed_df


def filter_low_abundance(abundance_df, min_prevalence=0.1, min_abundance=0.01):
    """
    Filter out low abundance and low prevalence taxa.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Taxa abundance DataFrame with taxa as index, samples as columns
    min_prevalence : float
        Minimum fraction of samples in which a taxon must be present
    min_abundance : float
        Minimum mean relative abundance a taxon must have
        
    Returns:
    --------
    pandas.DataFrame
        Filtered abundance DataFrame
    """
    # Calculate prevalence (fraction of samples where taxon is present)
    prevalence = (abundance_df > 0).mean(axis=1)
    
    # Calculate mean abundance
    mean_abundance = abundance_df.mean(axis=1)
    
    # Filter based on thresholds
    keep_taxa = (prevalence >= min_prevalence) & (mean_abundance >= min_abundance)
    
    print(f"Filtering from {len(abundance_df)} to {keep_taxa.sum()} taxa")
    print(f"  Prevalence threshold: {min_prevalence:.2f} (must be present in {min_prevalence*100:.1f}% of samples)")
    print(f"  Abundance threshold: {min_abundance:.4f} (must have mean abundance ≥ {min_abundance*100:.2f}%)")
    
    return abundance_df.loc[keep_taxa]


def create_abundance_summary(abundance_df, metadata_df=None, group_var=None, top_n=20):
    """
    Create a summary table of abundance data.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Taxa abundance DataFrame with taxa as index, samples as columns
    metadata_df : pandas.DataFrame, optional
        Metadata DataFrame with samples as index
    group_var : str, optional
        Metadata variable to group by
    top_n : int
        Number of most abundant taxa to include
        
    Returns:
    --------
    pandas.DataFrame
        Summary table of abundance data
    """
    # Get mean abundance and prevalence
    mean_abundance = abundance_df.mean(axis=1)
    prevalence = (abundance_df > 0).mean(axis=1) * 100
    
    # Create summary DataFrame
    summary = pd.DataFrame({
        'Mean Abundance (%)': mean_abundance,
        'Prevalence (%)': prevalence
    })
    
    # Add group-specific means if metadata and group variable are provided
    if metadata_df is not None and group_var is not None and group_var in metadata_df.columns:
        # Get common samples
        common_samples = list(set(abundance_df.columns).intersection(set(metadata_df.index)))
        
        # Calculate mean abundance by group
        for group, group_df in metadata_df.loc[common_samples].groupby(group_var):
            group_samples = group_df.index
            group_mean = abundance_df[group_samples].mean(axis=1)
            summary[f'Mean in {group} (%)'] = group_mean
    
    # Sort by mean abundance and get top N
    summary = summary.sort_values('Mean Abundance (%)', ascending=False)
    
    if top_n is not None:
        summary = summary.head(top_n)
    
    return summary


def export_biom_format(abundance_df, output_file):
    """
    Export abundance data to BIOM format.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Taxa abundance DataFrame with taxa as index, samples as columns
    output_file : str
        Path to save BIOM file
        
    Returns:
    --------
    bool
        True if successful, False otherwise
    """
    try:
        from biom import Table
        from biom.util import biom_open
        
        # Create BIOM table
        table = Table(
            abundance_df.values,
            observation_ids=abundance_df.index,
            sample_ids=abundance_df.columns
        )
        
        # Write to file
        with biom_open(output_file, 'w') as f:
            table.to_hdf5(f, "Sylph abundance data")
        
        print(f"Exported abundance data to BIOM format: {output_file}")
        return True
    
    except Exception as e:
        print(f"Error exporting to BIOM format: {str(e)}")
        return False


def add_taxonomy_metadata(abundance_df):
    """
    Extract taxonomy information from taxon names and add as metadata.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Taxa abundance DataFrame with taxa as index, samples as columns
        
    Returns:
    --------
    pandas.DataFrame
        DataFrame with taxonomy information added
    """
    # Initialize taxonomy columns
    taxonomy_df = pd.DataFrame(index=abundance_df.index)
    taxonomy_df['Domain'] = 'Unknown'
    taxonomy_df['Genus'] = 'Unknown'
    taxonomy_df['Species'] = 'Unknown'
    taxonomy_df['Type'] = 'Unknown'
    
    # Extract taxonomy from taxon names
    for taxon in abundance_df.index:
        # Check for marker type (Bacteria, Virus, Fungus)
        if taxon.startswith('Bacteria:'):
            taxonomy_df.loc[taxon, 'Type'] = 'Bacteria'
            taxonomy_df.loc[taxon, 'Domain'] = 'Bacteria'
            name_part = taxon.split('Bacteria:')[1].strip()
        elif taxon.startswith('Virus:'):
            taxonomy_df.loc[taxon, 'Type'] = 'Virus'
            taxonomy_df.loc[taxon, 'Domain'] = 'Virus'
            name_part = taxon.split('Virus:')[1].strip()
        elif taxon.startswith('Fungus:'):
            taxonomy_df.loc[taxon, 'Type'] = 'Fungus'
            taxonomy_df.loc[taxon, 'Domain'] = 'Fungi'
            name_part = taxon.split('Fungus:')[1].strip()
        else:
            name_part = taxon
        
        # Try to extract genus and species
        name_parts = name_part.split()
        if len(name_parts) >= 2:
            # Check if first part starts with capital letter (likely genus)
            if name_parts[0][0].isupper():
                taxonomy_df.loc[taxon, 'Genus'] = name_parts[0]
                # Check if second part is lowercase (likely species)
                if len(name_parts) >= 2 and name_parts[1][0].islower():
                    taxonomy_df.loc[taxon, 'Species'] = ' '.join(name_parts[:2])
    
    # Join with original data
    result_df = abundance_df.copy()
    result_df = pd.concat([result_df, taxonomy_df], axis=1)
    
    return result_df
