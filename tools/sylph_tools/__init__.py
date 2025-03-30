"""
Functions for calculating alpha and beta diversity metrics for Sylph microbiome data.
"""

import pandas as pd
import numpy as np
from scipy import stats
from skbio.diversity import alpha_diversity, beta_diversity
from skbio.stats.distance import DistanceMatrix


def calculate_alpha_diversity(abundance_df, metrics=None):
    """
    Calculate alpha diversity metrics for each sample.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Taxa abundance DataFrame with taxa as index, samples as columns
    metrics : list, optional
        List of diversity metrics to calculate
        Default: ['shannon', 'simpson', 'observed_otus']
        
    Returns:
    --------
    pandas.DataFrame
        DataFrame with alpha diversity metrics for each sample
    """
    if metrics is None:
        metrics = ['shannon', 'simpson', 'observed_otus']
    
    # Initialize results DataFrame
    alpha_div = pd.DataFrame(index=abundance_df.columns)
    
    # Convert to relative abundance for proper diversity calculations
    rel_abundance = abundance_df.copy()
    sample_sums = rel_abundance.sum()
    for col in rel_abundance.columns:
        if sample_sums[col] > 0:
            rel_abundance[col] = rel_abundance[col] / sample_sums[col]
    
    # Calculate each reque