"""
Statistical analysis functions for Sylph microbiome data.
"""

import pandas as pd
import numpy as np
from scipy import stats
from skbio.stats.distance import permanova
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns


def perform_permanova(distance_matrix, metadata_df, variable, permutations=999):
    """
    Perform PERMANOVA test to see if grouping variable explains community differences.
    
    Parameters:
    -----------
    distance_matrix : skbio.DistanceMatrix
        Beta diversity distance matrix
    metadata_df : pandas.DataFrame
        Metadata DataFrame with samples as index
    variable : str
        Grouping variable in metadata
    permutations : int
        Number of permutations to use
        
    Returns:
    --------
    dict
        PERMANOVA results
    """
    # Filter metadata to include only samples in distance matrix
    common_samples = list(set(distance_matrix.ids).intersection(set(metadata_df.index)))
    
    if len(common_samples) < 5:
        return {
            'test-statistic': np.nan,
            'p-value': np.nan,
            'sample size': len(common_samples),
            'note': 'Insufficient samples for PERMANOVA'
        }
    
    # Filter distance matrix and metadata
    filtered_dm = distance_matrix.filter(common_samples)
    filtered_metadata = metadata_df.loc[common_samples]
    
    # Get grouping vector
    grouping = filtered_metadata[variable].astype(str).values
    
    # Check that we have at least two groups with 2+ samples
    unique_groups = np.unique(grouping)
    if len(unique_groups) < 2:
        return {
            'test-statistic': np.nan,
            'p-value': np.nan,
            'sample size': len(common_samples),
            'note': f'Only one group found in {variable}'
        }
    
    # Count samples per group
    valid_test = True
    for group in unique_groups:
        if np.sum(grouping == group) < 2:
            valid_test = False
            break
    
    if not valid_test:
        return {
            'test-statistic': np.nan,
            'p-value': np.nan,
            'sample size': len(common_samples),
            'note': f'At least one group in {variable} has fewer than 2 samples'
        }
    
    try:
        # Perform PERMANOVA
        results = permanova(filtered_dm, grouping, permutations=permutations)
        
        return {
            'test-statistic': results['test statistic'],
            'p-value': results['p-value'],
            'sample size': len(common_samples),
            'note': 'Successful test'
        }
    except Exception as e:
        return {
            'test-statistic': np.nan,
            'p-value': np.nan,
            'sample size': len(common_samples),
            'note': f'Error: {str(e)}'
        }


def differential_abundance_analysis(abundance_df, metadata_df, variable, method='wilcoxon'):
    """
    Identify differentially abundant taxa between groups.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Taxa abundance DataFrame with taxa as index, samples as columns
    metadata_df : pandas.DataFrame
        Metadata DataFrame with samples as index
    variable : str
        Grouping variable in metadata
    method : str
        Statistical method to use ('wilcoxon' or 'kruskal')
        
    Returns:
    --------
    pandas.DataFrame
        Results of differential abundance testing
    """
    # Get common samples
    common_samples = list(set(abundance_df.columns).intersection(set(metadata_df.index)))
    
    if len(common_samples) < 5:
        print(f"Error: Not enough samples for differential abundance testing (found {len(common_samples)})")
        return pd.DataFrame()
    
    # Filter to common samples
    filtered_abundance = abundance_df[common_samples]
    
    # Get groups
    groups = metadata_df.loc[common_samples, variable]
    unique_groups = groups.unique()
    
    if len(unique_groups) < 2:
        print(f"Error: Need at least 2 groups for differential abundance testing (found {len(unique_groups)})")
        return pd.DataFrame()
    
    # Determine which method to use based on number of groups
    if len(unique_groups) == 2 or method.lower() == 'wilcoxon':
        # Use Mann-Whitney U test for two groups
        return _two_group_differential_testing(filtered_abundance, groups, unique_groups)
    else:
        # Use Kruskal-Wallis test for multiple groups
        return _multi_group_differential_testing(filtered_abundance, groups, unique_groups)


def _two_group_differential_testing(abundance_df, groups, unique_groups):
    """
    Run differential abundance testing for two groups using Mann-Whitney U test.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Taxa abundance DataFrame with taxa as index, samples as columns
    groups : pandas.Series
        Group assignment for each sample
    unique_groups : array-like
        Unique group values
        
    Returns:
    --------
    pandas.DataFrame
        DataFrame with differential abundance results
    """
    # Initialize results
    results = []
    
    # Get sample indices for each group
    group1_samples = groups[groups == unique_groups[0]].index
    group2_samples = groups[groups == unique_groups[1]].index
    
    # Loop through each taxon
    for taxon in abundance_df.index:
        # Get abundance values for each group
        values1 = abundance_df.loc[taxon, group1_samples]
        values2 = abundance_df.loc[taxon, group2_samples]
        
        # Calculate mean abundance in each group
        mean1 = values1.mean()
        mean2 = values2.mean()
        
        # Calculate fold change (log2)
        # Add small pseudocount to avoid division by zero
        pseudocount = 1e-5
        fold_change = np.log2((mean2 + pseudocount) / (mean1 + pseudocount))
        
        # Calculate absolute fold change
        abs_fold_change = np.abs(fold_change)
        
        # Perform Mann-Whitney U test
        try:
            stat, p_value = stats.mannwhitneyu(values1, values2, alternative='two-sided')
        except Exception as e:
            # If test fails, set p-value to 1
            p_value = 1.0
            stat = np.nan
        
        # Calculate effect size (Cliff's Delta)
        try:
            cliff_delta = _calculate_cliffs_delta(values1, values2)
        except:
            cliff_delta = np.nan
        
        results.append({
            'Species': taxon,
            'P-value': p_value,
            'Group1': unique_groups[0],
            'Group2': unique_groups[1],
            'Mean in Group1': mean1,
            'Mean in Group2': mean2,
            'Log2 Fold Change': fold_change,
            'Abs Fold Change': abs_fold_change,
            'Cliff Delta': cliff_delta,
            'Test': 'Mann-Whitney U'
        })
    
    # Create DataFrame from results
    results_df = pd.DataFrame(results)
    
    # Apply multiple testing correction
    if not results_df.empty and len(results_df) > 1:
        try:
            # Use Benjamini-Hochberg procedure for FDR control
            results_df['Adjusted P-value'] = multipletests(
                results_df['P-value'], method='fdr_bh'
            )[1]
        except Exception as e:
            print(f"Error applying multiple testing correction: {str(e)}")
            results_df['Adjusted P-value'] = results_df['P-value']
    else:
        results_df['Adjusted P-value'] = results_df['P-value']
    
    # Sort by adjusted p-value
    results_df = results_df.sort_values('Adjusted P-value')
    
    # Rename fold change column for consistency with MetaPhlAn analysis
    results_df['Fold Change'] = 2**results_df['Log2 Fold Change']
    
    return results_df


def _multi_group_differential_testing(abundance_df, groups, unique_groups):
    """
    Run differential abundance testing for multiple groups using Kruskal-Wallis test.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Taxa abundance DataFrame with taxa as index, samples as columns
    groups : pandas.Series
        Group assignment for each sample
    unique_groups : array-like
        Unique group values
        
    Returns:
    --------
    pandas.DataFrame
        DataFrame with differential abundance results
    """
    # Initialize results
    results = []
    
    # Loop through each taxon
    for taxon in abundance_df.index:
        # Initialize group values and means
        group_values = {}
        group_means = {}
        
        # Get values and calculate means for each group
        for group in unique_groups:
            group_samples = groups[groups == group].index
            group_values[group] = abundance_df.loc[taxon, group_samples]
            group_means[group] = group_values[group].mean()
        
        # Perform Kruskal-Wallis test
        try:
            # Prepare data for the test
            all_values = [group_values[group] for group in unique_groups]
            stat, p_value = stats.kruskal(*all_values)
        except Exception as e:
            # If test fails, set p-value to 1
            p_value = 1.0
            stat = np.nan
        
        # Calculate maximum fold change between any two groups
        max_fold_change = 0
        max_group_pair = None
        
        for i, group1 in enumerate(unique_groups):
            for group2 in unique_groups[i+1:]:
                # Add small pseudocount to avoid division by zero
                pseudocount = 1e-5
                fold_change = np.abs(np.log2((group_means[group2] + pseudocount) / (group_means[group1] + pseudocount)))
                
                if fold_change > max_fold_change:
                    max_fold_change = fold_change
                    max_group_pair = (group1, group2)
        
        # Create result entry
        result = {
            'Species': taxon,
            'P-value': p_value,
            'Test': 'Kruskal-Wallis',
            'Max Log2 Fold Change': max_fold_change,
        }
        
        # Add mean for each group
        for group in unique_groups:
            result[f'Mean in {group}'] = group_means[group]
        
        # Add fold change information
        if max_group_pair:
            result['Max Fold Change Groups'] = f'{max_group_pair[0]} vs {max_group_pair[1]}'
            # Add actual fold change (non-log) for consistency with MetaPhlAn analysis
            result['Fold Change'] = 2**max_fold_change
        
        results.append(result)
    
    # Create DataFrame from results
    results_df = pd.DataFrame(results)
    
    # Apply multiple testing correction
    if not results_df.empty and len(results_df) > 1:
        try:
            # Use Benjamini-Hochberg procedure for FDR control
            results_df['Adjusted P-value'] = multipletests(
                results_df['P-value'], method='fdr_bh'
            )[1]
        except Exception as e:
            print(f"Error applying multiple testing correction: {str(e)}")
            results_df['Adjusted P-value'] = results_df['P-value']
    else:
        results_df['Adjusted P-value'] = results_df['P-value']
    
    # Sort by adjusted p-value
    results_df = results_df.sort_values('Adjusted P-value')
    
    return results_df


def _calculate_cliffs_delta(group1, group2):
    """
    Calculate Cliff's Delta effect size.
    
    Parameters:
    -----------
    group1 : array-like
        Values for first group
    group2 : array-like
        Values for second group
        
    Returns:
    --------
    float
        Cliff's Delta effect size
    """
    # Convert to arrays
    x = np.asarray(group1)
    y = np.asarray(group2)
    
    # Calculate the Cliff's Delta
    count_greater = 0
    count_less = 0
    
    for i in x:
        for j in y:
            if i > j:
                count_greater += 1
            elif i < j:
                count_less += 1
    
    # Calculate delta
    delta = (count_greater - count_less) / (len(x) * len(y))
    
    return delta


def plot_abundance_boxplot(abundance_df, metadata_df, species, group_var):
    """
    Create a boxplot of a species abundance by group.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Taxa abundance DataFrame with taxa as index, samples as columns
    metadata_df : pandas.DataFrame
        Metadata DataFrame with samples as index
    species : str
        Taxon name to plot
    group_var : str
        Grouping variable from metadata
        
    Returns:
    --------
    matplotlib.figure.Figure
        Boxplot figure
    """
    # Get common samples
    common_samples = list(set(abundance_df.columns).intersection(set(metadata_df.index)))
    
    # Get species abundance
    if species in abundance_df.index:
        species_abundance = abundance_df.loc[species, common_samples]
    else:
        raise ValueError(f"Species '{species}' not found in abundance data")
    
    # Get group information
    group_info = metadata_df.loc[common_samples, group_var]
    
    # Create data frame for plotting
    plot_data = pd.DataFrame({
        'Abundance': species_abundance,
        group_var: group_info
    })
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Create boxplot
    sns.boxplot(x=group_var, y='Abundance', data=plot_data, ax=ax)
    
    # Add individual points
    sns.stripplot(x=group_var, y='Abundance', data=plot_data, 
                 color='black', size=4, alpha=0.5, ax=ax)
    
    # Set titles and labels
    ax.set_title(f'{species} Abundance by {group_var}')
    ax.set_xlabel(group_var)
    ax.set_ylabel('Relative Abundance')
    
    # Adjust layout
    plt.tight_layout()
    
    return fig


def plot_beta_diversity_ordination(beta_dm, metadata_df, variable, method='PCoA'):
    """
    Create ordination plot from beta diversity distance matrix.
    
    Parameters:
    -----------
    beta_dm : skbio.DistanceMatrix
        Beta diversity distance matrix
    metadata_df : pandas.DataFrame
        Metadata DataFrame with samples as index
    variable : str
        Metadata variable for coloring points
    method : str
        Ordination method ('PCoA' or 'NMDS')
        
    Returns:
    --------
    matplotlib.figure.Figure
        Ordination plot figure
    """
    if method.upper() == 'PCOA':
        # Use principal coordinates analysis
        from skbio.stats.ordination import pcoa
        
        # Perform PCoA
        pcoa_results = pcoa(beta_dm)
        
        # Get the first two principal coordinates
        pc1 = pcoa_results.samples['PC1']
        pc2 = pcoa_results.samples['PC2']
        
        # Get variance explained
        variance_explained = pcoa_results.proportion_explained
        pc1_var = variance_explained[0] * 100
        pc2_var = variance_explained[1] * 100
        
        # Create a DataFrame for plotting
        plot_df = pd.DataFrame({
            'PC1': pc1,
            'PC2': pc2,
            'Sample': pcoa_results.samples.index
        })
        
        # Filter metadata to only include samples in the distance matrix
        common_samples = list(set(beta_dm.ids).intersection(set(metadata_df.index)))
        
        # Join with metadata to get the grouping variable
        plot_df = plot_df.set_index('Sample')
        plot_df[variable] = metadata_df.loc[common_samples, variable]
        
        # Create plot
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Add scatter plot with group colors
        sns.scatterplot(
            data=plot_df.reset_index(), 
            x='PC1', 
            y='PC2', 
            hue=variable, 
            s=100, 
            ax=ax
        )
        
        # Add axis labels with variance explained
        ax.set_xlabel(f'PC1 ({pc1_var:.1f}% variance explained)')
        ax.set_ylabel(f'PC2 ({pc2_var:.1f}% variance explained)')
        
        # Add title
        ax.set_title(f'PCoA of Beta Diversity ({variable})')
        
        # Adjust legend
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        
    elif method.upper() == 'NMDS':
        # Use non-metric multidimensional scaling
        from sklearn.manifold import MDS
        
        # Convert distance matrix to numpy array
        dist_array = beta_dm.data
        
        # Create MDS with non-metric scaling (this is essentially NMDS)
        mds = MDS(n_components=2, dissimilarity='precomputed', random_state=42, 
                 metric=False, n_init=10, max_iter=500)
        
        # Fit the model and transform
        coords = mds.fit_transform(dist_array)
        
        # Create a DataFrame for plotting
        plot_df = pd.DataFrame({
            'NMDS1': coords[:, 0],
            'NMDS2': coords[:, 1],
            'Sample': beta_dm.ids
        })
        
        # Filter metadata to only include samples in the distance matrix
        common_samples = list(set(beta_dm.ids).intersection(set(metadata_df.index)))
        
        # Join with metadata to get the grouping variable
        plot_df = plot_df.set_index('Sample')
        plot_df[variable] = metadata_df.loc[common_samples, variable]
        
        # Create plot
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Add scatter plot with group colors
        sns.scatterplot(
            data=plot_df.reset_index(), 
            x='NMDS1', 
            y='NMDS2', 
            hue=variable, 
            s=100, 
            ax=ax
        )
        
        # Add title
        ax.set_title(f'NMDS of Beta Diversity ({variable})')
        
        # Adjust legend
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        
        # Add stress value if available
        stress = getattr(mds, 'stress_', None)
        if stress is not None:
            ax.text(0.02, 0.98, f"Stress: {stress:.3f}", 
                   transform=ax.transAxes, va='top', ha='left',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    else:
        raise ValueError(f"Unknown ordination method: {method}. Use 'PCoA' or 'NMDS'.")
    
    # Adjust layout
    plt.tight_layout()
    
    return fig