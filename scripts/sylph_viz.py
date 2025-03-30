"""
Visualization functions for Sylph microbiome data.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from skbio.stats.ordination import pcoa
from sklearn.manifold import MDS


def plot_alpha_diversity_boxplot(alpha_df, metadata_df, group_var, metric=None):
    """
    Create a boxplot of alpha diversity by group.
    
    Parameters:
    -----------
    alpha_df : pandas.DataFrame
        Alpha diversity DataFrame with samples as index
    metadata_df : pandas.DataFrame
        Metadata DataFrame with samples as index
    group_var : str
        Grouping variable from metadata
    metric : str, optional
        Alpha diversity metric to plot (if None, plots all metrics)
        
    Returns:
    --------
    matplotlib.figure.Figure or dict
        Boxplot figure(s)
    """
    # Get common samples
    common_samples = list(set(alpha_df.index).intersection(set(metadata_df.index)))
    
    # Filter to common samples
    alpha_subset = alpha_df.loc[common_samples]
    metadata_subset = metadata_df.loc[common_samples]
    
    if metric is None:
        # Plot all metrics
        metrics = alpha_df.columns
        figures = {}
        
        for m in metrics:
            fig = _create_diversity_boxplot(alpha_subset, metadata_subset, group_var, m)
            figures[m] = fig
        
        return figures
    else:
        # Plot specific metric
        if metric not in alpha_df.columns:
            raise ValueError(f"Metric '{metric}' not found in alpha diversity data")
        
        return _create_diversity_boxplot(alpha_subset, metadata_subset, group_var, metric)


def _create_diversity_boxplot(alpha_df, metadata_df, group_var, metric):
    """Helper function to create a diversity boxplot."""
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Prepare data for plotting
    plot_data = pd.DataFrame({
        metric: alpha_df[metric],
        group_var: metadata_df[group_var]
    })
    
    # Create boxplot
    sns.boxplot(x=group_var, y=metric, data=plot_data, ax=ax)
    
    # Add individual points
    sns.stripplot(x=group_var, y=metric, data=plot_data, 
                 color='black', size=4, alpha=0.5, ax=ax)
    
    # Set titles and labels
    ax.set_title(f'{metric} Diversity by {group_var}')
    ax.set_xlabel(group_var)
    ax.set_ylabel(f'{metric} Diversity')
    
    # Rotate x-axis labels if needed
    plt.xticks(rotation=45 if len(str(plot_data[group_var].iloc[0])) > 10 else 0)
    
    # Adjust layout
    plt.tight_layout()
    
    return fig


def plot_ordination(beta_dm, metadata_df, variable, method='PCoA'):
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
    try:
        if method.upper() == 'PCOA':
            # Use principal coordinates analysis
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
            
        elif method.upper() == 'NMDS':
            # Use non-metric multidimensional scaling
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
            
            # Add stress value if available
            stress = getattr(mds, 'stress_', None)
            if stress is not None:
                ax.text(0.02, 0.98, f"Stress: {stress:.3f}", 
                       transform=ax.transAxes, va='top', ha='left',
                       bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        else:
            raise ValueError(f"Unknown ordination method: {method}. Use 'PCoA' or 'NMDS'.")
        
        # Adjust legend
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        
        # Adjust layout
        plt.tight_layout()
        
        return fig
    
    except Exception as e:
        # Create a simple error message plot
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.text(0.5, 0.5, f"Error creating {method} plot:\n{str(e)}",
               ha='center', va='center', fontsize=12)
        ax.set_title(f'{method} of Beta Diversity ({variable})')
        ax.axis('off')
        
        return fig


def plot_stacked_bar(abundance_df, metadata_df, group_var, top_n=10, other_category=True):
    """
    Create a stacked bar plot of the most abundant taxa by group.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Taxa abundance DataFrame with taxa as index, samples as columns
    metadata_df : pandas.DataFrame
        Metadata DataFrame with samples as index
    group_var : str
        Grouping variable from metadata
    top_n : int
        Number of top taxa to include
    other_category : bool
        Whether to include an "Other" category for remaining taxa
        
    Returns:
    --------
    matplotlib.figure.Figure
        Stacked bar plot figure
    """
    try:
        # Get common samples
        common_samples = list(set(abundance_df.columns).intersection(set(metadata_df.index)))
        
        # Filter to common samples
        filtered_abundance = abundance_df[common_samples]
        
        # Select taxa to display based on mean abundance
        mean_abundance = filtered_abundance.mean(axis=1)
        top_taxa = mean_abundance.nlargest(top_n).index.tolist()
        
        # Filter to top taxa
        top_abundance = filtered_abundance.loc[top_taxa]
        
        # Create "Other" category for remaining taxa if requested
        if other_category:
            other_abundance = filtered_abundance.drop(top_taxa).sum(axis=0)
            
            # Add "Other" to the data
            plot_data = top_abundance.copy()
            plot_data.loc['Other'] = other_abundance
        else:
            plot_data = top_abundance.copy()
        
        # Get group information
        group_info = metadata_df.loc[common_samples, group_var]
        
        # Calculate mean abundance by group
        group_means = {}
        for group in group_info.unique():
            group_samples = group_info[group_info == group].index
            group_means[group] = plot_data[group_samples].mean(axis=1)
        
        # Convert to DataFrame for plotting
        plot_df = pd.DataFrame(group_means)
        
        # Convert to percentages
        for col in plot_df.columns:
            plot_df[col] = plot_df[col] * 100 / plot_df[col].sum()
        
        # Create plot
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Create stacked bar chart
        plot_df.plot(kind='bar', stacked=True, ax=ax, colormap='tab20')
        
        # Set titles and labels
        ax.set_title(f'Mean Taxa Abundance by {group_var}')
        ax.set_xlabel(group_var)
        ax.set_ylabel('Relative Abundance (%)')
        
        # Adjust legend
        ax.legend(title='Taxa', bbox_to_anchor=(1.05, 1), loc='upper left')
        
        # Adjust layout
        plt.tight_layout()
        
        return fig
    
    except Exception as e:
        # Create a simple error message plot
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.text(0.5, 0.5, f"Error creating stacked bar plot:\n{str(e)}",
               ha='center', va='center', fontsize=12)
        ax.set_title(f'Taxa Abundance by {group_var}')
        ax.axis('off')
        
        return fig


def plot_correlation_network(abundance_df, threshold=0.6, min_prevalence=0.3, top_n=50):
    """
    Create a network plot of correlated taxa.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Taxa abundance DataFrame with taxa as index, samples as columns
    threshold : float
        Correlation threshold for including edges
    min_prevalence : float
        Minimum prevalence to include a taxon
    top_n : int
        Maximum number of taxa to include
        
    Returns:
    --------
    matplotlib.figure.Figure
        Network plot figure
    """
    try:
        import networkx as nx
        from scipy.stats import spearmanr
        
        # Filter by prevalence
        prevalence = (abundance_df > 0).mean(axis=1)
        filtered_df = abundance_df.loc[prevalence >= min_prevalence]
        
        # Limit to top abundant taxa if there are too many
        if len(filtered_df) > top_n:
            mean_abundance = filtered_df.mean(axis=1)
            top_taxa = mean_abundance.nlargest(top_n).index.tolist()
            filtered_df = filtered_df.loc[top_taxa]
        
        # Calculate correlation matrix
        corr_matrix = filtered_df.T.corr(method='spearman')
        
        # Create network
        G = nx.Graph()
        
        # Add nodes (taxa)
        for taxon in filtered_df.index:
            # Use mean abundance for node size
            abundance = filtered_df.loc[taxon].mean()
            G.add_node(taxon, abundance=abundance)
        
        # Add edges (correlations above threshold)
        for taxon1 in filtered_df.index:
            for taxon2 in filtered_df.index:
                if taxon1 >= taxon2:  # Avoid duplicates and self-loops
                    continue
                
                corr = corr_matrix.loc[taxon1, taxon2]
                
                # Only add edges for correlations above threshold
                if abs(corr) >= threshold:
                    G.add_edge(taxon1, taxon2, weight=abs(corr), sign=np.sign(corr))
        
        # Create plot
        fig, ax = plt.subplots(figsize=(12, 12))
        
        # Set positions using spring layout
        pos = nx.spring_layout(G, k=0.3, seed=42)
        
        # Draw nodes
        node_sizes = [G.nodes[node]['abundance'] * 1000 for node in G.nodes]
        nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color='skyblue', alpha=0.8)
        
        # Draw edges with different colors for positive and negative correlations
        pos_edges = [(u, v) for u, v, d in G.edges(data=True) if d['sign'] > 0]
        neg_edges = [(u, v) for u, v, d in G.edges(data=True) if d['sign'] < 0]
        
        edge_weights_pos = [G[u][v]['weight'] * 2 for u, v in pos_edges]
        edge_weights_neg = [G[u][v]['weight'] * 2 for u, v in neg_edges]
        
        nx.draw_networkx_edges(G, pos, edgelist=pos_edges, width=edge_weights_pos, 
                             edge_color='green', alpha=0.6)
        nx.draw_networkx_edges(G, pos, edgelist=neg_edges, width=edge_weights_neg, 
                             edge_color='red', alpha=0.6)
        
        # Draw node labels
        label_dict = {}
        for node in G.nodes:
            # Truncate long names
            if len(node) > 20:
                label_dict[node] = node.split()[0][:10] + '...'
            else:
                label_dict[node] = node
        
        nx.draw_networkx_labels(G, pos, labels=label_dict, font_size=8)
        
        # Add title and legend
        plt.title('Taxa Correlation Network', fontsize=15)
        
        # Create legend
        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D([0], [0], color='green', lw=2, label='Positive correlation'),
            Line2D([0], [0], color='red', lw=2, label='Negative correlation')
        ]
        ax.legend(handles=legend_elements, loc='upper left')
        
        # Add text with network statistics
        plt.text(0.02, 0.02, f"Nodes: {G.number_of_nodes()}\nEdges: {G.number_of_edges()}\nThreshold: {threshold}",
                transform=ax.transAxes, fontsize=10,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        # Remove axis
        plt.axis('off')
        
        return fig
    
    except Exception as e:
        # Create a simple error message plot
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.text(0.5, 0.5, f"Error creating correlation network:\n{str(e)}",
               ha='center', va='center', fontsize=12)
        ax.set_title('Taxa Correlation Network')
        ax.axis('off')
        
        return fig


def plot_relative_abundance_heatmap(abundance_df, metadata_df, group_var=None, top_n=30, 
                               cluster_samples=True, cluster_taxa=True, cmap='YlGnBu',
                               figsize=(12, 10)):
    """
    Create a heatmap of taxa abundance.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Taxa abundance DataFrame with taxa as index, samples as columns
    metadata_df : pandas.DataFrame
        Metadata DataFrame with samples as index
    group_var : str, optional
        Grouping variable for annotating samples
    top_n : int, optional
        Number of top taxa to include (if None, use all)
    cluster_samples : bool
        Whether to cluster samples
    cluster_taxa : bool
        Whether to cluster taxa
    cmap : str
        Colormap for heatmap
    figsize : tuple
        Figure size
        
    Returns:
    --------
    matplotlib.figure.Figure
        Heatmap figure
    """
    try:
        # Get common samples
        common_samples = list(set(abundance_df.columns).intersection(set(metadata_df.index)))
        
        # Filter to common samples
        filtered_abundance = abundance_df[common_samples]
        
        # Select top taxa if specified
        if top_n is not None and top_n < len(filtered_abundance):
            mean_abundance = filtered_abundance.mean(axis=1)
            top_taxa = mean_abundance.nlargest(top_n).index.tolist()
            plot_data = filtered_abundance.loc[top_taxa]
        else:
            plot_data = filtered_abundance
        
        # Ensure data is normalized for visualization
        plot_data = plot_data.apply(lambda x: x / x.max() if x.max() > 0 else x, axis=0)
        
        # Set up group annotations if provided
        if group_var is not None and group_var in metadata_df.columns:
            # Get groups for each sample
            groups = metadata_df.loc[common_samples, group_var]
            unique_groups = groups.unique()
            
            # Define group colors
            group_colors = dict(zip(unique_groups, sns.color_palette("Set2", len(unique_groups))))
            col_colors = pd.Series(groups).map(group_colors)
            
            # Create clustergrid with annotations
            g = sns.clustermap(
                plot_data,
                cmap=cmap,
                figsize=figsize,
                row_cluster=cluster_taxa,
                col_cluster=cluster_samples,
                col_colors=col_colors,
                xticklabels=False,
                cbar_kws={"label": "Relative Abundance"},
                dendrogram_ratio=(0.1, 0.1),
                colors_ratio=0.02,
            )
            
            # Add legend for groups
            handles = [plt.Rectangle((0,0), 1, 1, color=color) for color in group_colors.values()]
            plt.legend(handles, group_colors.keys(), title=group_var, 
                     loc='upper left', bbox_to_anchor=(1.01, 1.02))
        else:
            # Create clustergrid without annotations
            g = sns.clustermap(
                plot_data,
                cmap=cmap,
                figsize=figsize,
                row_cluster=cluster_taxa,
                col_cluster=cluster_samples,
                xticklabels=False,
                cbar_kws={"label": "Relative Abundance"}
            )
        
        # Adjust rotation of y labels
        plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)
        
        # Add title
        if top_n is not None:
            plt.suptitle(f"Top {len(plot_data.index)} Taxa by Abundance", y=1.02)
        else:
            plt.suptitle(f"Taxa Abundance Heatmap", y=1.02)
        
        return g.fig
    
    except Exception as e:
        # Create a simple error message plot
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.text(0.5, 0.5, f"Error creating heatmap:\n{str(e)}",
               ha='center', va='center', fontsize=12)
        ax.set_title('Taxa Abundance Heatmap')
        ax.axis('off')
        
        return fig
