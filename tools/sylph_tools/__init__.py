"""
Sylph Tools package for microbiome data analysis.

This package provides functions for processing and analyzing Sylph profiling data,
including diversity metrics, visualization, differential abundance testing, and more.

Main modules:
- sylph_diversity: Functions for calculating alpha and beta diversity metrics
- sylph_stats: Statistical analysis functions
- sylph_utils: Utility functions for data loading and preprocessing
- sylph_viz: Visualization functions
- sylph_tools: Main entry point exposing all key functions

Usage:
    import sylph_tools
    from sylph_tools import load_metadata, calculate_alpha_diversity, ...
"""

# Import core modules
from .sylph_diversity import (
    calculate_alpha_diversity,
    calculate_beta_diversity,
    compare_alpha_diversity
)

from .sylph_stats import (
    perform_permanova,
    differential_abundance_analysis,
    plot_abundance_boxplot,
    plot_beta_diversity_ordination
)

from .sylph_viz import (
    plot_alpha_diversity_boxplot,
    plot_ordination,
    plot_stacked_bar,
    plot_correlation_network,
    plot_relative_abundance_heatmap
)

from .sylph_utils import (
    load_metadata,
    preprocess_abundance_data,
    filter_low_abundance,
    create_abundance_summary,
    export_biom_format,
    add_taxonomy_metadata
)

# Re-export from sylph_tools for convenience and consistency
from .sylph_tools import *

# Define package version
__version__ = '0.1.0'

# Define all exportable names
__all__ = [
    # From sylph_diversity
    'calculate_alpha_diversity',
    'calculate_beta_diversity',
    'compare_alpha_diversity',
    
    # From sylph_stats
    'perform_permanova',
    'differential_abundance_analysis',
    'plot_abundance_boxplot',
    'plot_beta_diversity_ordination',
    
    # From sylph_viz
    'plot_alpha_diversity_boxplot',
    'plot_ordination',
    'plot_stacked_bar',
    'plot_correlation_network',
    'plot_relative_abundance_heatmap',
    
    # From sylph_utils
    'load_metadata',
    'preprocess_abundance_data',
    'filter_low_abundance',
    'create_abundance_summary',
    'export_biom_format',
    'add_taxonomy_metadata'
]
