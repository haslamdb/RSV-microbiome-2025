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

# Define package version
__version__ = '0.1.0'

# Import functions directly - using try/except for each module for robustness
try:
    from sylph_diversity import (
        calculate_alpha_diversity,
        calculate_beta_diversity,
        compare_alpha_diversity
    )
except ImportError:
    # Provide fallback or placeholders if needed
    def calculate_alpha_diversity(*args, **kwargs):
        print("Warning: calculate_alpha_diversity not available")
        return None
    
    def calculate_beta_diversity(*args, **kwargs):
        print("Warning: calculate_beta_diversity not available")
        return None
    
    def compare_alpha_diversity(*args, **kwargs):
        print("Warning: compare_alpha_diversity not available")
        return None

try:
    from sylph_stats import (
        perform_permanova,
        differential_abundance_analysis,
        plot_abundance_boxplot,
        plot_beta_diversity_ordination
    )
except ImportError:
    # Provide fallback or placeholders
    def perform_permanova(*args, **kwargs): 
        print("Warning: perform_permanova not available")
        return None
    
    def differential_abundance_analysis(*args, **kwargs):
        print("Warning: differential_abundance_analysis not available")
        return None
    
    def plot_abundance_boxplot(*args, **kwargs):
        print("Warning: plot_abundance_boxplot not available")
        return None
    
    def plot_beta_diversity_ordination(*args, **kwargs):
        print("Warning: plot_beta_diversity_ordination not available")
        return None

try:
    from sylph_viz import (
        plot_alpha_diversity_boxplot,
        plot_ordination,
        plot_stacked_bar,
        plot_correlation_network,
        plot_relative_abundance_heatmap
    )
except ImportError:
    # Provide fallback or placeholders
    def plot_alpha_diversity_boxplot(*args, **kwargs):
        print("Warning: plot_alpha_diversity_boxplot not available")
        return None
    
    def plot_ordination(*args, **kwargs):
        print("Warning: plot_ordination not available")
        return None
    
    def plot_stacked_bar(*args, **kwargs):
        print("Warning: plot_stacked_bar not available")
        return None
    
    def plot_correlation_network(*args, **kwargs):
        print("Warning: plot_correlation_network not available")
        return None
    
    def plot_relative_abundance_heatmap(*args, **kwargs):
        print("Warning: plot_relative_abundance_heatmap not available")
        return None

try:
    from sylph_utils import (
        load_metadata,
        preprocess_abundance_data,
        filter_low_abundance,
        create_abundance_summary,
        export_biom_format,
        add_taxonomy_metadata
    )
except ImportError:
    # Provide fallback or placeholders
    def load_metadata(*args, **kwargs):
        print("Warning: load_metadata not available")
        return None
    
    def preprocess_abundance_data(*args, **kwargs):
        print("Warning: preprocess_abundance_data not available")
        return None
    
    def filter_low_abundance(*args, **kwargs):
        print("Warning: filter_low_abundance not available")
        return None
    
    def create_abundance_summary(*args, **kwargs):
        print("Warning: create_abundance_summary not available")
        return None
    
    def export_biom_format(*args, **kwargs):
        print("Warning: export_biom_format not available")
        return None
    
    def add_taxonomy_metadata(*args, **kwargs):
        print("Warning: add_taxonomy_metadata not available")
        return None

# Try to import from sylph_tools, but okay if it fails
try:
    from sylph_tools import *
except ImportError:
    pass

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
