"""
Sylph analysis toolkit for microbiome data analysis.

This module provides functions for processing and analyzing Sylph profiling data,
including diversity metrics, visualization, and differential abundance testing.

Usage:
    from sylph_tools import load_metadata, calculate_alpha_diversity, ...
"""

# Import core functions
from .diversity import (
    calculate_alpha_diversity,
    calculate_beta_diversity,
    compare_alpha_diversity
)

from .stats import (
    perform_permanova,
    differential_abundance_analysis
)

from .viz import (
    plot_alpha_diversity_boxplot,
    plot_ordination,
    plot_stacked_bar,
    plot_correlation_network,
    plot_relative_abundance_heatmap
)

from .utils import (
    load_metadata,
    preprocess_abundance_data,
    filter_low_abundance
)

__all__ = [
    'calculate_alpha_diversity',
    'calculate_beta_diversity',
    'compare_alpha_diversity',
    'perform_permanova',
    'differential_abundance_analysis',
    'plot_alpha_diversity_boxplot',
    'plot_ordination',
    'plot_stacked_bar',
    'plot_correlation_network',
    'plot_relative_abundance_heatmap',
    'load_metadata',
    'preprocess_abundance_data',
    'filter_low_abundance'
]