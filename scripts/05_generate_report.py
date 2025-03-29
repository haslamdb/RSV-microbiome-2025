#!/usr/bin/env python3
"""
Generate a comprehensive report of the microbiome analysis.

This script:
1. Compiles key findings from previous analyses
2. Creates summary visualizations
3. Generates a report with key statistics and figures
4. Summarizes diversity, differential abundance, and longitudinal results
5. Produces publication-ready figures and tables

Usage:
    python scripts/05_generate_report.py [--config CONFIG_FILE]
"""

import os
import sys
import argparse
import yaml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import datetime
import jinja2

# Add project root to Python path
project_root = Path(__file__).resolve().parents[1]
sys.path.append(str(project_root))

# Import functions from metaphlan_tools
from metaphlan_tools import (
    load_metadata,
    create_abundance_summary,
    plot_relative_abundance_heatmap,
    plot_stacked_bar
)

# Import for network plots
from metaphlan_tools.viz import plot_correlation_network


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Generate analysis report')
    parser.add_argument('--config', type=str, default='config/analysis_parameters.yml',
                       help='Path to configuration file')
    parser.add_argument('--html', action='store_true',
                       help='Generate HTML report')
    return parser.parse_args()


def load_analysis_results(results_dir):
    """
    Load previously generated analysis results.
    
    Parameters:
    -----------
    results_dir : Path
        Path to results directory
    
    Returns:
    --------
    dict
        Dictionary containing analysis results
    """
    results = {
        'alpha_diversity': None,
        'permanova': None,
        'differential_abundance': {},
        'temporal_patterns': None,
    }
    
    # Try to load alpha diversity
    alpha_file = results_dir / 'tables' / 'alpha_diversity.csv'
    if alpha_file.exists():
        results['alpha_diversity'] = pd.read_csv(alpha_file, index_col=0)
    
    # Try to load PERMANOVA results
    permanova_file = results_dir / 'tables' / 'permanova_results.csv'
    if permanova_file.exists():
        results['permanova'] = pd.read_csv(permanova_file, index_col=0)
    
    # Try to load differential abundance results
    diff_files = list(results_dir.glob('differential_abundance_*.csv'))
    for file in diff_files:
        var = file.stem.replace('differential_abundance_', '')
        results['differential_abundance'][var] = pd.read_csv(file)
    
    # Try to load temporal pattern results
    temporal_file = results_dir / 'tables' / 'temporal_pattern_stats.csv'
    if temporal_file.exists():
        results['temporal_patterns'] = pd.read_csv(temporal_file, index_col=0)
    
    return results


def create_summary_table(abundance_df, metadata_df, group_vars, top_n=20):
    """
    Create a comprehensive summary table of abundance data.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Species abundance DataFrame with species as index, samples as columns
    metadata_df : pandas.DataFrame
        Metadata DataFrame with samples as index
    group_vars : list
        List of metadata variables to group by
    top_n : int
        Number of most abundant species to include
        
    Returns:
    --------
    pandas.DataFrame
        Summary table of abundance data
    """
    # Create summary
    summary = create_abundance_summary(abundance_df, metadata_df, group_vars[0] if group_vars else None, top_n)
    
    # Add summary statistics for other group variables
    for var in group_vars[1:]:
        if var in metadata_df.columns:
            # Get common samples
            common_samples = set(abundance_df.columns).intersection(set(metadata_df.index))
            
            # Calculate mean abundance by group
            for group, group_df in metadata_df.loc[list(common_samples)].groupby(var):
                group_samples = group_df.index
                group_mean = abundance_df[group_samples].mean(axis=1)
                summary[f'Mean in {group} (%)'] = group_mean
    
    return summary


def generate_html_report(summary_stats, results, config, output_dir):
    """
    Generate an HTML report of the analysis results.
    
    Parameters:
    -----------
    summary_stats : dict
        Dictionary of summary statistics
    results : dict
        Dictionary of analysis results
    config : dict
        Configuration dictionary
    output_dir : Path
        Path to output directory
    
    Returns:
    --------
    None
    """
    # Create template loader
    template_loader = jinja2.FileSystemLoader(searchpath=project_root / 'templates')
    template_env = jinja2.Environment(loader=template_loader)
    
    # Check if templates directory exists, if not create it with a basic template
    templates_dir = project_root / 'templates'
    if not templates_dir.exists():
        templates_dir.mkdir(exist_ok=True)
        
        # Create a basic template
        with open(templates_dir / 'report_template.html', 'w') as f:
            f.write("""<!DOCTYPE html>
<html>
<head>
    <title>RSV Microbiome Analysis Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; line-height: 1.6; }
        h1 { color: #2c3e50; }
        h2 { color: #3498db; border-bottom: 1px solid #eee; padding-bottom: 5px; }
        table { border-collapse: collapse; width: 100%; margin-bottom: 20px; }
        th, td { text-align: left; padding: 8px; border: 1px solid #ddd; }
        th { background-color: #f2f2f2; }
        tr:nth-child(even) { background-color: #f9f9f9; }
        .figure { margin: 20px 0; text-align: center; }
        .figure img { max-width: 800px; border: 1px solid #ddd; }
        .figure-caption { font-style: italic; color: #666; }
        .metadata { background-color: #f8f9fa; padding: 15px; border-radius: 5px; margin-bottom: 20px; }
        .summary { background-color: #e8f4f8; padding: 15px; border-radius: 5px; margin-bottom: 20px; }
    </style>
</head>
<body>
    <h1>RSV Microbiome Analysis Report</h1>
    <p>Generated on {{ today }}</p>
    
    <div class="metadata">
        <h2>Dataset Summary</h2>
        <ul>
            <li><strong>Number of samples:</strong> {{ summary_stats.n_samples }}</li>
            <li><strong>Number of species:</strong> {{ summary_stats.n_species }}</li>
            <li><strong>Number of subjects:</strong> {{ summary_stats.n_subjects }}</li>
            <li><strong>Time points:</strong> {{ summary_stats.time_points|join(', ') }}</li>
        </ul>
    </div>
    
    <div class="summary">
        <h2>Key Findings</h2>
        <ul>
            {% for finding in summary_stats.key_findings %}
            <li>{{ finding }}</li>
            {% endfor %}
        </ul>
    </div>
    
    <h2>Diversity Analysis</h2>
    {% if 'alpha_diversity_findings' in summary_stats %}
    <h3>Alpha Diversity</h3>
    <p>{{ summary_stats.alpha_diversity_findings }}</p>
    {% endif %}
    
    {% if 'beta_diversity_findings' in summary_stats %}
    <h3>Beta Diversity</h3>
    <p>{{ summary_stats.beta_diversity_findings }}</p>
    {% endif %}
    
    <div class="figure">
        <img src="../figures/alpha_diversity_{{ config.metadata.group_variables[0] }}.png" alt="Alpha Diversity">
        <p class="figure-caption">Alpha diversity metrics across {{ config.metadata.group_variables[0] }} groups.</p>
    </div>
    
    <h2>Differential Abundance</h2>
    {% if 'differential_findings' in summary_stats %}
    <p>{{ summary_stats.differential_findings }}</p>
    {% endif %}
    
    <div class="figure">
        <img src="../figures/differential_species_heatmap.png" alt="Heatmap of differential species">
        <p class="figure-caption">Heatmap showing relative abundance of differentially abundant species.</p>
    </div>
    
    <h2>Longitudinal Analysis</h2>
    {% if 'longitudinal_findings' in summary_stats %}
    <p>{{ summary_stats.longitudinal_findings }}</p>
    {% endif %}
    
    <div class="figure">
        <img src="../figures/top_species_stacked_by_{{ config.metadata.time_variable }}.png" alt="Species over time">
        <p class="figure-caption">Changes in top species abundance over time.</p>
    </div>
    
    <h2>Appendix: Top Species</h2>
    <table>
        <tr>
            <th>Species</th>
            <th>Mean Abundance (%)</th>
            <th>Prevalence (%)</th>
            {% for group in summary_stats.groups %}
            <th>Mean in {{ group }} (%)</th>
            {% endfor %}
        </tr>
        {% for species, row in summary_stats.abundance_summary.iterrows() %}
        <tr>
            <td>{{ species }}</td>
            <td>{{ "%.2f"|format(row['Mean Abundance (%)']) }}</td>
            <td>{{ "%.1f"|format(row['Prevalence (%)']) }}</td>
            {% for group in summary_stats.groups %}
            <td>{{ "%.2f"|format(row['Mean in ' + group + ' (%)']) if 'Mean in ' + group + ' (%)' in row else '-' }}</td>
            {% endfor %}
        </tr>
        {% endfor %}
    </table>
</body>
</html>""")
    
    # Try to load the template
    try:
        template = template_env.get_template('report_template.html')
    except jinja2.exceptions.TemplateNotFound:
        print("Error: Template not found. Using simple report instead.")
        return
    
    # Render the template
    html_output = template.render(
        today=datetime.datetime.now().strftime("%Y-%m-%d"),
        summary_stats=summary_stats,
        config=config,
        results=results
    )
    
    # Write to file
    report_file = output_dir / 'report.html'
    with open(report_file, 'w') as f:
        f.write(html_output)
    
    print(f"HTML report saved to {report_file}")


def main():
    """Main function to generate analysis report."""
    # Parse arguments
    args = parse_args()
    
    # Load configuration
    config_path = project_root / args.config
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    # Set up paths
    processed_data_dir = project_root / 'data' / 'processed'
    results_dir = project_root / 'results'
    figures_dir = results_dir / 'figures'
    tables_dir = results_dir / 'tables'
    reports_dir = results_dir / 'reports'
    
    # Create directories if they don't exist
    figures_dir.mkdir(exist_ok=True, parents=True)
    tables_dir.mkdir(exist_ok=True, parents=True)
    reports_dir.mkdir(exist_ok=True, parents=True)
    
    # Load abundance data
    abundance_file = processed_data_dir / 'combined_abundance.csv'
    if not abundance_file.exists():
        print(f"Error: Abundance file not found at {abundance_file}")
        print("Please run 01_process_metaphlan_files.py first.")
        sys.exit(1)
    
    print(f"Loading abundance data from {abundance_file}")
    abundance_df = pd.read_csv(abundance_file, index_col=0)
    
    # Check if filtered abundance data exists
    filtered_file = processed_data_dir / 'filtered_abundance.csv'
    if filtered_file.exists():
        print(f"Loading filtered abundance data from {filtered_file}")
        filtered_abundance_df = pd.read_csv(filtered_file, index_col=0)
    else:
        filtered_abundance_df = abundance_df
    
    # Load metadata
    metadata_file = project_root / config['metadata']['filename']
    if not metadata_file.exists():
        print(f"Error: Metadata file not found at {metadata_file}")
        sys.exit(1)
    
    print(f"Loading metadata from {metadata_file}")
    metadata_df = load_metadata(metadata_file, config['metadata']['sample_id_column'])
    
    # Load previous analysis results
    print("\nLoading previous analysis results...")
    results = load_analysis_results(results_dir)
    
    # Create summary statistics
    print("\nGenerating summary statistics...")
    summary_stats = {}
    
    # Basic dataset statistics
    n_samples = len(set(abundance_df.columns).intersection(set(metadata_df.index)))
    n_species = len(abundance_df.index)
    
    subject_var = config['metadata']['subject_id_column']
    n_subjects = len(metadata_df[subject_var].unique()) if subject_var in metadata_df.columns else 'Unknown'
    
    time_var = config['metadata']['time_variable']
    time_points = list(metadata_df[time_var].unique()) if time_var in metadata_df.columns else ['Unknown']
    
    # Group variables
    group_vars = config['metadata']['group_variables']
    groups_summary = {}
    for var in group_vars:
        if var in metadata_df.columns:
            groups = list(metadata_df[var].unique())
            groups_summary[var] = {
                'groups': groups,
                'counts': metadata_df[var].value_counts().to_dict()
            }
    
    # Get significant species from differential abundance results
    sig_species = {}
    for var, df in results['differential_abundance'].items():
        if 'Adjusted P-value' in df.columns and 'Species' in df.columns:
            sig = df[df['Adjusted P-value'] < 0.05]['Species'].tolist()
            sig_species[var] = sig
    
    # Get temporal pattern species
    temporal_species = []
    if results['temporal_patterns'] is not None:
        temporal_df = results['temporal_patterns']
        if 'p-value' in temporal_df.columns:
            temporal_species = temporal_df[temporal_df['p-value'] < 0.05].index.tolist()
    
    # Create abundance summary
    print("Generating abundance summary...")
    abundance_summary = create_summary_table(filtered_abundance_df, metadata_df, group_vars, top_n=20)
    abundance_summary_file = tables_dir / 'abundance_summary.csv'
    abundance_summary.to_csv(abundance_summary_file)
    print(f"Abundance summary saved to {abundance_summary_file}")
    
    # Create correlation network
    print("\nGenerating species correlation network...")
    try:
        fig = plot_correlation_network(filtered_abundance_df, threshold=0.6, min_prevalence=0.3)
        network_file = figures_dir / 'correlation_network.png'
        fig.savefig(network_file, dpi=config['visualization']['figure_dpi'], bbox_inches='tight')
        plt.close(fig)
        print(f"Correlation network saved to {network_file}")
    except Exception as e:
        print(f"Error creating correlation network: {str(e)}")
    
    # Create global composition visualization
    print("\nGenerating global composition visualization...")
    try:
        # Create heatmap of top species
        print("Creating abundance heatmap...")
        top_n = min(30, len(filtered_abundance_df.index))
        mean_abundance = filtered_abundance_df.mean(axis=1)
        top_species = mean_abundance.nlargest(top_n).index.tolist()
        species_subset = filtered_abundance_df.loc[top_species]
        
        fig = plot_relative_abundance_heatmap(
            species_subset, 
            metadata_df, 
            group_vars[0] if group_vars else None,
            top_n=None,  # Already filtered to top species
            cluster_samples=True,
            cluster_taxa=True,
            cmap=config['visualization']['heatmap_colormap'],
            figsize=(12, 10)
        )
        
        heatmap_file = figures_dir / 'global_abundance_heatmap.png'
        fig.savefig(heatmap_file, dpi=config['visualization']['figure_dpi'], bbox_inches='tight')
        plt.close(fig)
        print(f"Global heatmap saved to {heatmap_file}")
        
        # Create stacked bar chart of composition by group
        for var in group_vars:
            if var in metadata_df.columns:
                print(f"Creating composition plot by {var}...")
                fig = plot_stacked_bar(
                    species_subset,
                    metadata_df,
                    var,
                    top_n=10,
                    other_category=True
                )
                
                bar_file = figures_dir / f'composition_by_{var}.png'
                fig.savefig(bar_file, dpi=config['visualization']['figure_dpi'], bbox_inches='tight')
                plt.close(fig)
                print(f"Composition plot saved to {bar_file}")
    except Exception as e:
        print(f"Error creating global composition visualization: {str(e)}")
    
    # Compile summary for report
    summary_stats = {
        'n_samples': n_samples,
        'n_species': n_species,
        'n_subjects': n_subjects,
        'time_points': time_points,
        'groups': groups_summary,
        'sig_species': sig_species,
        'temporal_species': temporal_species,
        'abundance_summary': abundance_summary.head(10),  # For HTML report
        'key_findings': []
    }
    
    # Add key findings
    if results['alpha_diversity'] is not None:
        summary_stats['alpha_diversity_findings'] = (
            f"Alpha diversity metrics were calculated for all samples. "
            f"The mean Shannon diversity index was {results['alpha_diversity']['Shannon'].mean():.2f}."
        )
        
        # Add to key findings if there are significant differences
        for var, groups in groups_summary.items():
            alpha_stats_file = tables_dir / f'alpha_diversity_{var}_stats.csv'
            if alpha_stats_file.exists():
                alpha_stats = pd.read_csv(alpha_stats_file, index_col=0)
                for metric in alpha_stats.index:
                    if 'p-value' in alpha_stats.columns and alpha_stats.loc[metric, 'p-value'] < 0.05:
                        finding = f"Significant differences in {metric} diversity were found between {var} groups (p={alpha_stats.loc[metric, 'p-value']:.3f})."
                        summary_stats['key_findings'].append(finding)
    
    if results['permanova'] is not None:
        summary_stats['beta_diversity_findings'] = (
            f"Beta diversity was calculated using the {config['diversity']['beta_metric']} distance metric. "
            f"PERMANOVA analysis was performed to test for differences in community composition."
        )
        
        # Add to key findings if there are significant differences
        for var in results['permanova'].index:
            if results['permanova'].loc[var, 'p-value'] < 0.05:
                finding = f"Significant differences in community composition were found between {var} groups (PERMANOVA p={results['permanova'].loc[var, 'p-value']:.3f})."
                summary_stats['key_findings'].append(finding)
    
    # Add differential abundance findings
    diff_counts = {var: len(species) for var, species in sig_species.items()}
    if any(diff_counts.values()):
        summary_stats['differential_findings'] = (
            f"Differential abundance testing identified species that differ significantly between clinical groups. "
            f"{', '.join([f'{count} species differed by {var}' for var, count in diff_counts.items() if count > 0])}."
        )
        
        # Add top species to key findings
        for var, species_list in sig_species.items():
            if species_list:
                top_species = species_list[:3]
                finding = f"Top differentially abundant species by {var}: {', '.join(top_species)}."
                summary_stats['key_findings'].append(finding)
    
    # Add longitudinal findings
    if len(temporal_species) > 0:
        summary_stats['longitudinal_findings'] = (
            f"Longitudinal analysis identified {len(temporal_species)} species with significant changes over {time_var}. "
            f"The most dynamic species included {', '.join(temporal_species[:3])}."
        )
        
        # Add to key findings
        finding = f"{len(temporal_species)} species showed significant temporal patterns across {time_var} points."
        summary_stats['key_findings'].append(finding)
    
    # Generate text report
    print("\nGenerating summary report...")
    report_file = reports_dir / 'summary_report.txt'
    with open(report_file, 'w') as f:
        f.write("======================================================\n")
        f.write(f"RSV MICROBIOME ANALYSIS SUMMARY REPORT\n")
        f.write(f"Generated on {datetime.datetime.now().strftime('%Y-%m-%d')}\n")
        f.write("======================================================\n\n")
        
        # Dataset summary
        f.write("DATASET SUMMARY\n")
        f.write("--------------\n")
        f.write(f"Number of samples: {n_samples}\n")
        f.write(f"Number of species: {n_species}\n")
        f.write(f"Number of subjects: {n_subjects}\n")
        f.write(f"Time points: {', '.join(str(t) for t in time_points)}\n\n")
        
        # Group summaries
        f.write("CLINICAL GROUPS\n")
        f.write("--------------\n")
        for var, info in groups_summary.items():
            f.write(f"{var} groups: {', '.join(str(g) for g in info['groups'])}\n")
            f.write(f"Distribution: {', '.join([f'{g}: {count}' for g, count in info['counts'].items()])}\n\n")
        
        # Key findings
        f.write("KEY FINDINGS\n")
        f.write("------------\n")
        for i, finding in enumerate(summary_stats['key_findings']):
            f.write(f"{i+1}. {finding}\n")
        f.write("\n")
        
        # Diversity analysis
        f.write("DIVERSITY ANALYSIS\n")
        f.write("-----------------\n")
        if 'alpha_diversity_findings' in summary_stats:
            f.write("Alpha Diversity:\n")
            f.write(f"{summary_stats['alpha_diversity_findings']}\n\n")
        
        if 'beta_diversity_findings' in summary_stats:
            f.write("Beta Diversity:\n")
            f.write(f"{summary_stats['beta_diversity_findings']}\n\n")
        
        # Differential abundance
        f.write("DIFFERENTIAL ABUNDANCE\n")
        f.write("---------------------\n")
        if 'differential_findings' in summary_stats:
            f.write(f"{summary_stats['differential_findings']}\n\n")
        
        f.write("Significantly different species by group:\n")
        for var, species_list in sig_species.items():
            if species_list:
                f.write(f"- {var}: {len(species_list)} species\n")
                for i, species in enumerate(species_list[:10]):  # Show top 10
                    df = results['differential_abundance'][var]
                    species_row = df[df['Species'] == species]
                    if not species_row.empty and 'Adjusted P-value' in species_row.columns:
                        p_val = species_row['Adjusted P-value'].iloc[0]
                        fold_change = species_row['Fold Change'].iloc[0] if 'Fold Change' in species_row.columns else 'N/A'
                        f.write(f"  {i+1}. {species} (adj. p={p_val:.3e}, fold change={fold_change})\n")
                    else:
                        f.write(f"  {i+1}. {species}\n")
                f.write("\n")
        
        # Longitudinal analysis
        f.write("LONGITUDINAL ANALYSIS\n")
        f.write("--------------------\n")
        if 'longitudinal_findings' in summary_stats:
            f.write(f"{summary_stats['longitudinal_findings']}\n\n")
        
        f.write("Species with significant temporal patterns:\n")
        if results['temporal_patterns'] is not None:
            temporal_df = results['temporal_patterns']
            if not temporal_df.empty and 'p-value' in temporal_df.columns:
                sig_temporal = temporal_df[temporal_df['p-value'] < 0.05].sort_values('p-value')
                for i, (species, row) in enumerate(sig_temporal.head(10).iterrows()):
                    f.write(f"  {i+1}. {species} (p={row['p-value']:.3e})\n")
        f.write("\n")
        
        # Top abundant species
        f.write("TOP ABUNDANT SPECIES\n")
        f.write("------------------\n")
        for i, (species, row) in enumerate(abundance_summary.head(15).iterrows()):
            f.write(f"{i+1}. {species} (mean: {row['Mean Abundance (%)']:.2f}%, prevalence: {row['Prevalence (%)']:.1f}%)\n")
        
        f.write("\n")
        f.write("======================================================\n")
        f.write("End of report\n")
    
    print(f"Summary report saved to {report_file}")
    
    # Generate HTML report if requested
    if args.html:
        print("\nGenerating HTML report...")
        generate_html_report(summary_stats, results, config, reports_dir)
    
    print("\nReport generation complete!")


if __name__ == "__main__":
    main()