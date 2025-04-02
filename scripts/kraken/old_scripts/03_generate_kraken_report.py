#!/usr/bin/env python3
"""
Generate comprehensive reports from Kraken2/Bracken analysis results.

This script:
1. Summarizes taxonomic profiles across samples and groups
2. Creates publication-ready tables and figures
3. Combines results from taxonomic and differential abundance analyses
4. Uses kraken_tools package for enhanced functionality

Usage:
    python scripts/kraken/03_generate_kraken_report.py [--results-dir RESULTS_DIR] [--output-dir OUTPUT_DIR]
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import glob

# Add project root to Python path
project_root = Path(__file__).resolve().parents[2]
sys.path.append(str(project_root))

# Add kraken_tools to Python path
kraken_tools_dir = Path.home() / "Documents" / "Code" / "kraken_tools"
sys.path.append(str(kraken_tools_dir))

# Import functions from kraken_tools
try:
    from kraken_tools.analysis.metadata import read_and_process_metadata
    from kraken_tools.analysis.taxonomy import read_and_process_taxonomy
    from kraken_tools.analysis.statistical import run_statistical_tests
    kraken_tools_available = True
except ImportError:
    kraken_tools_available = False
    print("Warning: kraken_tools module not available. Using built-in functions.")

# Set plotting style
plt.style.use('seaborn-v0_8-whitegrid')
sns.set(font_scale=1.1)


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Generate comprehensive Kraken/Bracken analysis reports')
    parser.add_argument('--results-dir', type=str, default='results/kraken_analysis',
                       help='Directory containing Kraken/Bracken analysis results')
    parser.add_argument('--diff-results-dir', type=str, default='results/kraken_differential',
                       help='Directory containing differential abundance results')
    parser.add_argument('--output-dir', type=str, default='results/kraken_report',
                       help='Directory to save generated reports')
    parser.add_argument('--abundance-file', type=str, default=None,
                       help='Path to abundance file (overrides automatic detection)')
    parser.add_argument('--metadata', type=str, default='metadata.csv',
                       help='Path to metadata file')
    parser.add_argument('--group-cols', type=str, default='Severity,Timing,Symptoms',
                       help='Comma-separated list of grouping variables')
    parser.add_argument('--top-taxa', type=int, default=20,
                       help='Number of top taxa to include in summary reports')
    parser.add_argument('--report-format', type=str, default='html',
                       choices=['html', 'pdf', 'both'],
                       help='Format for report output')
    return parser.parse_args()


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


def add_taxonomy_metadata(abundance_df):
    """
    Extract taxonomy information from Kraken/Bracken taxon names and add as metadata.
    
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
    taxonomy_df['Phylum'] = 'Unknown'
    taxonomy_df['Class'] = 'Unknown'
    taxonomy_df['Order'] = 'Unknown'
    taxonomy_df['Family'] = 'Unknown'
    taxonomy_df['Genus'] = 'Unknown'
    taxonomy_df['Species'] = 'Unknown'
    
    # Parse Kraken/Bracken taxonomy format
    for taxon in abundance_df.index:
        # For taxonomic path format (d__Bacteria|p__Firmicutes|...)
        if '|' in taxon:
            components = taxon.split('|')
            
            for component in components:
                # Extract rank and name
                if component.startswith('d__'):
                    taxonomy_df.loc[taxon, 'Domain'] = component[3:]
                elif component.startswith('p__'):
                    taxonomy_df.loc[taxon, 'Phylum'] = component[3:]
                elif component.startswith('c__'):
                    taxonomy_df.loc[taxon, 'Class'] = component[3:]
                elif component.startswith('o__'):
                    taxonomy_df.loc[taxon, 'Order'] = component[3:]
                elif component.startswith('f__'):
                    taxonomy_df.loc[taxon, 'Family'] = component[3:]
                elif component.startswith('g__'):
                    taxonomy_df.loc[taxon, 'Genus'] = component[3:]
                elif component.startswith('s__'):
                    taxonomy_df.loc[taxon, 'Species'] = component[3:]
        
        # For standard binomial format
        else:
            parts = taxon.split()
            if len(parts) >= 2:
                taxonomy_df.loc[taxon, 'Genus'] = parts[0]
                taxonomy_df.loc[taxon, 'Species'] = ' '.join(parts[:2])
    
    # Join with original data
    result_df = abundance_df.copy()
    result_df = pd.concat([result_df, taxonomy_df], axis=1)
    
    return result_df


def create_top_taxa_barplot(abundance_df, metadata_df=None, group_var=None, output_file=None, top_n=10):
    """
    Create a barplot of the top most abundant taxa.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Taxa abundance DataFrame with taxa as index, samples as columns
    metadata_df : pandas.DataFrame, optional
        Metadata DataFrame with samples as index
    group_var : str, optional
        Metadata variable to group by
    output_file : str, optional
        Path to save the plot
    top_n : int
        Number of top taxa to include
        
    Returns:
    --------
    matplotlib.figure.Figure
        Barplot figure
    """
    # Get the top N most abundant taxa
    mean_abundance = abundance_df.mean(axis=1)
    top_taxa = mean_abundance.nlargest(top_n).index.tolist()
    
    # Filter to top taxa
    top_abundance = abundance_df.loc[top_taxa]
    
    # Create plot based on whether metadata is available
    if metadata_df is not None and group_var is not None and group_var in metadata_df.columns:
        # Get common samples
        common_samples = list(set(abundance_df.columns).intersection(set(metadata_df.index)))
        
        # Calculate mean abundance by group
        group_means = {}
        for group, group_df in metadata_df.loc[common_samples].groupby(group_var):
            group_samples = group_df.index
            group_means[group] = top_abundance[group_samples].mean(axis=1)
        
        # Create DataFrame for plotting
        plot_data = pd.DataFrame(group_means)
        
        # Create stacked bar chart
        fig, ax = plt.subplots(figsize=(12, 8))
        plot_data.plot(kind='bar', stacked=False, ax=ax, width=0.7)
        
        # Set titles and labels
        ax.set_title(f'Top {top_n} Most Abundant Taxa by {group_var}')
        ax.set_ylabel('Mean Abundance (%)')
        ax.set_xlabel('Taxon')
        
        # Adjust legend and layout
        ax.legend(title=group_var, bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        
    else:
        # Just plot the overall mean abundance
        fig, ax = plt.subplots(figsize=(12, 6))
        mean_abundance.loc[top_taxa].plot(kind='bar', ax=ax)
        
        # Set titles and labels
        ax.set_title(f'Top {top_n} Most Abundant Taxa')
        ax.set_ylabel('Mean Abundance (%)')
        ax.set_xlabel('Taxon')
        
        # Adjust layout
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
    
    # Save figure if output file is provided
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
    
    return fig


def create_abundance_heatmap(abundance_df, metadata_df=None, group_var=None, output_file=None, top_n=30):
    """
    Create a heatmap of taxa abundance.
    
    Parameters:
    -----------
    abundance_df : pandas.DataFrame
        Taxa abundance DataFrame with taxa as index, samples as columns
    metadata_df : pandas.DataFrame, optional
        Metadata DataFrame with samples as index
    group_var : str, optional
        Metadata variable for grouping samples
    output_file : str, optional
        Path to save the plot
    top_n : int
        Number of top taxa to include
        
    Returns:
    --------
    matplotlib.figure.Figure
        Heatmap figure
    """
    # Get the top N most abundant taxa
    mean_abundance = abundance_df.mean(axis=1)
    top_taxa = mean_abundance.nlargest(top_n).index.tolist()
    
    # Filter to top taxa
    top_abundance = abundance_df.loc[top_taxa]
    
    # Create heatmap with sample clustering and optional group annotation
    if metadata_df is not None and group_var is not None and group_var in metadata_df.columns:
        # Get common samples
        common_samples = list(set(abundance_df.columns).intersection(set(metadata_df.index)))
        filtered_abundance = top_abundance[common_samples]
        
        # Get group information for sorting and annotation
        groups = metadata_df.loc[common_samples, group_var]
        
        # Define color map for groups
        unique_groups = sorted(groups.unique())
        group_colors = dict(zip(unique_groups, sns.color_palette("Set2", len(unique_groups))))
        
        # Create column colors for annotation
        col_colors = pd.Series(
            groups.map(group_colors),
            index=common_samples
        )
        
        # Create clustermap with group annotation
        g = sns.clustermap(
            filtered_abundance.T,  # Transpose so samples are rows
            cmap="YlGnBu",
            figsize=(14, 10),
            col_cluster=True,  # Cluster taxa (columns)
            row_cluster=True,  # Cluster samples (rows)
            col_colors=None,
            row_colors=col_colors,
            xticklabels=True,
            yticklabels=False,  # Too many samples usually
            cbar_kws={"label": "Abundance (%)"}
        )
        
        # Add title
        plt.suptitle(f"Abundance Heatmap of Top {top_n} Taxa", y=1.02)
        
        # Add group legend
        handles = [plt.Rectangle((0,0), 1, 1, color=color) for color in group_colors.values()]
        plt.legend(
            handles, 
            group_colors.keys(), 
            title=group_var, 
            bbox_to_anchor=(1, 1), 
            bbox_transform=plt.gcf().transFigure, 
            loc='upper right'
        )
        
        fig = g.fig
        
    else:
        # Simple heatmap without grouping
        fig, ax = plt.subplots(figsize=(14, 10))
        sns.heatmap(
            top_abundance,
            cmap="YlGnBu",
            ax=ax,
            cbar_kws={"label": "Abundance (%)"}
        )
        
        # Set titles and labels
        ax.set_title(f"Abundance Heatmap of Top {top_n} Taxa")
        ax.set_xlabel("Samples")
        ax.set_ylabel("Taxa")
        
        # Adjust layout
        plt.xticks(rotation=90)
        plt.tight_layout()
    
    # Save figure if output file is provided
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
    
    return fig


def summarize_differential_abundance(diff_dirs, min_significance=0.05):
    """
    Summarize results from differential abundance analyses.
    
    Parameters:
    -----------
    diff_dirs : list
        List of directories containing differential abundance results
    min_significance : float
        Significance threshold (adjusted p-value)
        
    Returns:
    --------
    pandas.DataFrame
        Summary of significant differential taxa
    """
    all_results = []
    
    # Process each directory
    for diff_dir in diff_dirs:
        # Find all TSV files with differential abundance results
        result_files = glob.glob(os.path.join(diff_dir, "**", "diff_abundance_*.tsv"), recursive=True)
        
        for file_path in result_files:
            try:
                # Extract variable name from filename
                variable = os.path.basename(file_path).replace("diff_abundance_", "").replace(".tsv", "")
                
                # Extract timing from directory path if available
                timing = "Overall"
                if "timing_" in file_path:
                    timing_match = re.search(r"timing_([^/]+)", file_path)
                    if timing_match:
                        timing = timing_match.group(1)
                
                # Load results
                results = pd.read_csv(file_path, sep='\t')
                
                # Filter to significant taxa
                if 'Significant' in results.columns:
                    sig_results = results[results['Significant'] == True]
                elif 'Adjusted P-value' in results.columns:
                    sig_results = results[results['Adjusted P-value'] < min_significance]
                else:
                    # Try other column names for p-value
                    p_value_cols = [col for col in results.columns if 'p' in col.lower() and 'value' in col.lower()]
                    if p_value_cols:
                        sig_results = results[results[p_value_cols[0]] < min_significance]
                    else:
                        print(f"Warning: Could not find p-value column in {file_path}")
                        continue
                
                # Add variable and timing information
                sig_results['Variable'] = variable
                sig_results['Timing'] = timing
                
                # Keep only essential columns
                keep_cols = ['Taxon', 'Variable', 'Timing', 'Adjusted P-value']
                
                # Add fold change column if available
                fold_change_cols = [col for col in sig_results.columns if 'fold' in col.lower() and 'change' in col.lower()]
                if fold_change_cols:
                    keep_cols.append(fold_change_cols[0])
                
                # Add test type if available
                if 'Test' in sig_results.columns:
                    keep_cols.append('Test')
                
                # Keep only columns that exist
                keep_cols = [col for col in keep_cols if col in sig_results.columns]
                
                # Append to all results
                all_results.append(sig_results[keep_cols])
                
            except Exception as e:
                print(f"Error processing {file_path}: {str(e)}")
    
    if not all_results:
        return pd.DataFrame()
    
    # Combine all results
    combined_results = pd.concat(all_results, ignore_index=True)
    
    # Sort by p-value
    if 'Adjusted P-value' in combined_results.columns:
        combined_results = combined_results.sort_values('Adjusted P-value')
    
    return combined_results


def create_venn_diagram(sig_taxa_dict, output_file=None):
    """
    Create a Venn diagram of significant taxa across groups.
    
    Parameters:
    -----------
    sig_taxa_dict : dict
        Dictionary mapping group names to sets of significant taxa
    output_file : str, optional
        Path to save the plot
        
    Returns:
    --------
    matplotlib.figure.Figure
        Venn diagram figure
    """
    try:
        from matplotlib_venn import venn2, venn3
        
        # Create figure
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Check number of groups
        if len(sig_taxa_dict) == 2:
            # Two-group Venn diagram
            groups = list(sig_taxa_dict.keys())
            venn = venn2(
                [set(sig_taxa_dict[groups[0]]), set(sig_taxa_dict[groups[1]])],
                set_labels=groups,
                ax=ax
            )
            
            # Set title
            ax.set_title(f"Overlap of Significant Taxa Between {groups[0]} and {groups[1]}")
            
        elif len(sig_taxa_dict) == 3:
            # Three-group Venn diagram
            groups = list(sig_taxa_dict.keys())
            venn = venn3(
                [set(sig_taxa_dict[groups[0]]), set(sig_taxa_dict[groups[1]]), set(sig_taxa_dict[groups[2]])],
                set_labels=groups,
                ax=ax
            )
            
            # Set title
            ax.set_title(f"Overlap of Significant Taxa Across {', '.join(groups)}")
            
        else:
            # Cannot create Venn diagram for more than 3 groups or fewer than 2
            if len(sig_taxa_dict) < 2:
                ax.text(0.5, 0.5, "Need at least 2 groups for Venn diagram", 
                       ha='center', va='center', fontsize=12)
            else:
                ax.text(0.5, 0.5, f"Cannot create Venn diagram for {len(sig_taxa_dict)} groups", 
                       ha='center', va='center', fontsize=12)
        
        # Save figure if output file is provided
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
        
        return fig
        
    except ImportError:
        print("Warning: matplotlib_venn not available. Cannot create Venn diagram.")
        
        # Create a simple error message plot
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.text(0.5, 0.5, "matplotlib_venn not available. Cannot create Venn diagram.",
               ha='center', va='center', fontsize=12)
        ax.axis('off')
        
        # Save figure if output file is provided
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
        
        return fig


def main():
    """Main function to generate comprehensive Kraken/Bracken analysis reports."""
    # Parse arguments
    args = parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    figures_dir = os.path.join(args.output_dir, "figures")
    tables_dir = os.path.join(args.output_dir, "tables")
    
    for dir_path in [figures_dir, tables_dir]:
        os.makedirs(dir_path, exist_ok=True)
    
    # Find abundance file if not specified
    abundance_file = args.abundance_file
    if not abundance_file:
        # Look for filtered abundance file first, then any abundance file
        filtered_files = glob.glob(os.path.join(args.results_dir, "filtered_kraken_*.tsv"))
        if filtered_files:
            abundance_file = filtered_files[0]
        else:
            abundance_files = glob.glob(os.path.join(args.results_dir, "*abundance*.tsv"))
            if abundance_files:
                abundance_file = abundance_files[0]
    
    if not abundance_file or not os.path.exists(abundance_file):
        print(f"Error: Could not find abundance file in {args.results_dir}")
        sys.exit(1)
    
    print(f"Using abundance file: {abundance_file}")
    
    # Load abundance data
    abundance_df = pd.read_csv(abundance_file, sep='\t', index_col=0)
    print(f"Loaded abundance data with {abundance_df.shape[0]} taxa and {abundance_df.shape[1]} samples")
    
    # Load metadata
    print(f"Loading metadata from {args.metadata}")
    if kraken_tools_available:
        metadata_df = read_and_process_metadata(args.metadata)
    else:
        metadata_df = pd.read_csv(args.metadata)
        metadata_df.set_index(args.metadata.split('=')[1] if '=' in args.metadata else 'SampleID', inplace=True)
    
    print(f"Loaded metadata with {metadata_df.shape[0]} samples and {metadata_df.shape[1]} columns")
    
    # Get grouping variables
    group_vars = args.group_cols.split(',')
    
    # Create abundance summaries
    print("\nGenerating abundance summaries...")
    overall_summary = create_abundance_summary(abundance_df, top_n=args.top_taxa)
    overall_summary.to_csv(os.path.join(tables_dir, "overall_abundance_summary.tsv"), sep='\t')
    print(f"Saved overall abundance summary to {os.path.join(tables_dir, 'overall_abundance_summary.tsv')}")
    
    # Create summaries for each grouping variable
    for group_var in group_vars:
        if group_var in metadata_df.columns:
            try:
                group_summary = create_abundance_summary(
                    abundance_df, 
                    metadata_df=metadata_df, 
                    group_var=group_var, 
                    top_n=args.top_taxa
                )
                
                # Save summary
                summary_file = os.path.join(tables_dir, f"abundance_summary_by_{group_var}.tsv")
                group_summary.to_csv(summary_file, sep='\t')
                print(f"Saved {group_var} abundance summary to {summary_file}")
                
                # Create barplot
                barplot_file = os.path.join(figures_dir, f"top_taxa_by_{group_var}.pdf")
                create_top_taxa_barplot(
                    abundance_df,
                    metadata_df=metadata_df,
                    group_var=group_var,
                    output_file=barplot_file,
                    top_n=10
                )
                print(f"Saved {group_var} barplot to {barplot_file}")
                
                # Create heatmap
                heatmap_file = os.path.join(figures_dir, f"abundance_heatmap_{group_var}.pdf")
                create_abundance_heatmap(
                    abundance_df,
                    metadata_df=metadata_df,
                    group_var=group_var,
                    output_file=heatmap_file,
                    top_n=30
                )
                print(f"Saved {group_var} heatmap to {heatmap_file}")
                
            except Exception as e:
                print(f"Error creating summaries for {group_var}: {str(e)}")
    
    # If differential abundance results are available, summarize them
    print("\nSummarizing differential abundance results...")
    if os.path.exists(args.diff_results_dir):
        # Find all subdirectories containing differential abundance results
        diff_dirs = [
            os.path.join(args.diff_results_dir, d) 
            for d in os.listdir(args.diff_results_dir) 
            if os.path.isdir(os.path.join(args.diff_results_dir, d))
        ]
        
        if diff_dirs:
            # Summarize differential abundance results
            diff_summary = summarize_differential_abundance(diff_dirs)
            
            if not diff_summary.empty:
                # Save summary
                diff_summary_file = os.path.join(tables_dir, "differential_abundance_summary.tsv")
                diff_summary.to_csv(diff_summary_file, sep='\t', index=False)
                print(f"Saved differential abundance summary to {diff_summary_file}")
                
                # Group significant taxa by variable
                sig_taxa_by_var = {}
                for var, group_df in diff_summary.groupby('Variable'):
                    sig_taxa_by_var[var] = group_df['Taxon'].tolist()
                
                # Create Venn diagram if multiple variables available
                if len(sig_taxa_by_var) >= 2:
                    venn_file = os.path.join(figures_dir, "significant_taxa_overlap.pdf")
                    create_venn_diagram(sig_taxa_by_var, output_file=venn_file)
                    print(f"Saved Venn diagram to {venn_file}")
            else:
                print("No significant differential abundance results found")
        else:
            print(f"No differential abundance result directories found in {args.diff_results_dir}")
    else:
        print(f"Differential abundance results directory not found: {args.diff_results_dir}")
    
    # Generate HTML report if requested
    if args.report_format in ['html', 'both']:
        try:
            # Use Jinja2 if available
            import jinja2
            print("\nGenerating HTML report...")
            
            # Create context data for template
            context = {
                'abundance_file': abundance_file,
                'n_taxa': abundance_df.shape[0],
                'n_samples': abundance_df.shape[1],
                'grouping_variables': group_vars,
                'top_taxa': overall_summary.index.tolist(),
                'report_date': pd.Timestamp.now().strftime('%Y-%m-%d %H:%M'),
                'has_diff_results': os.path.exists(args.diff_results_dir) and 'diff_summary' in locals()
            }
            
            # Add differential results if available
            if 'diff_summary' in locals() and not diff_summary.empty:
                context['n_diff_taxa'] = len(diff_summary['Taxon'].unique())
                context['n_diff_variables'] = len(diff_summary['Variable'].unique())
                
                # Get top differentially abundant taxa
                top_diff_taxa = diff_summary.head(10)['Taxon'].tolist()
                context['top_diff_taxa'] = top_diff_taxa
            
            # Create HTML report
            env = jinja2.Environment(
                loader=jinja2.FileSystemLoader("templates") if os.path.exists("templates") else jinja2.DictLoader({"report.html": DEFAULT_HTML_TEMPLATE})
            )
            template = env.get_template("report.html" if os.path.exists("templates/report.html") else "report.html")
            html_content = template.render(**context)
            
            # Save HTML report
            html_file = os.path.join(args.output_dir, "kraken_analysis_report.html")
            with open(html_file, 'w') as f:
                f.write(html_content)
            
            print(f"Saved HTML report to {html_file}")
            
        except ImportError:
            print("Warning: jinja2 not available. Cannot generate HTML report.")
    
    # Generate PDF report if requested
    if args.report_format in ['pdf', 'both']:
        try:
            # Use reportlab if available
            from reportlab.lib.pagesizes import letter
            from reportlab.pdfgen import canvas
            from reportlab.lib.styles import getSampleStyleSheet
            from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle
            
            print("\nGenerating PDF report...")
            
            # Create PDF report
            pdf_file = os.path.join(args.output_dir, "kraken_analysis_report.pdf")
            doc = SimpleDocTemplate(pdf_file, pagesize=letter)
            
            # Create content
            styles = getSampleStyleSheet()
            content = []
            
            # Add title
            content.append(Paragraph("Kraken2/Bracken Analysis Report", styles['Title']))
            content.append(Spacer(1, 12))
            
            # Add summary information
            content.append(Paragraph(f"Abundance File: {os.path.basename(abundance_file)}", styles['Normal']))
            content.append(Paragraph(f"Number of Taxa: {abundance_df.shape[0]}", styles['Normal']))
            content.append(Paragraph(f"Number of Samples: {abundance_df.shape[1]}", styles['Normal']))
            content.append(Spacer(1, 12))
            
            # Add top taxa
            content.append(Paragraph("Top Most Abundant Taxa:", styles['Heading2']))
            
            # Create a table of top taxa
            top_taxa_data = [["Taxon", "Mean Abundance (%)", "Prevalence (%)"]]
            for taxon in overall_summary.index[:10]:  # Top 10
                top_taxa_data.append([
                    taxon,
                    f"{overall_summary.loc[taxon, 'Mean Abundance (%)']:.2f}",
                    f"{overall_summary.loc[taxon, 'Prevalence (%)']:.1f}"
                ])
            
            # Add table to content
            top_taxa_table = Table(top_taxa_data)
            top_taxa_table.setStyle(TableStyle([
                ('BACKGROUND', (0, 0), (-1, 0), '#CCCCCC'),
                ('TEXTCOLOR', (0, 0), (-1, 0), '#000000'),
                ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
                ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
                ('BACKGROUND', (0, 1), (-1, -1), '#FFFFFF'),
                ('GRID', (0, 0), (-1, -1), 1, '#000000')
            ]))
            
            content.append(top_taxa_table)
            content.append(Spacer(1, 12))
            
            # Add figures if available
            for group_var in group_vars:
                barplot_file = os.path.join(figures_dir, f"top_taxa_by_{group_var}.pdf")
                if os.path.exists(barplot_file):
                    # Add section heading
                    content.append(Paragraph(f"Top Taxa by {group_var}", styles['Heading2']))
                    content.append(Spacer(1, 6))
                    
                    # Add image
                    img = Image(barplot_file, width=400, height=300)
                    content.append(img)
                    content.append(Spacer(1, 12))
            
            # Add differential abundance results if available
            if 'diff_summary' in locals() and not diff_summary.empty:
                content.append(Paragraph("Differential Abundance Results", styles['Heading2']))
                content.append(Spacer(1, 6))
                
                # Add summary text
                content.append(Paragraph(
                    f"Found {len(diff_summary['Taxon'].unique())} significantly different taxa "
                    f"across {len(diff_summary['Variable'].unique())} variables.",
                    styles['Normal']
                ))
                content.append(Spacer(1, 12))
                
                # Add Venn diagram if available
                venn_file = os.path.join(figures_dir, "significant_taxa_overlap.pdf")
                if os.path.exists(venn_file):
                    img = Image(venn_file, width=350, height=350)
                    content.append(img)
                    content.append(Spacer(1, 12))
            
            # Add footer with date
            content.append(Spacer(1, 30))
            content.append(Paragraph(
                f"Report generated on {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}",
                styles['Normal']
            ))
            
            # Build PDF document
            doc.build(content)
            
            print(f"Saved PDF report to {pdf_file}")
            
        except ImportError:
            print("Warning: reportlab not available. Cannot generate PDF report.")
    
    print("\nReport generation complete!")


# Default HTML template for reports
DEFAULT_HTML_TEMPLATE = """
<!DOCTYPE html>
<html>
<head>
    <title>Kraken2/Bracken Analysis Report</title>
    <style>
        body { font-family: Arial, sans-serif; max-width: 1200px; margin: 0 auto; padding: 20px; }
        h1, h2, h3 { color: #2C3E50; }
        .container { margin-top: 20px; }
        table { border-collapse: collapse; width: 100%; margin-bottom: 20px; }
        th, td { padding: 8px; text-align: left; border-bottom: 1px solid #ddd; }
        th { background-color: #f2f2f2; }
        tr:hover { background-color: #f5f5f5; }
    </style>
</head>
<body>
    <h1>Kraken2/Bracken Analysis Report</h1>
    <div class="container">
        <h2>Overview</h2>
        <p>File: {{ abundance_file }}</p>
        <p>Taxa: {{ n_taxa }}</p>
        <p>Samples: {{ n_samples }}</p>
        <p>Grouping Variables: {{ grouping_variables|join(', ') }}</p>
        <p>Report Date: {{ report_date }}</p>
    </div>
    
    <div class="container">
        <h2>Top Taxa</h2>
        <p>The most abundant taxa found in the samples:</p>
        <ul>
            {% for taxon in top_taxa[:10] %}
            <li>{{ taxon }}</li>
            {% endfor %}
        </ul>
    </div>
    
    {% if has_diff_results %}
    <div class="container">
        <h2>Differential Abundance</h2>
        <p>Found {{ n_diff_taxa }} taxa with significant differential abundance across {{ n_diff_variables }} variables.</p>
        
        {% if top_diff_taxa %}
        <h3>Top Differential Taxa</h3>
        <ul>
            {% for taxon in top_diff_taxa %}
            <li>{{ taxon }}</li>
            {% endfor %}
        </ul>
        {% endif %}
    </div>
    {% endif %}
    
    <div class="container">
        <h2>Figures</h2>
        <p>See the figures directory for detailed visualizations.</p>
    </div>
    
    <div class="container">
        <h2>Tables</h2>
        <p>See the tables directory for detailed data tables.</p>
    </div>
</body>
</html>
"""


if __name__ == "__main__":
    main()