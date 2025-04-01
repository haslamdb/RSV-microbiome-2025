#!/usr/bin/env python
# scripts/kraken/kraken_differential_abundance.py

import argparse
import os
import sys
import pandas as pd
import logging

# Add kraken_tools to path
kraken_tools_path = os.path.expanduser("~/Documents/Code/kraken_tools")
if os.path.exists(kraken_tools_path):
    sys.path.append(kraken_tools_path)
else:
    print(f"Error: kraken_tools directory not found at {kraken_tools_path}")
    sys.exit(1)

try:
    from kraken_tools.analysis.differential import run_differential_abundance_analysis
    from kraken_tools.logger import setup_logger, log_print
except ImportError as e:
    print(f"Error importing from kraken_tools: {e}")
    print("Make sure kraken_tools is installed and accessible")
    sys.exit(1)

def parse_arguments():
    parser = argparse.ArgumentParser(description="Run differential abundance analysis on Kraken2/Bracken data")
    
    parser.add_argument(
        "--abundance-file", 
        required=True,
        help="Path to abundance file (from process_kraken_data.py)"
    )
    
    parser.add_argument(
        "--metadata", 
        help="Path to metadata CSV file"
    )
    
    parser.add_argument(
        "--output-dir",
        default="results/differential_abundance",
        help="Directory for output files (default: results/differential_abundance)"
    )
    
    parser.add_argument(
        "--group-col", 
        default="Timing",
        help="Column name for grouping samples (default: Timing)"
    )
    
    parser.add_argument(
        "--methods",
        default="aldex2,ancom,ancom-bc",
        help="Comma-separated list of methods: aldex2, ancom, ancom-bc (default: all three)"
    )
    
    parser.add_argument(
        "--filter-groups",
        help="Comma-separated list of groups to include (default: use all groups)"
    )
    
    parser.add_argument(
        "--min-abundance", 
        type=float, 
        default=0.01, 
        help="Minimum relative abundance threshold (default: 0.01)"
    )
    
    parser.add_argument(
        "--min-prevalence", 
        type=float, 
        default=0.1, 
        help="Minimum prevalence threshold (default: 0.1)"
    )
    
    parser.add_argument(
        "--log-file",
        default=None,
        help="Path to log file (default: log to console only)"
    )
    
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Logging level (default: INFO)"
    )
    
    return parser.parse_args()

def main():
    # Parse command line arguments
    args = parse_arguments()
    
    # Setup logging
    logger = setup_logger(log_file=args.log_file, log_level=getattr(logging, args.log_level))
    log_print("Starting differential abundance analysis", level="info")
    
    # Check if abundance file exists
    if not os.path.exists(args.abundance_file):
        log_print(f"Error: Abundance file not found: {args.abundance_file}", level="error")
        sys.exit(1)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Parse methods
    methods = args.methods.split(',')
    log_print(f"Using methods: {', '.join(methods)}", level="info")
    
    # Parse filter groups if provided
    filter_groups = None
    if args.filter_groups:
        filter_groups = args.filter_groups.split(',')
        log_print(f"Filtering to groups: {', '.join(filter_groups)}", level="info")
    
    # Run differential abundance analysis
    try:
        # Load the abundance with metadata file which already contains both abundance and metadata
        abundance_data = pd.read_csv(args.abundance_file, sep='\t')
        
        # For our data, we know SampleID is the sample identifier, and other columns like Timing, Severity, Symptoms are metadata
        # Define known metadata columns
        metadata_cols = ['SampleID', 'SubjectID', 'CollectionDate', 'Timing', 'Severity', 'Symptoms']
        available_metadata_cols = [col for col in metadata_cols if col in abundance_data.columns]
        
        # Debugging info
        print(f"Columns in abundance_data: {', '.join(abundance_data.columns.tolist())}")
        print(f"Available metadata columns: {', '.join(available_metadata_cols)}")
        
        # Set SampleID as the index in both datasets to ensure proper matching
        abundance_data.set_index('SampleID', inplace=True)
        
        # All columns that are not metadata are taxonomy columns
        taxonomy_cols = [col for col in abundance_data.columns if col not in available_metadata_cols]
        
        # Create abundance dataframe with taxa as columns
        abundance_df = abundance_data[taxonomy_cols].copy()
        
        # Create metadata dataframe with only metadata columns, excluding SampleID which is now the index
        metadata_cols_minus_id = [col for col in available_metadata_cols if col != 'SampleID']
        metadata_df = abundance_data[metadata_cols_minus_id].copy()
        
        # Print sample index check
        print(f"Abundance samples: {len(abundance_df.index)}")
        print(f"Metadata samples: {len(metadata_df.index)}")
        
        log_print(f"Loaded abundance data with {len(abundance_df.index)} samples and {len(taxonomy_cols)} taxa", level="info")
        
        # Call the function with the correct signature
        results = run_differential_abundance_analysis(
            abundance_df=abundance_df,
            metadata_df=metadata_df,
            output_dir=args.output_dir,
            group_col=args.group_col,
            methods=methods,
            filter_groups=filter_groups
        )
        
        if results:
            log_print("Differential abundance analysis completed successfully", level="info")
            
            # Print summary of results
            for method, result in results.items():
                if isinstance(result, pd.DataFrame):
                    sig_count = sum(result.get('adjusted_p_value', result.get('q_value', result)) < 0.05)
                    log_print(f"Method {method}: {sig_count} significant taxa", level="info")
                else:
                    log_print(f"Method {method}: results available", level="info")
        else:
            log_print("Warning: No results produced from differential abundance analysis", level="warning")
    
    except Exception as e:
        log_print(f"Error during differential abundance analysis: {str(e)}", level="error")
        import traceback
        log_print(traceback.format_exc(), level="debug")
        sys.exit(1)
    
    log_print(f"Results saved to {args.output_dir}", level="info")

if __name__ == "__main__":
    main()