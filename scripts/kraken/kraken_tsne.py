#!/usr/bin/env python
# scripts/kraken/kraken_tsne.py

import argparse
import os
import sys
import logging

# Add kraken_tools to path
kraken_tools_path = os.path.expanduser("~/Documents/Code/kraken_tools")
if os.path.exists(kraken_tools_path):
    sys.path.append(kraken_tools_path)
else:
    print(f"Error: kraken_tools directory not found at {kraken_tools_path}")
    sys.exit(1)

try:
    from kraken_tools.analysis.tsne import run_tsne_analysis
    from kraken_tools.logger import setup_logger, log_print
except ImportError as e:
    print(f"Error importing from kraken_tools: {e}")
    print("Make sure kraken_tools is installed and accessible")
    sys.exit(1)

def parse_arguments():
    parser = argparse.ArgumentParser(description="Run t-SNE analysis on Kraken2/Bracken data")
    
    parser.add_argument(
        "--abundance-file", 
        required=True,
        help="Path to abundance file (from process_kraken_data.py)"
    )
    
    parser.add_argument(
        "--metadata", 
        required=True,
        help="Path to metadata CSV file"
    )
    
    parser.add_argument(
        "--output-dir",
        default="results/tsne",
        help="Directory for output files (default: results/tsne)"
    )
    
    parser.add_argument(
        "--target-taxa", 
        help="Comma-separated list of taxa to visualize (default: top 10 by abundance)"
    )
    
    parser.add_argument(
        "--categorical-vars", 
        default="Timing,Severity,Symptoms",
        help="Comma-separated list of categorical variables to visualize (default: Timing,Severity,Symptoms)"
    )
    
    parser.add_argument(
        "--group-col", 
        default="Timing",
        help="Primary grouping column for visualization (default: Timing)"
    )
    
    parser.add_argument(
        "--transform", 
        default="clr", 
        choices=["clr", "hellinger", "log", "none"],
        help="Transformation to apply to abundance data (default: clr)"
    )
    
    parser.add_argument(
        "--perplexity", 
        type=int, 
        default=30, 
        help="Perplexity parameter for t-SNE (default: 30)"
    )
    
    parser.add_argument(
        "--n-iter", 
        type=int, 
        default=1000, 
        help="Number of iterations for t-SNE (default: 1000)"
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
    log_print("Starting t-SNE analysis", level="info")
    
    # Check if abundance file exists
    if not os.path.exists(args.abundance_file):
        log_print(f"Error: Abundance file not found: {args.abundance_file}", level="error")
        sys.exit(1)
    
    # Check if metadata file exists
    if not os.path.exists(args.metadata):
        log_print(f"Error: Metadata file not found: {args.metadata}", level="error")
        sys.exit(1)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Run t-SNE analysis
    try:
        output_dir = run_tsne_analysis(
            abundance_file=args.abundance_file,
            metadata_file=args.metadata,
            output_dir=args.output_dir,
            target_taxa=args.target_taxa,
            categorical_vars=args.categorical_vars,
            group_col=args.group_col,
            transform=args.transform,
            perplexity=args.perplexity,
            n_iter=args.n_iter
        )
        
        if output_dir:
            log_print("t-SNE analysis completed successfully", level="info")
            log_print(f"Results saved to {output_dir}", level="info")
        else:
            log_print("Warning: No results produced from t-SNE analysis", level="warning")
    
    except Exception as e:
        log_print(f"Error during t-SNE analysis: {str(e)}", level="error")
        import traceback
        log_print(traceback.format_exc(), level="debug")
        sys.exit(1)

if __name__ == "__main__":
    main()