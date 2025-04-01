#!/usr/bin/env python
# scripts/kraken/kraken_permanova.py

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
    from kraken_tools.analysis.permanova import run_permanova_analysis
    from kraken_tools.logger import setup_logger, log_print
except ImportError as e:
    print(f"Error importing from kraken_tools: {e}")
    print("Make sure kraken_tools is installed and accessible")
    sys.exit(1)

def parse_arguments():
    parser = argparse.ArgumentParser(description="Run PERMANOVA analysis on Kraken2/Bracken data")
    
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
        default="results/permanova",
        help="Directory for output files (default: results/permanova)"
    )
    
    parser.add_argument(
        "--categorical-vars", 
        default="Timing,Severity,Symptoms",
        help="Comma-separated list of categorical variables to test (default: Timing,Severity,Symptoms)"
    )
    
    parser.add_argument(
        "--distance-metric", 
        default="bray", 
        choices=["bray", "jaccard", "euclidean"],
        help="Distance metric to use (default: bray)"
    )
    
    parser.add_argument(
        "--transform", 
        default="clr", 
        choices=["clr", "hellinger", "log", "none"],
        help="Transformation to apply to abundance data (default: clr)"
    )
    
    parser.add_argument(
        "--permutations", 
        type=int, 
        default=999, 
        help="Number of permutations for significance testing (default: 999)"
    )
    
    parser.add_argument(
        "--min-group-size", 
        type=int, 
        default=3, 
        help="Minimum number of samples per group (default: 3)"
    )
    
    parser.add_argument(
        "--make-pcoa", 
        action="store_true", 
        default=True, 
        help="Generate PCoA plots (default: True)"
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
    log_print("Starting PERMANOVA analysis", level="info")
    
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
    
    # Parse categorical variables
    categorical_vars = args.categorical_vars.split(',')
    log_print(f"Testing categorical variables: {', '.join(categorical_vars)}", level="info")
    
    # Run PERMANOVA analysis
    try:
        results = run_permanova_analysis(
            abundance_file=args.abundance_file,
            metadata_file=args.metadata,
            output_dir=args.output_dir,
            categorical_vars=args.categorical_vars,
            distance_metric=args.distance_metric,
            transform=args.transform,
            permutations=args.permutations,
            min_group_size=args.min_group_size,
            make_pcoa=args.make_pcoa
        )
        
        if results:
            log_print("PERMANOVA analysis completed successfully", level="info")
            
            # Print summary of results
            for var, result in results.items():
                if 'p_value' in result:
                    sig = "significant" if result['p_value'] < 0.05 else "not significant"
                    log_print(f"Variable {var}: R2={result.get('r2', 'NA'):.3f}, p-value={result['p_value']:.4f} ({sig})", level="info")
        else:
            log_print("Warning: No results produced from PERMANOVA analysis", level="warning")
    
    except Exception as e:
        log_print(f"Error during PERMANOVA analysis: {str(e)}", level="error")
        import traceback
        log_print(traceback.format_exc(), level="debug")
        sys.exit(1)
    
    log_print(f"Results saved to {args.output_dir}", level="info")

if __name__ == "__main__":
    main()