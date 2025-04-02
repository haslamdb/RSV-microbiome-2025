#!/usr/bin/env python
# scripts/kraken/kraken_rf_shap.py

import argparse
import os
import sys
import logging

# Set up logging first
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler()]
)
logger = logging.getLogger('kraken_analysis')

# Add kraken_tools to path
kraken_tools_path = os.path.expanduser("~/Documents/Code/kraken_tools")
if os.path.exists(kraken_tools_path):
    sys.path.append(kraken_tools_path)
else:
    logger.error(f"Error: kraken_tools directory not found at {kraken_tools_path}")
    sys.exit(1)

try:
    from kraken_tools.analysis.rf_shap import run_rf_shap_analysis
    from kraken_tools.logger import setup_logger, log_print
except ImportError as e:
    logger.error(f"Error importing from kraken_tools: {e}")
    logger.error("Make sure kraken_tools is installed and accessible")
    sys.exit(1)

def parse_arguments():
    parser = argparse.ArgumentParser(description="Run Random Forest with SHAP analysis on Kraken2/Bracken data")
    
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
        default="results/rf_shap",
        help="Directory for output files (default: results/rf_shap)"
    )
    
    parser.add_argument(
        "--target-taxa", 
        help="Comma-separated list of taxa to analyze (default: top 10 most abundant)"
    )
    
    parser.add_argument(
        "--predictors", 
        default="Timing,Severity,Symptoms",
        help="Comma-separated list of predictor variables (default: Timing,Severity,Symptoms)"
    )
    
    parser.add_argument(
        "--random-effects", 
        default="SubjectID",
        help="Comma-separated list of random effect variables for mixed models (default: SubjectID)"
    )
    
    parser.add_argument(
        "--transform", 
        default="clr", 
        choices=["clr", "hellinger", "log", "none"],
        help="Transformation to apply to abundance data (default: clr)"
    )
    
    parser.add_argument(
        "--n-estimators", 
        type=int, 
        default=100, 
        help="Number of trees in Random Forest (default: 100)"
    )
    
    parser.add_argument(
        "--test-size", 
        type=float, 
        default=0.2, 
        help="Proportion of data for testing (default: 0.2)"
    )
    
    parser.add_argument(
        "--top-n", 
        type=int, 
        default=10, 
        help="Number of top taxa to analyze if target-taxa not specified (default: 10)"
    )
    
    parser.add_argument(
        "--mixed-model", 
        default="lmer", 
        choices=["lmer", "glmm"],
        help="Type of mixed model to use (default: lmer)"
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

def patch_rf_shap_module():
    """Monkey patch the run_rf_shap_analysis function to fix the logging issue"""
    import types
    import inspect
    from kraken_tools.analysis.rf_shap import run_rf_shap_analysis as orig_func
    
    source = inspect.getsource(orig_func)
    # Fix the line that's causing issues
    patched_source = source.replace(
        "logger = logging.getLogger('kraken_analysis')", 
        "import logging as logging_module\nlogger = logging_module.getLogger('kraken_analysis')"
    )
    
    # Create namespace for exec
    namespace = {}
    exec(patched_source, globals(), namespace)
    patched_func = namespace['run_rf_shap_analysis']
    
    # Replace the original function
    import kraken_tools.analysis.rf_shap
    kraken_tools.analysis.rf_shap.run_rf_shap_analysis = patched_func
    
    return patched_func

def main():
    # Parse command line arguments
    args = parse_arguments()
    
    # Setup logging with custom logger
    custom_logger = setup_logger(log_file=args.log_file, log_level=getattr(logging, args.log_level))
    log_print("Starting Random Forest with SHAP analysis", level="info")
    
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
    
    # Parse target taxa if provided
    target_taxa = None
    if args.target_taxa:
        target_taxa = args.target_taxa
        log_print(f"Analyzing specific taxa: {target_taxa}", level="info")
    else:
        log_print(f"Analyzing top {args.top_n} taxa by abundance", level="info")
    
    # Patch the rf_shap module to fix the logging issue
    patched_run_rf_shap = patch_rf_shap_module()
    
    # Run RF-SHAP analysis
    try:
        results = patched_run_rf_shap(
            abundance_file=args.abundance_file,
            metadata_file=args.metadata,
            output_dir=args.output_dir,
            target_taxa=target_taxa,
            predictors=args.predictors,
            random_effects=args.random_effects,
            transform=args.transform,
            n_estimators=args.n_estimators,
            test_size=args.test_size,
            top_n=args.top_n,
            mixed_model=args.mixed_model,
            log_file=args.log_file
        )
        
        if results:
            log_print("Random Forest with SHAP analysis completed successfully", level="info")
            log_print(f"Analyzed {len(results)} taxa", level="info")
        else:
            log_print("Warning: No results produced from RF-SHAP analysis", level="warning")
    
    except Exception as e:
        log_print(f"Error during RF-SHAP analysis: {str(e)}", level="error")
        import traceback
        log_print(traceback.format_exc(), level="debug")
        sys.exit(1)
    
    log_print(f"Results saved to {args.output_dir}", level="info")

if __name__ == "__main__":
    main()