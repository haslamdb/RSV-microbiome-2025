#!/usr/bin/env python
# scripts/kraken/run_analyses.py

import argparse
import os
import sys
import yaml
import logging
from pathlib import Path

# Add the kraken_tools repo to sys.path
kraken_tools_path = os.path.expanduser("~/Documents/Code/kraken_tools")
if os.path.exists(kraken_tools_path):
    sys.path.append(kraken_tools_path)
else:
    print(f"Error: kraken_tools directory not found at {kraken_tools_path}")
    sys.exit(1)

# Import functions from kraken_tools
try:
    from kraken_tools.analysis.differential import run_differential_abundance_analysis
    from kraken_tools.analysis.glmm_analysis import run_glmm_analysis
    from kraken_tools.analysis.permanova import run_permanova_analysis
    from kraken_tools.analysis.rf_shap import run_rf_shap_analysis
    from kraken_tools.analysis.tsne import run_tsne_analysis
    from kraken_tools.logger import setup_logger, log_print
except ImportError as e:
    print(f"Error importing from kraken_tools: {e}")
    print("Make sure kraken_tools is installed and accessible")
    sys.exit(1)

def parse_arguments():
    parser = argparse.ArgumentParser(description="Run microbiome analyses for RSV study")
    
    parser.add_argument(
        "--config", 
        required=True,
        help="Path to configuration file (YAML)"
    )
    
    parser.add_argument(
        "--abundance-file", 
        required=True,
        help="Path to taxonomic abundance file"
    )
    
    parser.add_argument(
        "--output-dir",
        default="results",
        help="Directory for output files (default: results)"
    )
    
    parser.add_argument(
        "--analysis",
        choices=["all", "diff-abundance", "glmm", "permanova", "rf-shap", "tsne"],
        default="all",
        help="Analysis to run (default: all)"
    )
    
    parser.add_argument(
        "--metadata",
        default=None,
        help="Path to metadata file (default: from config file)"
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

def load_config(config_path):
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
    return config

def run_diff_abundance(abundance_file, metadata_file, output_dir, config, logger):
    """Run differential abundance analysis"""
    log_print("Running differential abundance analysis...", level="info")
    
    # Create output directory
    diff_abundance_dir = os.path.join(output_dir, "differential_abundance")
    os.makedirs(diff_abundance_dir, exist_ok=True)
    
    # Get parameters from config
    min_abundance = config["differential_abundance"]["min_abundance"]
    min_prevalence = config["differential_abundance"]["min_prevalence"]
    group_col = config["metadata"]["group_variables"][0]  # Default to first group variable
    
    # Run differential abundance analysis
    try:
        results = run_differential_abundance_analysis(
            abundance_file=abundance_file,
            metadata_file=metadata_file,
            output_dir=diff_abundance_dir,
            group_col=group_col,
            methods=["aldex2", "ancom", "ancom-bc"],
            min_abundance=min_abundance,
            min_prevalence=min_prevalence
        )
        
        if results:
            log_print(f"Differential abundance analysis completed successfully. Results in {diff_abundance_dir}", level="info")
            return True
        else:
            log_print("Differential abundance analysis failed to produce results", level="error")
            return False
    except Exception as e:
        log_print(f"Error running differential abundance analysis: {str(e)}", level="error")
        return False

def run_glmm(abundance_file, metadata_file, output_dir, config, logger):
    """Run GLMM analysis"""
    log_print("Running GLMM analysis...", level="info")
    
    # Create output directory
    glmm_dir = os.path.join(output_dir, "glmm")
    os.makedirs(glmm_dir, exist_ok=True)
    
    # Get parameters from config
    min_abundance = config["differential_abundance"]["min_abundance"]
    min_prevalence = config["differential_abundance"]["min_prevalence"]
    subject_id_col = config["metadata"]["subject_id_column"]
    time_var = config["metadata"]["time_variable"]
    
    # Formula for GLMM
    formula = f"Count ~ {time_var} + (1|{subject_id_col})"
    
    # Run GLMM analysis
    try:
        success = run_glmm_analysis(
            abundance_file=abundance_file,
            metadata_file=metadata_file,
            output_dir=glmm_dir,
            formula=formula,
            model="negbin",
            min_abundance=min_abundance,
            min_prevalence=min_prevalence,
            logger=logger
        )
        
        if success:
            log_print(f"GLMM analysis completed successfully. Results in {glmm_dir}", level="info")
            return True
        else:
            log_print("GLMM analysis failed to produce results", level="error")
            return False
    except Exception as e:
        log_print(f"Error running GLMM analysis: {str(e)}", level="error")
        return False

def run_permanova(abundance_file, metadata_file, output_dir, config, logger):
    """Run PERMANOVA analysis"""
    log_print("Running PERMANOVA analysis...", level="info")
    
    # Create output directory
    permanova_dir = os.path.join(output_dir, "permanova")
    os.makedirs(permanova_dir, exist_ok=True)
    
    # Get parameters from config
    group_vars = config["metadata"]["group_variables"]
    categorical_vars = ",".join(group_vars)
    permutations = config["diversity"]["permanova_permutations"]
    beta_metric = config["diversity"]["beta_metric"]
    
    # Map config's beta_metric to PERMANOVA's distance_metric
    distance_metric_map = {
        "braycurtis": "bray",
        "jaccard": "jaccard",
        "euclidean": "euclidean"
    }
    distance_metric = distance_metric_map.get(beta_metric, "bray")
    
    # Run PERMANOVA analysis
    try:
        results = run_permanova_analysis(
            abundance_file=abundance_file,
            metadata_file=metadata_file,
            output_dir=permanova_dir,
            categorical_vars=categorical_vars,
            distance_metric=distance_metric,
            transform="clr",
            permutations=permutations,
            min_group_size=3,
            make_pcoa=True
        )
        
        if results:
            log_print(f"PERMANOVA analysis completed successfully. Results in {permanova_dir}", level="info")
            return True
        else:
            log_print("PERMANOVA analysis failed to produce results", level="error")
            return False
    except Exception as e:
        log_print(f"Error running PERMANOVA analysis: {str(e)}", level="error")
        return False

def run_rf_shap(abundance_file, metadata_file, output_dir, config, logger):
    """Run Random Forest with SHAP analysis"""
    log_print("Running Random Forest with SHAP analysis...", level="info")
    
    # Create output directory
    rf_shap_dir = os.path.join(output_dir, "rf_shap")
    os.makedirs(rf_shap_dir, exist_ok=True)
    
    # Get parameters from config
    subject_id_col = config["metadata"]["subject_id_column"]
    group_vars = config["metadata"]["group_variables"]
    predictors = ",".join(group_vars)
    
    # Run RF-SHAP analysis
    try:
        results = run_rf_shap_analysis(
            abundance_file=abundance_file,
            metadata_file=metadata_file,
            output_dir=rf_shap_dir,
            target_taxa=None,  # Will use top_n taxa
            predictors=predictors,
            random_effects=subject_id_col,
            transform="clr",
            n_estimators=100,
            test_size=0.2,
            top_n=10,
            mixed_model="lmer"
        )
        
        if results:
            log_print(f"RF-SHAP analysis completed successfully. Results in {rf_shap_dir}", level="info")
            return True
        else:
            log_print("RF-SHAP analysis failed to produce results", level="error")
            return False
    except Exception as e:
        log_print(f"Error running RF-SHAP analysis: {str(e)}", level="error")
        return False

def run_tsne(abundance_file, metadata_file, output_dir, config, logger):
    """Run t-SNE visualization"""
    log_print("Running t-SNE visualization...", level="info")
    
    # Create output directory
    tsne_dir = os.path.join(output_dir, "tsne")
    os.makedirs(tsne_dir, exist_ok=True)
    
    # Get parameters from config
    group_vars = config["metadata"]["group_variables"]
    categorical_vars = ",".join(group_vars)
    
    # Run t-SNE analysis
    try:
        output_path = run_tsne_analysis(
            abundance_file=abundance_file,
            metadata_file=metadata_file,
            output_dir=tsne_dir,
            target_taxa=None,  # Will use top taxa by abundance
            categorical_vars=categorical_vars,
            group_col=group_vars[0],  # Default to first group variable
            transform="clr",
            perplexity=30,
            n_iter=1000
        )
        
        if output_path:
            log_print(f"t-SNE visualization completed successfully. Results in {tsne_dir}", level="info")
            return True
        else:
            log_print("t-SNE visualization failed to produce results", level="error")
            return False
    except Exception as e:
        log_print(f"Error running t-SNE visualization: {str(e)}", level="error")
        return False

def main():
    # Parse command-line arguments
    args = parse_arguments()
    
    # Setup logging
    logger = setup_logger(log_file=args.log_file, log_level=getattr(logging, args.log_level))
    log_print("Starting RSV microbiome analyses", level="info")
    
    # Load configuration
    try:
        config = load_config(args.config)
    except Exception as e:
        log_print(f"Error loading configuration: {str(e)}", level="error")
        sys.exit(1)
    
    # Get metadata file path
    metadata_file = args.metadata if args.metadata else os.path.join(
        os.path.dirname(args.config), 
        "..", 
        config["metadata"]["filename"]
    )
    
    if not os.path.exists(metadata_file):
        log_print(f"Metadata file not found: {metadata_file}", level="error")
        sys.exit(1)
    
    # Check abundance file exists
    if not os.path.exists(args.abundance_file):
        log_print(f"Abundance file not found: {args.abundance_file}", level="error")
        sys.exit(1)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Run requested analyses
    analysis_functions = {
        "diff-abundance": run_diff_abundance,
        "glmm": run_glmm,
        "permanova": run_permanova,
        "rf-shap": run_rf_shap,
        "tsne": run_tsne
    }
    
    if args.analysis == "all":
        # Run all analyses
        for analysis_name, analysis_func in analysis_functions.items():
            success = analysis_func(
                args.abundance_file, 
                metadata_file, 
                args.output_dir, 
                config, 
                logger
            )
            if not success:
                log_print(f"Warning: {analysis_name} analysis failed or produced no results", level="warning")
    else:
        # Run the specific analysis
        success = analysis_functions[args.analysis](
            args.abundance_file, 
            metadata_file, 
            args.output_dir, 
            config, 
            logger
        )
        if not success:
            log_print(f"Analysis {args.analysis} failed or produced no results", level="error")
            sys.exit(1)
    
    log_print("RSV microbiome analyses completed", level="info")

if __name__ == "__main__":
    main()