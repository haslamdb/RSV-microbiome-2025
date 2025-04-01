#!/bin/bash
# scripts/kraken/run_all_analyses.sh
# Master script to run all Kraken2/Bracken analyses

# Set base directories
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
CONFIG_FILE="$PROJECT_DIR/config/analysis_parameters.yml"
METADATA_FILE="$PROJECT_DIR/metadata.csv"
RESULTS_DIR="$PROJECT_DIR/results"
LOG_DIR="$RESULTS_DIR/logs"

# Create directories if they don't exist
mkdir -p "$RESULTS_DIR"
mkdir -p "$LOG_DIR"

# Generate timestamp for log files
timestamp=$(date +"%Y%m%d_%H%M%S")

# Function to print section headers
print_header() {
    echo
    echo "========================================"
    echo "  $1"
    echo "========================================"
    echo
}

# Check if the kraken_tools directory exists
KRAKEN_TOOLS_DIR="$HOME/Documents/Code/kraken_tools"
if [ ! -d "$KRAKEN_TOOLS_DIR" ]; then
    echo "Error: kraken_tools directory not found at $KRAKEN_TOOLS_DIR"
    exit 1
fi

# Parse command line arguments
KREPORT_DIR=""
BRACKEN_DIR=""
SKIP_PROCESSING=false
TAXONOMIC_LEVEL="S"
ANALYSES=()
NORMALIZE=true
NORMALIZATION_METHOD="clr"

print_usage() {
    echo "Usage: $0 [OPTIONS]"
    echo "Run Kraken2/Bracken analyses for RSV microbiome study"
    echo 
    echo "Options:"
    echo "  -k, --kreport-dir DIR     Directory containing Kraken2 report files"
    echo "  -b, --bracken-dir DIR     Directory containing Bracken abundance files"
    echo "  -s, --skip-processing     Skip the data processing step"
    echo "  -t, --taxonomic-level LVL Taxonomic level (D,P,C,O,F,G,S) (default: S)"
    echo "  -a, --analyses LIST       Comma-separated list of analyses to run:"
    echo "                           diff-abundance,glmm,permanova,rf-shap,tsne,all"
    echo "                           (default: all analyses)"
    echo "  -n, --no-normalize        Disable normalization"
    echo "  -N, --norm-method METHOD  Normalization method: relabundance,cpm,log10,clr"
    echo "                           (default: clr)"
    echo "  -m, --metadata FILE       Path to metadata file (default: project_dir/metadata.csv)"
    echo "  -h, --help                Display this help message and exit"
    echo
    echo "Examples:"
    echo "  $0 --kreport-dir data/kreports --bracken-dir data/bracken"
    echo "  $0 --skip-processing --analyses diff-abundance,permanova"
    echo "  $0 --norm-method log10 --analyses tsne,rf-shap"
}

# Process command line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -k|--kreport-dir)
            KREPORT_DIR="$2"
            shift 2
            ;;
        -b|--bracken-dir)
            BRACKEN_DIR="$2"
            shift 2
            ;;
        -s|--skip-processing)
            SKIP_PROCESSING=true
            shift
            ;;
        -t|--taxonomic-level)
            TAXONOMIC_LEVEL="$2"
            shift 2
            ;;
        -a|--analyses)
            IFS=',' read -r -a ANALYSES <<< "$2"
            shift 2
            ;;
        -n|--no-normalize)
            NORMALIZE=false
            shift
            ;;
        -N|--norm-method)
            NORMALIZATION_METHOD="$2"
            shift 2
            ;;
        -m|--metadata)
            METADATA_FILE="$2"
            shift 2
            ;;
        -h|--help)
            print_usage
            exit 0
            ;;
        *)
            echo "Error: Unknown option $1"
            print_usage
            exit 1
            ;;
    esac
done

# Print run information
print_header "RSV Microbiome Analysis with Kraken2/Bracken"
echo "Project directory: $PROJECT_DIR"
echo "Configuration file: $CONFIG_FILE"
echo "Metadata file: $METADATA_FILE"
echo "Results directory: $RESULTS_DIR"
echo "Log directory: $LOG_DIR"
echo "Timestamp: $timestamp"

if [ -n "$KREPORT_DIR" ]; then
    echo "Kraken2 report directory: $KREPORT_DIR"
fi
if [ -n "$BRACKEN_DIR" ]; then
    echo "Bracken directory: $BRACKEN_DIR"
fi
echo "Taxonomic level: $TAXONOMIC_LEVEL"
if [ "$NORMALIZE" = true ]; then
    echo "Normalization: Enabled (method: $NORMALIZATION_METHOD)"
else
    echo "Normalization: Disabled"
fi

# Check if at least one of kreport_dir or bracken_dir is provided when processing is needed
if [ "$SKIP_PROCESSING" = false ] && [ -z "$KREPORT_DIR" ] && [ -z "$BRACKEN_DIR" ]; then
    echo "Error: At least one of --kreport-dir or --bracken-dir must be provided"
    exit 1
fi

# Set up abundance file path
if [ "$NORMALIZE" = true ]; then
    ABUNDANCE_FILE="$RESULTS_DIR/kraken_analysis/normalized_abundance.tsv"
else
    ABUNDANCE_FILE="$RESULTS_DIR/kraken_analysis/filtered_abundance.tsv"
fi

# Step 1: Process Kraken/Bracken data
if [ "$SKIP_PROCESSING" = false ]; then
    print_header "Step 1: Processing Kraken2/Bracken Data"
    
    PROCESS_LOG="$LOG_DIR/process_kraken_${timestamp}.log"
    
    # Build command
    PROCESS_CMD="python $SCRIPT_DIR/process_kraken_data.py --config $CONFIG_FILE"
    
    if [ -n "$KREPORT_DIR" ]; then
        PROCESS_CMD="$PROCESS_CMD --kreport-dir $KREPORT_DIR"
    fi
    
    if [ -n "$BRACKEN_DIR" ]; then
        PROCESS_CMD="$PROCESS_CMD --bracken-dir $BRACKEN_DIR"
    fi
    
    PROCESS_CMD="$PROCESS_CMD --output-dir $RESULTS_DIR/kraken_analysis"
    PROCESS_CMD="$PROCESS_CMD --taxonomic-level $TAXONOMIC_LEVEL"
    PROCESS_CMD="$PROCESS_CMD --metadata $METADATA_FILE"
    
    # Add normalization options
    if [ "$NORMALIZE" = true ]; then
        PROCESS_CMD="$PROCESS_CMD --normalize --normalization-method $NORMALIZATION_METHOD"
    else
        PROCESS_CMD="$PROCESS_CMD --no-normalize"
    fi
    
    PROCESS_CMD="$PROCESS_CMD --log-file $PROCESS_LOG"
    
    echo "Running: $PROCESS_CMD"
    $PROCESS_CMD
    
    if [ $? -ne 0 ]; then
        echo "Error: Processing failed. Check log at $PROCESS_LOG"
        exit 1
    fi
    
    echo "Processing completed successfully."
else
    echo "Skipping data processing step as requested."
    
    # Check if abundance file exists
    if [ ! -f "$ABUNDANCE_FILE" ]; then
        echo "Warning: Abundance file $ABUNDANCE_FILE not found."
        
        # Try to find alternative files
        if [ "$NORMALIZE" = true ]; then
            # Look for normalized files first
            NORM_METHOD_FILE="$RESULTS_DIR/kraken_analysis/normalized_${NORMALIZATION_METHOD}_${TAXONOMIC_LEVEL}_abundance.tsv"
            if [ -f "$NORM_METHOD_FILE" ]; then
                ABUNDANCE_FILE="$NORM_METHOD_FILE"
                echo "Using alternative normalized abundance file: $ABUNDANCE_FILE"
            else
                # Try generic normalized file
                GENERIC_NORM_FILE="$RESULTS_DIR/kraken_analysis/normalized_abundance.tsv"
                if [ -f "$GENERIC_NORM_FILE" ]; then
                    ABUNDANCE_FILE="$GENERIC_NORM_FILE"
                    echo "Using generic normalized abundance file: $ABUNDANCE_FILE"
                fi
            fi
        fi
        
        # If still not found, try filtered files
        if [ ! -f "$ABUNDANCE_FILE" ]; then
            LEVEL_ABUNDANCE_FILE="$RESULTS_DIR/kraken_analysis/filtered_${TAXONOMIC_LEVEL}_abundance.tsv"
            if [ -f "$LEVEL_ABUNDANCE_FILE" ]; then
                ABUNDANCE_FILE="$LEVEL_ABUNDANCE_FILE"
                echo "Using filtered abundance file: $ABUNDANCE_FILE"
            else
                GENERIC_FILTERED_FILE="$RESULTS_DIR/kraken_analysis/filtered_abundance.tsv"
                if [ -f "$GENERIC_FILTERED_FILE" ]; then
                    ABUNDANCE_FILE="$GENERIC_FILTERED_FILE"
                    echo "Using generic filtered abundance file: $ABUNDANCE_FILE"
                else
                    echo "Error: No suitable abundance file found. Please run the processing step first."
                    exit 1
                fi
            fi
        fi
    fi
fi

# Determine which analyses to run
if [ ${#ANALYSES[@]} -eq 0 ]; then
    ANALYSES=("diff-abundance" "permanova" "rf-shap" "tsne")
elif [ "${ANALYSES[0]}" = "all" ]; then
    ANALYSES=("diff-abundance" "permanova" "rf-shap" "tsne")
fi

# Step 2: Run analyses
print_header "Step 2: Running Analyses"

# Run differential abundance analysis
if [[ " ${ANALYSES[*]} " =~ " diff-abundance " ]]; then
    print_header "Running Differential Abundance Analysis"
    
    DIFF_LOG="$LOG_DIR/diff_abundance_${timestamp}.log"
    
    DIFF_CMD="python $SCRIPT_DIR/kraken_differential_abundance.py"
    DIFF_CMD="$DIFF_CMD --abundance-file $ABUNDANCE_FILE"
    DIFF_CMD="$DIFF_CMD --metadata $METADATA_FILE"
    DIFF_CMD="$DIFF_CMD --output-dir $RESULTS_DIR/differential_abundance"
    DIFF_CMD="$DIFF_CMD --log-file $DIFF_LOG"
    
    echo "Running: $DIFF_CMD"
    $DIFF_CMD
    
    if [ $? -ne 0 ]; then
        echo "Warning: Differential abundance analysis failed. Check log at $DIFF_LOG"
    else
        echo "Differential abundance analysis completed successfully."
    fi
fi

# Run PERMANOVA analysis
if [[ " ${ANALYSES[*]} " =~ " permanova " ]]; then
    print_header "Running PERMANOVA Analysis"
    
    PERMANOVA_LOG="$LOG_DIR/permanova_${timestamp}.log"
    
    PERMANOVA_CMD="python $SCRIPT_DIR/kraken_permanova.py"
    PERMANOVA_CMD="$PERMANOVA_CMD --abundance-file $ABUNDANCE_FILE"
    PERMANOVA_CMD="$PERMANOVA_CMD --metadata $METADATA_FILE"
    PERMANOVA_CMD="$PERMANOVA_CMD --output-dir $RESULTS_DIR/permanova"
    PERMANOVA_CMD="$PERMANOVA_CMD --log-file $PERMANOVA_LOG"
    
    echo "Running: $PERMANOVA_CMD"
    $PERMANOVA_CMD
    
    if [ $? -ne 0 ]; then
        echo "Warning: PERMANOVA analysis failed. Check log at $PERMANOVA_LOG"
    else
        echo "PERMANOVA analysis completed successfully."
    fi
fi

# Run RF-SHAP analysis
if [[ " ${ANALYSES[*]} " =~ " rf-shap " ]]; then
    print_header "Running Random Forest with SHAP Analysis"
    
    RF_SHAP_LOG="$LOG_DIR/rf_shap_${timestamp}.log"
    
    RF_SHAP_CMD="python $SCRIPT_DIR/kraken_rf_shap.py"
    RF_SHAP_CMD="$RF_SHAP_CMD --abundance-file $ABUNDANCE_FILE"
    RF_SHAP_CMD="$RF_SHAP_CMD --metadata $METADATA_FILE"
    RF_SHAP_CMD="$RF_SHAP_CMD --output-dir $RESULTS_DIR/rf_shap"
    RF_SHAP_CMD="$RF_SHAP_CMD --log-file $RF_SHAP_LOG"
    
    echo "Running: $RF_SHAP_CMD"
    $RF_SHAP_CMD
    
    if [ $? -ne 0 ]; then
        echo "Warning: RF-SHAP analysis failed. Check log at $RF_SHAP_LOG"
    else
        echo "RF-SHAP analysis completed successfully."
    fi
fi

# Run t-SNE analysis
if [[ " ${ANALYSES[*]} " =~ " tsne " ]]; then
    print_header "Running t-SNE Analysis"
    
    TSNE_LOG="$LOG_DIR/tsne_${timestamp}.log"
    
    TSNE_CMD="python $SCRIPT_DIR/kraken_tsne.py"
    TSNE_CMD="$TSNE_CMD --abundance-file $ABUNDANCE_FILE"
    TSNE_CMD="$TSNE_CMD --metadata $METADATA_FILE"
    TSNE_CMD="$TSNE_CMD --output-dir $RESULTS_DIR/tsne"
    TSNE_CMD="$TSNE_CMD --log-file $TSNE_LOG"
    
    echo "Running: $TSNE_CMD"
    $TSNE_CMD
    
    if [ $? -ne 0 ]; then
        echo "Warning: t-SNE analysis failed. Check log at $TSNE_LOG"
    else
        echo "t-SNE analysis completed successfully."
    fi
fi

print_header "Analysis Pipeline Completed"
echo "All analyses have been completed."
echo "Results are available in: $RESULTS_DIR"
echo "Log files are available in: $LOG_DIR"