#!/bin/bash
# scripts/kraken/run_analyses.sh

# Set base directories
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
CONFIG_FILE="$PROJECT_DIR/config/analysis_parameters.yml"
RESULTS_DIR="$PROJECT_DIR/results"

# Create directories if they don't exist
mkdir -p "$RESULTS_DIR"
mkdir -p "$RESULTS_DIR/kraken_analysis"

# Parse command line arguments
ANALYSIS="all"
ABUNDANCE_FILE=""
CUSTOM_METADATA=""

print_usage() {
    echo "Usage: $0 [OPTIONS]"
    echo "Run microbiome analyses for the RSV study"
    echo 
    echo "Options:"
    echo "  -a, --analysis TYPE      Type of analysis to run: diff-abundance, glmm, permanova,"
    echo "                           rf-shap, tsne, or all (default: all)"
    echo "  -f, --abundance-file FILE Path to abundance file (required if not using default)"
    echo "  -m, --metadata FILE      Path to metadata file (default from config)"
    echo "  -h, --help               Display this help message and exit"
    echo
    echo "Examples:"
    echo "  $0 --analysis diff-abundance --abundance-file results/kraken_analysis/filtered_abundance.tsv"
    echo "  $0 --analysis all"
}

# Process command line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -a|--analysis)
            ANALYSIS="$2"
            shift 2
            ;;
        -f|--abundance-file)
            ABUNDANCE_FILE="$2"
            shift 2
            ;;
        -m|--metadata)
            CUSTOM_METADATA="$2"
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

# Check if the kraken_tools directory exists
KRAKEN_TOOLS_DIR="$HOME/Documents/Code/kraken_tools"
if [ ! -d "$KRAKEN_TOOLS_DIR" ]; then
    echo "Error: kraken_tools directory not found at $KRAKEN_TOOLS_DIR"
    exit 1
fi

# Set default abundance file if not provided
if [ -z "$ABUNDANCE_FILE" ]; then
    ABUNDANCE_FILE="$RESULTS_DIR/kraken_analysis/filtered_abundance.tsv"
    echo "No abundance file provided, using default: $ABUNDANCE_FILE"
    
    # Check if the abundance file exists, if not, try to generate it
    if [ ! -f "$ABUNDANCE_FILE" ]; then
        echo "Warning: Default abundance file not found. Please run the processing step first or specify an abundance file."
        echo "Attempting to find alternative abundance files..."
        
        # Check if there are any abundance files in the results directory
        FOUND_FILES=$(find "$RESULTS_DIR" -name "*abundance*.tsv" -o -name "*abundance*.csv" | head -1)
        
        if [ -n "$FOUND_FILES" ]; then
            ABUNDANCE_FILE="$FOUND_FILES"
            echo "Found alternative abundance file: $ABUNDANCE_FILE"
        else
            echo "Error: No abundance files found. Please run the processing step first or specify an abundance file."
            exit 1
        fi
    fi
fi

# Set metadata parameter
METADATA_PARAM=""
if [ -n "$CUSTOM_METADATA" ]; then
    METADATA_PARAM="--metadata $CUSTOM_METADATA"
fi

# Check if the analysis type is valid
valid_analyses=("all" "diff-abundance" "glmm" "permanova" "rf-shap" "tsne")
if [[ ! " ${valid_analyses[*]} " =~ " $ANALYSIS " ]]; then
    echo "Error: Invalid analysis type: $ANALYSIS"
    print_usage
    exit 1
fi

# Create a timestamped log file
LOG_DIR="$RESULTS_DIR/logs"
mkdir -p "$LOG_DIR"
timestamp=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="$LOG_DIR/analysis_${ANALYSIS}_${timestamp}.log"

# Print run information
echo "=== RSV Microbiome Analysis ==="
echo "Analysis:       $ANALYSIS"
echo "Config file:    $CONFIG_FILE"
echo "Abundance file: $ABUNDANCE_FILE"
echo "Results dir:    $RESULTS_DIR"
echo "Log file:       $LOG_FILE"
echo "============================"

# Run the analysis
echo "Starting analysis... (check $LOG_FILE for progress)"
python "$SCRIPT_DIR/run_analyses.py" \
    --config "$CONFIG_FILE" \
    --abundance-file "$ABUNDANCE_FILE" \
    --output-dir "$RESULTS_DIR" \
    --analysis "$ANALYSIS" \
    $METADATA_PARAM \
    --log-file "$LOG_FILE" \
    --log-level "INFO"

exit_status=$?

if [ $exit_status -eq 0 ]; then
    echo "Analysis completed successfully. Results are in $RESULTS_DIR/$ANALYSIS"
    echo "Check the log file for details: $LOG_FILE"
else
    echo "Analysis failed with exit code $exit_status"
    echo "Check the log file for errors: $LOG_FILE"
    exit $exit_status
fi