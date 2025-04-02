#!/usr/bin/env python
# run_permanova.py

import os
import subprocess
import sys

def main():
    """
    Run PERMANOVA analysis using the updated implementation with Hellinger transformation
    on filtered_S_abundance.tsv.
    """
    print("Running PERMANOVA analysis...")
    
    # Define paths
    abundance_file = "processed_data/filtered_S_abundance.tsv"
    metadata_file = "metadata.csv"
    output_dir = "results/permanova_hellinger"
    
    # Ensure the abundance file exists
    if not os.path.exists(abundance_file):
        print(f"Error: Abundance file not found: {abundance_file}")
        return 1
    
    # Ensure the metadata file exists
    if not os.path.exists(metadata_file):
        print(f"Error: Metadata file not found: {metadata_file}")
        return 1
    
    # Create command to run PERMANOVA
    cmd = [
        "python", "scripts/kraken/kraken_permanova.py",
        "--abundance-file", abundance_file,
        "--metadata", metadata_file,
        "--output-dir", output_dir,
        "--transform", "hellinger",
        "--distance-metric", "bray",
        "--categorical-vars", "Timing,Severity,Symptoms",
        "--permutations", "999"
    ]
    
    # Run the command
    try:
        print(f"Executing: {' '.join(cmd)}")
        result = subprocess.run(cmd, check=True)
        print(f"PERMANOVA analysis completed successfully. Results saved to {output_dir}")
        return 0
    except subprocess.CalledProcessError as e:
        print(f"Error running PERMANOVA analysis: {e}")
        return 1
    except Exception as e:
        print(f"Unexpected error: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main())