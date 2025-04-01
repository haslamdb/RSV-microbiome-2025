# RSV-microbiome-2025 Guidelines

## Build/Run/Test Commands
- Install dependencies: `pip install -r requirements.txt`
- Run scripts sequentially: `bash scripts/kraken/kraken_rsv_microbiome.sh`
- Process Kraken data: `python scripts/kraken/01_process_kraken_data.py --config config/analysis_parameters.yml`
- Run differential abundance analysis: `python scripts/kraken/02_kraken_differential_abundance.py --abundance-file results/kraken_analysis/filtered_kraken_s_abundance.tsv`

## Code Style Guidelines
- PEP 8 style guide for Python
- Use docstrings for all functions and modules
- Import order: standard library, third-party, local modules
- Exception handling with specific exception types
- Variable naming: snake_case for variables/functions, PascalCase for classes
- Prefer pathlib over os.path for file operations
- Include robust error handling with informative messages
- Add type hints for function parameters when appropriate
- Use NumPy style docstrings for functions
- Wrap command-line tools in argparse interface