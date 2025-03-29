import sys
import numpy as np
import pandas as pd
import sklearn
import metaphlan_tools

print(f"Python version: {sys.version}")
print(f"NumPy version: {np.__version__}")
print(f"pandas version: {pd.__version__}")
print(f"scikit-learn version: {sklearn.__version__}")
print(f"metaphlan_tools location: {metaphlan_tools.__file__}")

# Test basic functionality
print("\nTesting basic functionality...")
try:
    # Import key functions (adjust these to match your actual functions)
    from metaphlan_tools import parse_metaphlan_file, combine_samples, load_metadata
    print("Successfully imported core functions")
except Exception as e:
    print(f"Error importing functions: {e}")
