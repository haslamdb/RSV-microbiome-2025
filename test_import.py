import sys
from pathlib import Path

# Print current path
print("Current Python path:")
for p in sys.path:
    print(f"  - {p}")

# Try different import approaches
print("\nTrying imports:")
try:
    import metaphlan_tools
    print("SUCCESS: import metaphlan_tools")
    print(f"Module location: {metaphlan_tools.__file__}")
except ImportError as e:
    print(f"FAILED: import metaphlan_tools - {e}")

project_root = Path(__file__).resolve().parent
tools_dir = project_root / 'tools'
sys.path.append(str(tools_dir))

try:
    import metaphlan_tools
    print("SUCCESS: import metaphlan_tools (after adding tools dir)")
    print(f"Module location: {metaphlan_tools.__file__}")
except ImportError as e:
    print(f"FAILED: import metaphlan_tools (after adding tools dir) - {e}")

# Check if the module exists as a directory
tools_path = project_root / 'tools' / 'metaphlan_tools'
print(f"\nChecking if module directory exists: {tools_path}")
print(f"Directory exists: {tools_path.exists()}")
if tools_path.exists():
    print(f"Contents: {list(tools_path.glob('*.py'))}")
