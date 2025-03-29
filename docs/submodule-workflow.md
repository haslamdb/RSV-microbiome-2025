# Workflow for Updating metaphlan_tools Submodule

This document outlines the steps to update the `metaphlan_tools` submodule in the `RSV-microbiome-2025` project when changes are made to the original `metaphlan_tools` repository.

## Initial Setup (One-time only)

If you haven't set up the submodule yet:

```bash
# Navigate to your main project
cd path/to/RSV-microbiome-2025

# Add the submodule
git submodule add https://github.com/yourusername/metaphlan_tools.git tools/metaphlan_tools

# Commit the addition
git commit -m "Add metaphlan_tools as submodule for microbiome analysis"
```

## Workflow for Updating the Submodule

### 1. Make Changes to metaphlan_tools

First, make your updates to the original `metaphlan_tools` repository:

```bash
# Navigate to the original repository
cd path/to/metaphlan_tools

# Make your changes
# ...

# Commit and push changes
git add .
git commit -m "Description of your changes"
git push origin main  # or whatever your branch name is
```

### 2. Update the Submodule in RSV-microbiome-2025

Now, update the submodule in your main project:

```bash
# Navigate to your main project
cd path/to/RSV-microbiome-2025

# Go to the submodule directory
cd tools/metaphlan_tools

# Fetch the latest changes
git fetch

# Check out the commit you want (usually the latest on main)
git checkout main  # or a specific tag/commit
git pull

# Go back to your main project root
cd ../..

# Commit the updated submodule reference
git add tools/metaphlan_tools
git commit -m "Update metaphlan_tools submodule to latest version"
git push
```

### Alternative Approach: Update from Main Project

You can also update the submodule directly from your main project:

```bash
# Navigate to your main project
cd path/to/RSV-microbiome-2025

# Update all submodules to latest commits on their tracked branches
git submodule update --remote --merge

# Commit the updated submodule reference
git add tools/metaphlan_tools
git commit -m "Update metaphlan_tools submodule to latest version"
git push
```

## For Collaborators: Cloning a Project with Submodules

When someone clones your RSV-microbiome-2025 repository, they need to initialize and update the submodules:

```bash
# Clone the main repository with submodules
git clone --recurse-submodules https://github.com/yourusername/RSV-microbiome-2025.git

# OR, if they already cloned without submodules:
git submodule init
git submodule update
```

## Testing Submodule Updates

After updating the submodule, it's good practice to test that everything works correctly:

```bash
# Navigate to your main project
cd path/to/RSV-microbiome-2025

# Install the submodule in development mode
pip install -e tools/metaphlan_tools

# Run tests if available
# ...

# Try a simple analysis to verify functionality
python -c "import metaphlan_tools; print(metaphlan_tools.__version__)"
```

## Updating to a Specific Version

If you need to update to a specific version (tag, branch, or commit):

```bash
cd path/to/RSV-microbiome-2025/tools/metaphlan_tools
git fetch
git checkout v1.0.0  # or a specific commit hash or branch name
cd ../..
git add tools/metaphlan_tools
git commit -m "Update metaphlan_tools submodule to version v1.0.0"
git push
```

## Common Issues and Troubleshooting

### Detached HEAD State

Submodules are often in a "detached HEAD" state, which is normal. If you need to make changes directly in the submodule:

1. First checkout a branch: `git checkout main`
2. Make your changes
3. Commit and push from within the submodule
4. Update the reference in the main project as described above

### Conflicts During Submodule Update

If you encounter conflicts during a submodule update:

1. Navigate to the submodule directory
2. Resolve conflicts normally using git
3. Commit the resolved conflicts
4. Continue the update process

### Forgot to Push Submodule Changes

If you updated the submodule locally but forgot to push its changes:

```bash
cd path/to/RSV-microbiome-2025/tools/metaphlan_tools
git push origin main  # or whatever branch you're on
cd ../..
git push
```
