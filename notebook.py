{
 "cells": [
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "# Nasal Microbiome Diversity Analysis\n",
    "\n",
    "This notebook analyzes the microbial diversity of nasal microbiome samples in the RSV study.\n",
    "\n",
    "## Analysis steps:\n",
    "1. Load processed abundance data and metadata\n",
    "2. Calculate alpha diversity metrics\n",
    "3. Compare alpha diversity between clinical groups\n",
    "4. Calculate beta diversity and perform ordination\n",
    "5. Test for community composition differences using PERMANOVA\n",
    "6. Visualize diversity patterns\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "source": [
    "# Import necessary libraries\n",
    "import os\n",
    "import sys\n",
    "import pandas as pd\n",
    "import numpy as np\n",
    "import matplotlib.pyplot as plt\n",
    "import seaborn as sns\n",
    "import yaml\n",
    "from pathlib import Path\n",
    "\n",
    "# Set plotting style\n",
    "sns.set(style=\"whitegrid\")\n",
    "plt.rcParams['figure.figsize'] = (12, 8)\n",
    "\n",
    "# Add project root to Python path\n",
    "project_root = Path().resolve().parents[0]\n",
    "sys.path.append(str(project_root))\n",
    "\n",
    "# Import functions from metaphlan_tools\n",
    "from metaphlan_tools import (\n",
    "    load_metadata,\n",
    "    calculate_alpha_diversity,\n",
    "    compare_alpha_diversity,\n",
    "    calculate_beta_diversity,\n",
    "    perform_permanova,\n",
    "    plot_alpha_diversity_boxplot\n",
    ")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "source": [
    "# Load configuration\n",
    "config_path = project_root / 'config' / 'analysis_parameters.yml'\n",
    "with open(config_path, 'r') as f:\n",
    "    config = yaml.safe_load(f)\n",
    "\n",
    "# Set up paths\n",
    "processed_data_dir = project_root / 'data' / 'processed'\n",
    "results_dir = project_root / 'results'\n",
    "figures_dir = results_dir / 'figures'\n",
    "tables_dir = results_dir / 'tables'\n",
    "\n",
    "# Create directories if they don't exist\n",
    "figures_dir.mkdir(exist_ok=True, parents=True)\n",
    "tables_dir.mkdir(exist_ok=True, parents=True)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## 1. Load Data"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "source": [
    "# Load the combined abundance table\n",
    "abundance_file = processed_data_dir / 'combined_abundance.csv'\n",
    "abundance_df = pd.read_csv(abundance_file, index_col=0)\n",
    "\n",
    "# Load metadata\n",
    "metadata_file = project_root / config['metadata']['filename']\n",
    "metadata_df = load_metadata(metadata_file, config['metadata']['sample_id_column'])\n",
    "\n",
    "# Check data dimensions\n",
    "print(f\"Abundance data: {abundance_df.shape[0]} species, {abundance_df.shape[1]} samples\")\n",
    "print(f\"Metadata: {metadata_df.shape[0]} samples, {metadata_df.shape[1]} variables\")\n",
    "\n",
    "# Check sample overlap\n",
    "common_samples = set(abundance_df.columns).intersection(set(metadata_df.index))\n",
    "print(f\"Samples with both abundance and metadata: {len(common_samples)}\")\n",
    "\n",
    "# Filter to common samples if needed\n",
    "if len(common_samples) < len(abundance_df.columns):\n",
    "    abundance_df = abundance_df[list(common_samples)]"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## 2. Calculate Alpha Diversity"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "source": [
    "# Calculate alpha diversity metrics\n",
    "alpha_metrics = config['diversity']['alpha_metrics']\n",
    "alpha_df = calculate_alpha_diversity(abundance_df, metrics=alpha_metrics)\n",
    "\n",
    "# Display the first few rows\n",
    "alpha_df.head()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "source": [
    "# Save alpha diversity results\n",
    "alpha_file = tables_dir / 'alpha_diversity.csv'\n",
    "alpha_df.to_csv(alpha_file)\n",
    "print(f\"Alpha diversity saved to {alpha_file}\")"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## 3. Compare Alpha Diversity Between Clinical Groups"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "source": [
    "# Get group variables from config\n",
    "group_vars = config['metadata']['group_variables']\n",
    "\n",
    "# Test for differences in alpha diversity between groups\n",
    "for var in group_vars:\n",
    "    if var in metadata_df.columns:\n",
    "        print(f\"\\nAnalyzing differences in alpha diversity by {var}:\")\n",
    "        results = compare_alpha_diversity(alpha_df, metadata_df, var)\n",
    "        \n",
    "        # Print results\n",
    "        for metric, stats in results.items():\n",
    "            print(f\"  {metric}: {stats['test']} p-value = {stats['p-value']:.4f}\")\n",
    "        \n",
    "        # Create and save boxplot\n",
    "        fig = plot_alpha_diversity_boxplot(alpha_df, metadata_df, var)\n",
    "        boxplot_file = figures_dir / f\"alpha_diversity_{var}.png\"\n",
    "        fig.savefig(boxplot_file, dpi=config['visualization']['figure_dpi'], bbox_inches='tight')\n",
    "        plt.close(fig)\n",
    "        print(f\"  Boxplot saved to {boxplot_file}\")\n",
    "    else:\n",
    "        print(f\"Warning: Variable '{var}' not found in metadata\")"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## 4. Calculate Beta Diversity"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "source": [
    "# Calculate beta diversity\n",
    "beta_metric = config['diversity']['beta_metric']\n",
    "beta_dm = calculate_beta_diversity(abundance_df, metric=beta_metric)\n",
    "\n",
    "# Show a sample of the distance matrix\n",
    "print(f\"Beta diversity distance matrix ({beta_metric}):\")\n",
    "print(beta_dm.data[0:5, 0:5])\n",
    "print(f\"Size: {beta_dm.shape}\")"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## 5. Test for Community Differences with PERMANOVA"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "source": [
    "# Perform PERMANOVA tests\n",
    "permanova_results = {}\n",
    "\n",
    "for var in group_vars:\n",
    "    if var in metadata_df.columns:\n",
    "        print(f\"\\nPerforming PERMANOVA for {var}:\")\n",
    "        result = perform_permanova(beta_dm, metadata_df, var)\n",
    "        permanova_results[var] = result\n",
    "        \n",
    "        # Print results\n",
    "        print(f\"  Test statistic: {result['test-statistic']:.4f}\")\n",
    "        print(f\"  p-value: {result['p-value']:.4f}\")\n",
    "        print(f\"  Sample size: {result['sample size']}\")\n",
    "    else:\n",
    "        print(f\"Warning: Variable '{var}' not found in metadata\")\n",
    "\n",
    "# Save PERMANOVA results\n",
    "permanova_file = tables_dir / 'permanova_results.csv'\n",
    "permanova_df = pd.DataFrame(permanova_results).T\n",
    "permanova_df.to_csv(permanova_file)\n",
    "print(f\"\\nPERMANOVA results saved to {permanova_file}\")"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## 6. Visualize Beta Diversity with Ordination"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "source": [
    "# Import ordination function\n",
    "from metaphlan_tools.stats import plot_beta_diversity_ordination\n",
    "\n",
    "# Create ordination plots for each group variable\n",
    "for var in group_vars:\n",
    "    if var in metadata_df.columns:\n",
    "        print(f\"Creating ordination plot for {var}...\")\n",
    "        fig = plot_beta_diversity_ordination(beta_dm, metadata_df, var, method='PCoA')\n",
    "        \n",
    "        # Save the figure\n",
    "        ordination_file = figures_dir / f\"beta_diversity_pcoa_{var}.png\"\n",
    "        fig.savefig(ordination_file, dpi=config['visualization']['figure_dpi'], bbox_inches='tight')\n",
    "        plt.close(fig)\n",
    "        print(f\"Ordination plot saved to {ordination_file}\")"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Summary of Findings\n",
    "\n",
    "**Alpha Diversity:**\n",
    "- [Add your observations about alpha diversity results here]\n",
    "- Were there significant differences between clinical groups?\n",
    "- Which diversity metrics showed the strongest patterns?\n",
    "\n",
    "**Beta Diversity:**\n",
    "- [Add your observations about beta diversity and community composition here]\n",
    "- Did PERMANOVA indicate significant differences between groups?\n",
    "- What patterns are visible in the ordination plots?\n",
    "\n",
    "**Next Steps:**\n",
    "- Identify which specific species are driving these community differences\n",
    "- Analyze longitudinal changes in diversity metrics\n",
    "- Correlate diversity with clinical variables"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.8.5"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 4
}
