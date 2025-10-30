# EC Proteomics Analysis Pipeline

## Overview

This code supports the analysis presented in:
> **"Multi-omics analysis of endothelial cells reveals the metabolic diversity that underlies endothelial cell functions"**  
> *bioRxiv*: https://doi.org/10.1101/2025.03.03.641143

![Python](https://img.shields.io/badge/python-3.8+-blue.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)
![Status](https://img.shields.io/badge/status-stable-brightgreen.svg)

## Features

- **Object-oriented design** with the `ProteomicsAnalyzer` class
- **Principal Component Analysis (PCA)** with explained variance plots
- **Partial Least Squares Discriminant Analysis (PLS-DA)** with permutation testing
- **Interactive plots** using Plotly and publication-quality static plots
- **Differential expression analysis** with multiple testing correction
- **Quality control metrics** and diagnostic plots
- **Automated HTML reports** and data export capabilities

## Installation

### Automated Setup (Recommended)

```bash
git clone https://github.com/sdurot/EC_proteomics.git
cd EC_proteomics

# Setup virtual environment and launch Jupyter. UV is recommended.
python setup_environment.py --method uv --jupyter
```

### Manual Virtual Environment

```bash
python -m venv .venv
.venv\Scripts\activate  # Windows
source .venv/bin/activate  # Unix/Linux/macOS
pip install -r requirements.txt
```

### Alternative Methods

- **venv**: `python setup_environment.py --method venv --jupyter`
- **Conda**: `python setup_environment.py --method conda --jupyter`
- **Global Python**: If packages already available, run directly

## Requirements

- **Python** ≥ 3.8
- **Core**: NumPy, Pandas, SciPy, Scikit-learn
- **Visualization**: Matplotlib, Seaborn, Plotly
- **Statistics**: Statsmodels, Multipy
- **Files**: OpenPyXL, XlsxWriter

## Quick Start

### Basic Analysis

```python
from proteomics_analysis import ProteomicsAnalyzer

# Initialize and run analysis
analyzer = ProteomicsAnalyzer()
analyzer.run_analysis()

# Or step by step
analyzer.perform_pca(time_points=['D2', 'D5'])
analyzer.perform_plsda(threshold=0.9)
```

### Jupyter Notebook

Open `EC_proteomics_analysis.ipynb` for an interactive tutorial with all analysis steps.

### Command Line

```bash
python proteomics_analysis.py
```

## Data Structure

### Input Format
Excel file with protein intensities and sample names in format: `CELLLINE_TIMEPOINT_REPLICATE`

### File Structure
```
EC_proteomics/
├── proteomics_analysis.py          # Main analysis class
├── EC_proteomics_analysis.ipynb    # Jupyter tutorial
├── setup_environment.py            # Environment setup
├── Suppl_table_1_Px_data.xlsx      # Proteomics data
├── requirements.txt                 # Dependencies
└── results/                        # Output directory
    ├── *_associated_proteins_90pct.csv
    └── *.png
```

## Analysis Workflow

1. **Data Loading**: Excel data with error handling and quality control
2. **Normalization**: Z-score normalization and missing value analysis  
3. **PCA**: Dimensionality reduction with explained variance analysis
4. **PLS-DA**: Supervised classification with permutation testing
5. **Differential Expression**: Statistical testing with multiple testing correction
6. **Reporting**: HTML reports, interactive plots, and data export

## Usage Examples

```python
# Cell type classification
plsda_results = analyzer.perform_plsda(threshold=0.9)

# PCA with specific time points
pca_results = analyzer.perform_pca(time_points=['D2', 'D5'])

# Time-course analysis
for cell_line in ['HDBEC', 'iLEC']:
    diff_results = analyzer.differential_analysis(
        reference_condition=f"{cell_line}_D2",
        comparison_conditions=[f"{cell_line}_D5"]
    )
```

## Advanced Usage

```python
# Custom parameters
analyzer = ProteomicsAnalyzer(output_dir="custom_results")
pca_results = analyzer.perform_pca(n_components=15)

# Export for pathway analysis
protein_lists = analyzer.results['plsda']['significant_proteins']
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.


## Support

- Issues: [GitHub Issues](https://github.com/sdurot/EC_proteomics/issues)