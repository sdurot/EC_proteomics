# EC Proteomics Analysis Pipeline

A comprehensive, modern data science pipeline for analyzing endothelial cell proteomics data with advanced statistical methods, interactive visualizations, and automated reporting.

![Python](https://img.shields.io/badge/python-3.8+-blue.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)
![Status](https://img.shields.io/badge/status-stable-brightgreen.svg)

## 📖 Overview

This repository provides a complete data science framework for LC-MS/MS proteomics analysis, specifically designed for endothelial cell biology research. The pipeline transforms raw proteomics data into publication-ready insights through modern statistical methods and interactive visualizations.

### 🔬 Publication

This code supports the analysis presented in:
> **"Multi-omics analysis of endothelial cells reveals the metabolic diversity that underlies endothelial cell functions"**  
> *bioRxiv*: https://doi.org/10.1101/2025.03.03.641143

## ✨ Features

### 🚀 Modern Data Science Approach
- **Object-oriented design** with the `ProteomicsAnalyzer` class
- **Method chaining** for streamlined analysis workflows
- **Type hints** and comprehensive documentation
- **Logging** for transparent analysis tracking

### 📊 Statistical Analysis
- **Principal Component Analysis (PCA)** with explained variance plots
- **Partial Least Squares Discriminant Analysis (PLS-DA)** with permutation testing
- **Differential expression analysis** with multiple testing correction
- **Quality control metrics** and diagnostic plots

### 🎨 Advanced Visualizations
- **Interactive plots** using Plotly for web-based exploration
- **Publication-quality** static plots with matplotlib/seaborn
- **Volcano plots** for differential expression results
- **Correlation heatmaps** and clustering dendrograms

### 📈 Automated Reporting
- **HTML reports** with comprehensive analysis summaries
- **Export capabilities** for further analysis in other tools
- **Protein lists** ready for pathway enrichment analysis

## 🛠️ Installation

### Option 1: Direct Installation (Recommended)

```bash
# Clone the repository
git clone https://github.com/sdurot/EC_proteomics.git
cd EC_proteomics

# Install dependencies
pip install -r requirements.txt

# Optional: Install in development mode
pip install -e .
```

### Option 2: Conda Environment

```bash
# Create conda environment
conda create -n ec-proteomics python=3.9
conda activate ec-proteomics

# Install dependencies
pip install -r requirements.txt
```

### Option 3: Docker (Coming Soon)

```bash
# Build Docker image
docker build -t ec-proteomics .

# Run analysis
docker run -v /path/to/data:/data ec-proteomics
```

## 📋 Requirements

### Core Dependencies
- **Python** ≥ 3.8
- **NumPy** ≥ 1.21.0 - Numerical computing
- **Pandas** ≥ 1.3.0 - Data manipulation
- **SciPy** ≥ 1.7.0 - Statistical functions
- **Scikit-learn** ≥ 1.0.0 - Machine learning algorithms

### Visualization
- **Matplotlib** ≥ 3.4.0 - Static plotting
- **Seaborn** ≥ 0.11.0 - Statistical visualizations
- **Plotly** ≥ 5.0.0 - Interactive plots

### Statistical Analysis
- **Statsmodels** ≥ 0.12.0 - Advanced statistics
- **Multipy** ≥ 0.15 - Multiple testing correction

### File Handling
- **OpenPyXL** ≥ 3.0.0 - Excel file support
- **XlsxWriter** - Excel export functionality

## 🚀 Quick Start

### 1. Basic Analysis

```python
from proteomics_analysis import ProteomicsAnalyzer

# Initialize analyzer
analyzer = ProteomicsAnalyzer(data_path=".", output_dir="results")

# Complete analysis pipeline
analyzer.load_data("Suppl_table_1.xlsx") \
        .preprocess_data() \
        .normalize_data() \
        .quality_control() \
        .perform_pca() \
        .perform_plsda() \
        .generate_report()
```

### 2. Jupyter Notebook

Open `proteomics_analysis_example.ipynb` for an interactive tutorial with detailed explanations and visualizations.

### 3. Command Line

```bash
# Run complete analysis
python proteomics_analysis.py

# Or if installed as package
ec-proteomics
```

## 📁 Data Structure

### Input Data Format

Your Excel file should contain:

```
Sheet 1: 'Px_LFQ_data'
- Column 1: Gene names
- Columns 2-N: Sample intensities (format: CELLLINE_TIMEPOINT_REPLICATE)
- Last columns: Protein metadata

Sheet 2: 'protein_genes_data'  
- Gene names and protein information
```

### Example File Structure
```
EC_proteomics/
├── proteomics_analysis.py          # Main analysis class
├── proteomics_analysis_example.ipynb # Jupyter tutorial
├── Suppl_table_1.xlsx              # Your data file
├── requirements.txt                 # Dependencies
├── setup.py                        # Package installation
├── README.md                       # This file
└── results/                        # Output directory
    ├── analysis_report.html         # Comprehensive report
    ├── quality_control.png          # QC plots
    ├── pca_analysis.html            # Interactive PCA
    ├── plsda_cell_type.png         # PLS-DA results
    ├── volcano_plot_*.html         # Differential analysis
    └── *.csv                       # Data tables
```

## 🔬 Analysis Workflow

### 1. Data Loading & Preprocessing
- Load Excel data with proper error handling
- Remove proteins with excessive missing values
- Generate metadata from sample names
- Quality control assessment

### 2. Normalization & Quality Control
- Z-score normalization (or min-max, robust scaling)
- Missing value analysis
- Sample correlation assessment
- Distribution plots and diagnostics

### 3. Dimensionality Reduction
- **PCA**: Identify major sources of variation
- **t-SNE/UMAP**: Non-linear dimensionality reduction (optional)
- Explained variance analysis
- Loading plots for feature interpretation

### 4. Statistical Classification
- **PLS-DA**: Supervised classification (BEC vs LEC)
- Permutation testing for significance
- Cross-validation for model validation
- Feature importance ranking

### 5. Differential Expression
- Statistical testing (t-test, Mann-Whitney U)
- Multiple testing correction (q-values, FDR)
- Effect size calculation (log2 fold change)
- Volcano plots and MA plots

### 6. Reporting & Export
- Automated HTML report generation
- Interactive visualizations
- Publication-ready figures
- Data export for downstream analysis

## 📊 Analysis Examples

### Cell Type Classification

```python
# Perform PLS-DA to identify BEC vs LEC markers
plsda_results = analyzer.perform_plsda(
    target_variable='cell_type',
    n_permutations=1000
)

# Get significant proteins
lec_markers = plsda_results['significant_proteins']['high_weights']
bec_markers = plsda_results['significant_proteins']['low_weights']

print(f"Found {len(lec_markers)} LEC-specific proteins")
print(f"Found {len(bec_markers)} BEC-specific proteins")
```

### Time-Course Analysis

```python
# Compare D5 vs D2 for each cell line
for cell_line in ['HDBEC', 'HDLEC', 'HUVEC', 'iLEC']:
    diff_results = analyzer.differential_analysis(
        reference_condition=f"{cell_line}_D2",
        comparison_conditions=[f"{cell_line}_D5"]
    )
```

### Custom Protein Visualization

```python
# Plot specific protein across conditions
def plot_protein(protein_name):
    # Implementation in the example notebook
    pass

plot_protein('PECAM1')  # Endothelial marker
```

## 🔧 Advanced Usage

### Custom Analysis Parameters

```python
# Initialize with custom settings
analyzer = ProteomicsAnalyzer(
    data_path="/path/to/data",
    output_dir="custom_results"
)

# Use different normalization
analyzer.normalize_data(method='robust')  # or 'minmax'

# Adjust PCA components
pca_results = analyzer.perform_pca(n_components=15)

# Custom differential analysis
diff_results = analyzer.differential_analysis(
    reference_condition="HUVEC_D2",
    comparison_conditions=["HUVEC_D5", "HUVEC_D7"],
    statistical_test='mannwhitney',
    multiple_testing_correction='fdr_bh'
)
```

### Integration with Other Tools

```python
# Export for pathway analysis
protein_lists = analyzer.results['plsda']['significant_proteins']

# Save for GSEA, STRING, etc.
with open('lec_proteins.txt', 'w') as f:
    f.write('\n'.join(protein_lists['high_weights'].index))
```

## 🎯 Key Improvements Over Original Code

### 🏗️ Code Organization
- **Object-oriented design** vs. procedural scripting
- **Modular functions** for better maintainability
- **Error handling** and input validation
- **Comprehensive logging** and progress tracking

### 📈 Statistical Enhancements
- **Permutation testing** for PLS-DA significance
- **Cross-validation** for model validation
- **Multiple testing corrections** (q-values, FDR)
- **Effect size calculations** and confidence intervals

### 🎨 Visualization Improvements
- **Interactive plots** with Plotly
- **Publication-quality** static figures
- **Consistent styling** and color schemes
- **Automated plot saving** in multiple formats

### 📊 Analysis Features
- **Quality control metrics** and diagnostics
- **Automated report generation** 
- **Batch effect detection** and correction options
- **Missing value analysis** and imputation

### 🔧 Usability Enhancements
- **Jupyter notebook tutorial** with examples
- **Command-line interface** for automation
- **Comprehensive documentation** and docstrings
- **Flexible input/output formats**

## 🤝 Contributing

We welcome contributions! Please see our contributing guidelines:

1. **Fork** the repository
2. **Create** a feature branch (`git checkout -b feature/amazing-feature`)
3. **Commit** your changes (`git commit -m 'Add amazing feature'`)
4. **Push** to the branch (`git push origin feature/amazing-feature`)
5. **Open** a Pull Request

### Development Setup

```bash
# Clone for development
git clone https://github.com/sdurot/EC_proteomics.git
cd EC_proteomics

# Install in development mode
pip install -e ".[dev]"

# Run tests
pytest tests/

# Format code
black proteomics_analysis.py
flake8 proteomics_analysis.py
```

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📚 Citation

If you use this code in your research, please cite:

```bibtex
@article{durot2025multiomics,
  title={Multi-omics analysis of endothelial cells reveals the metabolic diversity that underlies endothelial cell functions},
  author={Durot, S. and colleagues},
  journal={bioRxiv},
  year={2025},
  doi={10.1101/2025.03.03.641143}
}
```

## 🆘 Support

- **Issues**: [GitHub Issues](https://github.com/sdurot/EC_proteomics/issues)
- **Discussions**: [GitHub Discussions](https://github.com/sdurot/EC_proteomics/discussions)
- **Email**: [Contact the authors](mailto:your-email@example.com)

## 🗺️ Roadmap

### Upcoming Features
- [ ] **Docker containerization** for reproducible environments
- [ ] **Pathway enrichment analysis** integration
- [ ] **Network analysis** capabilities (heat diffusion, RWR)
- [ ] **Multi-omics integration** tools
- [ ] **Batch effect correction** methods
- [ ] **Machine learning classifiers** (Random Forest, SVM)
- [ ] **Time-series analysis** for longitudinal data
- [ ] **GUI interface** for non-programmers

### Long-term Goals
- [ ] **Web application** for online analysis
- [ ] **Database integration** (UniProt, STRING)
- [ ] **Automated report generation** with AI insights
- [ ] **Real-time collaboration** features

## 🏆 Acknowledgments

- Original analysis concepts from the bioRxiv publication
- Inspiration from modern data science practices
- Community feedback and contributions
- Open-source libraries that make this possible

---

**Happy analyzing!** 🧬✨

For questions, issues, or contributions, please don't hesitate to reach out. We're excited to see how you use this pipeline for your proteomics research!