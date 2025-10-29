"""
Configuration file for EC Proteomics Analysis Pipeline

This file contains default settings and parameters for the analysis pipeline.
You can modify these values to customize the analysis for your specific needs.
"""

# Data Processing Configuration
DATA_CONFIG = {
    "excel_file": "Suppl_table_1.xlsx",
    "intensity_sheet": "Px_LFQ_data",
    "protein_info_sheet": "protein_genes_data",
    "gene_name_column": "Gene names",
    "remove_outliers": True,
    "min_valid_values": 0.7,  # Minimum fraction of valid values per protein
}

# Normalization Configuration
NORMALIZATION_CONFIG = {
    "method": "zscore",  # Options: "zscore", "minmax", "robust"
    "center": True,
    "scale": True,
}

# PCA Configuration
PCA_CONFIG = {
    "n_components": 10,
    "center_data": True,
    "scale_data": False,  # Already done in normalization
}

# PLS-DA Configuration
PLSDA_CONFIG = {
    "n_components": 2,
    "scale": True,
    "n_permutations": 1000,
    "significance_level": 0.05,
    "cross_validation_folds": 5,
}

# Differential Analysis Configuration
DIFF_ANALYSIS_CONFIG = {
    "statistical_test": "ttest",  # Options: "ttest", "mannwhitney"
    "multiple_testing_correction": "qvalue",  # Options: "qvalue", "fdr_bh", "bonferroni"
    "log2fc_threshold": 0.5,
    "pvalue_threshold": 0.05,
    "qvalue_threshold": 0.05,
}

# Plotting Configuration
PLOT_CONFIG = {
    "figure_dpi": 300,
    "figure_format": ["png", "svg"],
    "color_palette": "husl",
    "font_size": 12,
    "font_family": "Arial",
    "interactive_plots": True,
    "save_plots": True,
}

# Output Configuration
OUTPUT_CONFIG = {
    "output_directory": "results",
    "save_intermediate_results": True,
    "generate_html_report": True,
    "export_protein_lists": True,
    "log_level": "INFO",  # Options: "DEBUG", "INFO", "WARNING", "ERROR"
}

# Cell Type Configuration
CELL_TYPE_CONFIG = {
    "bec_cell_lines": ["HDBEC", "HUVEC"],
    "lec_cell_lines": ["HDLEC", "iLEC"],
    "time_points": ["D2", "D5", "D7"],
    "sample_name_separator": "_",
}

# Quality Control Configuration
QC_CONFIG = {
    "correlation_method": "spearman",
    "missing_value_threshold": 0.3,  # Maximum fraction of missing values allowed
    "outlier_detection_method": "iqr",  # Options: "iqr", "zscore", "isolation_forest"
    "outlier_threshold": 3.0,
}

# Advanced Analysis Configuration
ADVANCED_CONFIG = {
    "pathway_analysis": {
        "enabled": False,
        "database": "KEGG",  # Options: "KEGG", "GO", "Reactome"
        "organism": "hsa",  # Human
    },
    "network_analysis": {
        "enabled": False,
        "string_db_version": "11.5",
        "confidence_threshold": 0.4,
    },
    "batch_correction": {
        "enabled": False,
        "method": "combat",  # Options: "combat", "limma"
    },
}

# File paths (will be set dynamically)
PATHS = {
    "data_path": ".",
    "output_path": "results",
    "log_file": "analysis.log",
}

# Sample metadata parsing configuration
METADATA_CONFIG = {
    "sample_column_format": "{cell_line}_{time_point}_{replicate}",
    "intensity_column_start": 2,  # 0-based index
    "intensity_column_end": -2,   # Negative index from end
}

# Reproducibility settings
REPRODUCIBILITY_CONFIG = {
    "random_seed": 42,
    "set_numpy_seed": True,
    "set_sklearn_seed": True,
}

# Performance settings
PERFORMANCE_CONFIG = {
    "n_jobs": -1,  # Number of parallel jobs (-1 for all available cores)
    "memory_limit": "8GB",  # Maximum memory usage
    "chunk_size": 1000,  # For processing large datasets in chunks
}

def get_config():
    """Return the complete configuration dictionary."""
    return {
        "data": DATA_CONFIG,
        "normalization": NORMALIZATION_CONFIG,
        "pca": PCA_CONFIG,
        "plsda": PLSDA_CONFIG,
        "differential": DIFF_ANALYSIS_CONFIG,
        "plotting": PLOT_CONFIG,
        "output": OUTPUT_CONFIG,
        "cell_types": CELL_TYPE_CONFIG,
        "quality_control": QC_CONFIG,
        "advanced": ADVANCED_CONFIG,
        "paths": PATHS,
        "metadata": METADATA_CONFIG,
        "reproducibility": REPRODUCIBILITY_CONFIG,
        "performance": PERFORMANCE_CONFIG,
    }

def update_config(config_updates):
    """
    Update configuration with custom values.
    
    Parameters:
    -----------
    config_updates : dict
        Dictionary with configuration updates
    """
    config = get_config()
    
    def update_nested_dict(original, updates):
        for key, value in updates.items():
            if isinstance(value, dict) and key in original:
                update_nested_dict(original[key], value)
            else:
                original[key] = value
    
    update_nested_dict(config, config_updates)
    return config