#!/usr/bin/env python3
"""
Endothelial Cell Proteomics Analysis Pipeline

A comprehensive data science pipeline for analyzing proteomics data from endothelial cells.
This module provides tools for data preprocessing, normalization, statistical analysis,
and visualization of LC-MS/MS proteomics datasets.

Authors: sdurot
Created: 2025
Last Modified: October 2025

Publication: "Multi-omics analysis of endothelial cells reveals the metabolic diversity 
that underlies endothelial cell functions" (bioRxiv: https://doi.org/10.1101/2025.03.03.641143)
"""

import os
import warnings
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Union
import logging

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.cross_decomposition import PLSRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score
from multipy.fdr import qvalue
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning)

class ProteomicsAnalyzer:
    """
    A comprehensive class for endothelial cell proteomics data analysis.
    
    This class provides methods for loading, preprocessing, normalizing, and analyzing
    LC-MS/MS proteomics data with a focus on endothelial cell biology.
    """
    
    def __init__(self, data_path: str = None, output_dir: str = "results"):
        """
        Initialize the ProteomicsAnalyzer.
        
        Parameters:
        -----------
        data_path : str, optional
            Path to the directory containing proteomics data files
        output_dir : str, default "results"
            Directory to save analysis results and plots
        """
        self.data_path = Path(data_path) if data_path else Path.cwd()
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # Initialize data containers
        self.raw_data = None
        self.protein_info = None
        self.clean_data = None
        self.normalized_data = None
        self.metadata = None
        self.results = {}
        
        # Configure matplotlib for publication-quality figures
        self._setup_plotting()
        
        logger.info(f"ProteomicsAnalyzer initialized with data path: {self.data_path}")
        logger.info(f"Results will be saved to: {self.output_dir}")
        
        # Try to load data automatically if the default file exists
        self._try_auto_load_data()
    
    @property
    def data(self):
        """
        Get the most processed version of the data available.
        
        Returns:
        --------
        pd.DataFrame
            The normalized data if available, otherwise clean data, otherwise raw data
        """
        if self.normalized_data is not None:
            return self.normalized_data
        elif self.clean_data is not None:
            return self.clean_data
        elif self.raw_data is not None:
            return self.raw_data
        else:
            logger.warning("No data loaded. Please use load_data() method first.")
            return None
    
    def _try_auto_load_data(self):
        """Try to automatically load data if the default file exists."""
        default_files = ["Suppl_table_1_Px_data.xlsx", "Suppl_table_1.xlsx"]
        
        for filename in default_files:
            file_path = self.data_path / filename
            if file_path.exists():
                try:
                    logger.info(f"Found data file: {filename}. Loading automatically...")
                    self.load_data(filename)
                    return
                except Exception as e:
                    logger.debug(f"Failed to auto-load {filename}: {e}")
                    continue
        
        logger.info("No default data file found. Use load_data() method to load your data.")
    
    def _setup_plotting(self):
        """Configure matplotlib and seaborn for publication-quality plots."""
        plt.style.use('default')
        plt.rcParams.update({
            'font.size': 12,
            'font.family': 'Arial',
            'figure.dpi': 300,
            'savefig.dpi': 300,
            'savefig.bbox': 'tight',
            'axes.spines.top': False,
            'axes.spines.right': False,
            'axes.grid': True,
            'grid.alpha': 0.3
        })
        sns.set_palette("husl")
    
    def load_data(self, excel_file: str = "Suppl_table_1.xlsx") -> 'ProteomicsAnalyzer':
        """
        Load proteomics data from Excel file.
        
        Parameters:
        -----------
        excel_file : str
            Name of the Excel file containing proteomics data
            
        Returns:
        --------
        self : ProteomicsAnalyzer
            Returns self for method chaining
        """
        try:
            file_path = self.data_path / excel_file
            logger.info(f"Loading data from {file_path}")
            
            # Load main proteomics data
            self.raw_data = pd.read_excel(
                file_path, 
                sheet_name='Px_LFQ_data'
            ).rename(columns={'newColName': 'Gene names'})
            
            # Load protein information
            self.protein_info = pd.read_excel(
                file_path, 
                sheet_name='protein_genes_data'
            ).rename(columns={'Genes': 'Gene names'})
            
            # Merge datasets
            self.raw_data = self.raw_data.merge(
                self.protein_info, 
                on='Gene names'
            ).rename(columns={'First.Protein:Description': 'Protein description'})
            
            logger.info(f"Loaded {self.raw_data.shape[0]} proteins across {self.raw_data.shape[1]-3} samples")
            logger.debug(f"Column names: {list(self.raw_data.columns)}")
            
            # Create genes_prot_data for compatibility with original code
            self._create_protein_metadata()
            
        except FileNotFoundError:
            logger.error(f"File {excel_file} not found in {self.data_path}")
            raise
        except Exception as e:
            logger.error(f"Error loading data: {e}")
            raise
            
        return self
    
    def _create_protein_metadata(self):
        """Create protein metadata DataFrame for compatibility."""
        # Find available columns
        desc_cols = ['Protein description', 'First.Protein.Description', 'First.Protein:Description', 'Protein.Description']
        id_cols = ['Protein.Ids', 'Protein IDs', 'Protein_IDs', 'ProteinIDs']
        
        desc_col = None
        id_col = None
        
        for col in desc_cols:
            if col in self.raw_data.columns:
                desc_col = col
                break
        
        for col in id_cols:
            if col in self.raw_data.columns:
                id_col = col
                break
        
        # Create metadata DataFrame
        metadata_dict = {'Gene names': self.raw_data['Gene names']}
        
        if desc_col:
            metadata_dict['Protein description'] = self.raw_data[desc_col]
        
        if id_col:
            metadata_dict['Protein IDs'] = self.raw_data[id_col]
        
        self.genes_prot_data = pd.DataFrame(metadata_dict).reset_index(drop=True)
    
    def preprocess_data(self, remove_outliers: bool = True) -> 'ProteomicsAnalyzer':
        """
        Preprocess and clean the proteomics data.
        
        Parameters:
        -----------
        remove_outliers : bool, default True
            Whether to remove identified outlier samples
            
        Returns:
        --------
        self : ProteomicsAnalyzer
            Returns self for method chaining
        """
        logger.info("Preprocessing proteomics data...")
        
        # Create clean dataset (remove proteins with missing values)
        self.clean_data = self.raw_data.dropna().copy()
        self.clean_data.index = self.clean_data.iloc[:, 1]  # Set gene names as index
        
        # Generate metadata from column names
        intensity_columns = self.clean_data.columns[2:-2]  # Exclude metadata columns
        metadata_split = pd.DataFrame([col.split('_') for col in intensity_columns])
        
        self.metadata = pd.DataFrame({
            'sample_id': intensity_columns,
            'cell_line': metadata_split.iloc[:, 0],
            'time_point': metadata_split.iloc[:, 1],
            'replicate': metadata_split.iloc[:, 2] if metadata_split.shape[1] > 2 else 1,
            'condition': metadata_split.iloc[:, 0] + "_" + metadata_split.iloc[:, 1]
        })
        
        # Add cell type classification
        self.metadata['cell_type'] = self.metadata['cell_line'].apply(
            lambda x: 'BEC' if x in ['HDBEC', 'HUVEC'] else 'LEC'
        )
        
        logger.info(f"Clean dataset: {self.clean_data.shape[0]} proteins, {len(intensity_columns)} samples")
        logger.info(f"Cell lines: {self.metadata['cell_line'].unique()}")
        logger.info(f"Time points: {self.metadata['time_point'].unique()}")
        
        return self
    
    def normalize_data(self, method: str = 'zscore') -> 'ProteomicsAnalyzer':
        """
        Normalize proteomics data using specified method.
        
        Parameters:
        -----------
        method : str, default 'zscore'
            Normalization method ('zscore', 'minmax', 'robust')
            
        Returns:
        --------
        self : ProteomicsAnalyzer
            Returns self for method chaining
        """
        logger.info(f"Normalizing data using {method} method...")
        
        intensity_data = self.clean_data.iloc[:, 2:-2]
        
        if method == 'zscore':
            self.normalized_data = intensity_data.T.apply(stats.zscore).T
        elif method == 'minmax':
            from sklearn.preprocessing import MinMaxScaler
            scaler = MinMaxScaler()
            self.normalized_data = pd.DataFrame(
                scaler.fit_transform(intensity_data.T).T,
                index=intensity_data.index,
                columns=intensity_data.columns
            )
        elif method == 'robust':
            from sklearn.preprocessing import RobustScaler
            scaler = RobustScaler()
            self.normalized_data = pd.DataFrame(
                scaler.fit_transform(intensity_data.T).T,
                index=intensity_data.index,
                columns=intensity_data.columns
            )
        else:
            raise ValueError(f"Unknown normalization method: {method}")
        
        logger.info("Data normalization completed")
        return self
    
    def quality_control(self, save_plots: bool = True) -> Dict:
        """
        Perform quality control analysis and generate diagnostic plots.
        
        Parameters:
        -----------
        save_plots : bool, default True
            Whether to save plots to output directory
            
        Returns:
        --------
        dict : Quality control metrics
        """
        logger.info("Performing quality control analysis...")
        
        qc_results = {}
        
        # 1. Data completeness
        completeness = (1 - self.raw_data.iloc[:, 2:-2].isnull().sum() / len(self.raw_data)) * 100
        qc_results['completeness'] = completeness.describe()
        
        # 2. Distribution plots
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # Raw data distribution
        sample_data = self.raw_data.iloc[:, 2:8].values.flatten()  # First few samples
        sample_data = sample_data[~np.isnan(sample_data)]
        axes[0, 0].hist(np.log10(sample_data), bins=50, alpha=0.7)
        axes[0, 0].set_title('Raw Data Distribution (log10)')
        axes[0, 0].set_xlabel('log10(Intensity)')
        axes[0, 0].set_ylabel('Frequency')
        
        # Normalized data distribution
        norm_sample = self.normalized_data.iloc[:, :6].values.flatten()
        norm_sample = norm_sample[~np.isnan(norm_sample)]
        axes[0, 1].hist(norm_sample, bins=50, alpha=0.7)
        axes[0, 1].set_title('Normalized Data Distribution')
        axes[0, 1].set_xlabel('Z-score')
        axes[0, 1].set_ylabel('Frequency')
        
        # Sample correlation heatmap
        correlation_matrix = self.normalized_data.corr(method='spearman')
        sns.heatmap(correlation_matrix, ax=axes[1, 0], cmap='RdBu_r', 
                   center=0, square=True, cbar_kws={'shrink': 0.8})
        axes[1, 0].set_title('Sample Correlation Matrix')
        
        # Missing values heatmap
        missing_data = self.raw_data.iloc[:, 2:-2].isnull()
        sns.heatmap(missing_data.iloc[:100], ax=axes[1, 1], cmap='viridis',
                   cbar_kws={'shrink': 0.8})
        axes[1, 1].set_title('Missing Values Pattern (First 100 proteins)')
        
        plt.tight_layout()
        if save_plots:
            plt.savefig(self.output_dir / 'quality_control.png')
        plt.show()
        
        # 3. Sample-wise statistics
        qc_results['sample_stats'] = {
            'median_intensity': self.normalized_data.median(axis=0),
            'protein_detection': (self.raw_data.iloc[:, 2:-2] > 0).sum(axis=0)
        }
        
        logger.info("Quality control analysis completed")
        return qc_results
    
    def perform_pca(self, n_components: int = 10, save_plots: bool = True, 
                    time_points: List[str] = None) -> Dict:
        """
        Perform Principal Component Analysis.
        
        Parameters:
        -----------
        n_components : int, default 10
            Number of principal components to calculate
        save_plots : bool, default True
            Whether to save plots to output directory
        time_points : list, optional
            List of time points to include (e.g., ['D2', 'D5'])
            
        Returns:
        --------
        dict : PCA results including components, explained variance, and scores
        """
        # Ensure data is preprocessed and normalized
        if self.normalized_data is None:
            logger.info("Data not normalized. Preprocessing and normalizing data first...")
            self.preprocess_data()
            self.normalize_data()
        
        # Filter data by time points if specified
        if time_points:
            logger.info(f"Performing PCA with {n_components} components for time points: {time_points}")
            # Filter metadata for specified time points
            time_mask = self.metadata['time_point'].isin(time_points)
            filtered_metadata = self.metadata[time_mask]
            filtered_samples = filtered_metadata['sample_id']
            
            # Filter normalized data
            pca_data = self.normalized_data[filtered_samples].T  # Samples as rows, proteins as columns
            metadata_for_plot = filtered_metadata
        else:
            logger.info(f"Performing PCA with {n_components} components...")
            # Prepare data for PCA
            pca_data = self.normalized_data.T  # Samples as rows, proteins as columns
            metadata_for_plot = self.metadata
        
        # Perform PCA
        pca = PCA(n_components=n_components)
        pca_scores = pca.fit_transform(pca_data)
        
        # Create results dictionary
        pca_results = {
            'pca_model': pca,
            'scores': pd.DataFrame(
                pca_scores, 
                index=pca_data.index,  # Use the actual filtered data index
                columns=[f'PC{i+1}' for i in range(n_components)]
            ),
            'explained_variance_ratio': pca.explained_variance_ratio_,
            'cumulative_variance': np.cumsum(pca.explained_variance_ratio_)
        }
        
        # Create PCA plots
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=('Explained Variance', 'PCA Score Plot', 
                          'PC Loadings', 'Cumulative Variance'),
            specs=[[{"secondary_y": False}, {"secondary_y": False}],
                   [{"secondary_y": False}, {"secondary_y": False}]]
        )
        
        # Explained variance plot
        fig.add_trace(
            go.Bar(x=list(range(1, n_components+1)), 
                  y=pca.explained_variance_ratio_,
                  name='Explained Variance'),
            row=1, col=1
        )
        
        # PCA score plot with metadata
        pca_df = pca_results['scores'].copy()
        pca_df = pca_df.merge(metadata_for_plot.set_index('sample_id'), 
                             left_index=True, right_index=True)
        
        for cell_type in pca_df['cell_type'].unique():
            mask = pca_df['cell_type'] == cell_type
            fig.add_trace(
                go.Scatter(
                    x=pca_df.loc[mask, 'PC1'],
                    y=pca_df.loc[mask, 'PC2'],
                    mode='markers',
                    name=cell_type,
                    text=pca_df.loc[mask, 'condition']
                ),
                row=1, col=2
            )
        
        # PC loadings plot
        fig.add_trace(
            go.Scatter(
                x=pca.components_[0, :],
                y=pca.components_[1, :],
                mode='markers',
                name='Protein Loadings',
                opacity=0.6
            ),
            row=2, col=1
        )
        
        # Cumulative variance plot
        fig.add_trace(
            go.Scatter(
                x=list(range(1, n_components+1)),
                y=pca_results['cumulative_variance'],
                mode='lines+markers',
                name='Cumulative Variance'
            ),
            row=2, col=2
        )
        
        fig.update_layout(height=800, showlegend=True, 
                         title_text="Principal Component Analysis Results")
        fig.update_xaxes(title_text="Principal Component", row=1, col=1)
        fig.update_yaxes(title_text="Explained Variance Ratio", row=1, col=1)
        fig.update_xaxes(title_text=f"PC1 ({pca.explained_variance_ratio_[0]:.2%})", row=1, col=2)
        fig.update_yaxes(title_text=f"PC2 ({pca.explained_variance_ratio_[1]:.2%})", row=1, col=2)
        
        if save_plots:
            fig.write_html(self.output_dir / 'pca_analysis.html')
        fig.show()
        
        self.results['pca'] = pca_results
        logger.info(f"PCA completed. PC1 and PC2 explain {pca.explained_variance_ratio_[:2].sum():.2%} of variance")
        
        return pca_results
    
    def perform_plsda(self, target_variable: str = 'cell_type', n_components: int = 2, 
                      n_permutations: int = 1000, save_plots: bool = True, 
                      threshold_percentile: float = 90.0) -> Dict:
        """
        Perform Partial Least Squares Discriminant Analysis (PLS-DA).
        
        Parameters:
        -----------
        target_variable : str, default 'cell_type'
            Target variable for classification ('cell_type', 'cell_line', etc.)
        n_components : int, default 2
            Number of PLS components
        n_permutations : int, default 10000
            Number of permutations for significance testing
        save_plots : bool, default True
            Whether to save plots to output directory
        threshold_percentile : float, default 90.0
            Percentile threshold for significant weights
            
        Returns:
        --------
        dict : PLS-DA results including model, weights, and significance
        """
        # Ensure data is preprocessed and normalized
        if self.normalized_data is None:
            logger.info("Data not normalized. Preprocessing and normalizing data first...")
            self.preprocess_data()
            self.normalize_data()
        
        logger.info(f"Performing PLS-DA for {target_variable}...")
        
        # Prepare data using the same approach as original code
        # Use raw data with NaN replaced by 0, similar to original
        X = self.raw_data.iloc[:, 2:-2].replace(np.nan, 0).T  # Samples x Proteins
        
        # Encode cell type: BEC=0, LEC=1 (same as original)
        if target_variable == 'cell_type':
            celltype = np.where(self.metadata['cell_line'].isin(['HDBEC', 'HUVEC']), 0, 1)
        else:
            celltype = pd.get_dummies(self.metadata[target_variable]).iloc[:, 0].values
        
        # Fit PLS-DA model
        plsda = PLSRegression(n_components=n_components, scale=True)
        plsda.fit(X, celltype)
        
        # Get scores and weights (using x_scores_ and x_weights_ for consistency)
        scores = pd.DataFrame(plsda.x_scores_, 
                            index=X.index, 
                            columns=[f'LV{i+1}' for i in range(n_components)])
        weights = pd.DataFrame(plsda.x_weights_, 
                             index=self.raw_data.iloc[:, 1],  # Gene names as index
                             columns=[f'LV{i+1}' for i in range(n_components)])
        
        # Bootstrapping for significance testing (same as original)
        logger.info(f"Running {n_permutations} permutations for significance testing...")
        pls_boot_90_percentile = []
        pls_boot_10_percentile = []
        
        for _ in range(n_permutations):
            celltype_shuffle = np.random.permutation(celltype)
            plsda_perm = PLSRegression(n_components=n_components, scale=True)
            plsda_perm.fit(X, celltype_shuffle)
            weights_perm = pd.DataFrame(plsda_perm.x_weights_)[0]  # LV1 weights
            pls_boot_90_percentile.append(weights_perm.quantile(threshold_percentile/100))
            pls_boot_10_percentile.append(weights_perm.quantile((100-threshold_percentile)/100))
        
        # Calculate thresholds
        pls_top_threshold = np.mean(pls_boot_90_percentile)
        pls_bottom_threshold = np.mean(pls_boot_10_percentile)
        
        # Identify significant proteins
        lv1_weights = weights['LV1']
        bec_associated = lv1_weights[lv1_weights <= pls_bottom_threshold]
        lec_associated = lv1_weights[lv1_weights >= pls_top_threshold]
        
        # Additional bootstrapping to test significance of number of proteins
        pls_boot_number_top = []
        pls_boot_number_bottom = []
        pls_boot_number_combined = []
        
        logger.info("Running secondary bootstrap for protein count significance...")
        for _ in range(1000):
            celltype_shuffle = np.random.permutation(celltype)
            plsda_perm = PLSRegression(n_components=n_components, scale=True)
            plsda_perm.fit(X, celltype_shuffle)
            weights_lv1_perm = pd.DataFrame(plsda_perm.x_weights_).iloc[:, 0]
            
            number_top = (weights_lv1_perm >= pls_top_threshold).sum()
            number_bottom = (weights_lv1_perm <= pls_bottom_threshold).sum()
            
            pls_boot_number_top.append(number_top)
            pls_boot_number_bottom.append(number_bottom)
            pls_boot_number_combined.append(number_top + number_bottom)
        
        # Perform cross-validation (manual approach to avoid feature dimension mismatch)
        from sklearn.model_selection import StratifiedKFold
        from sklearn.metrics import roc_auc_score
        
        logger.info("Performing 5-fold cross-validation...")
        skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
        cv_scores = []
        
        for fold, (train_idx, test_idx) in enumerate(skf.split(X.values, celltype)):
            X_train, X_test = X.values[train_idx], X.values[test_idx]
            y_train, y_test = celltype[train_idx], celltype[test_idx]
            
            # Create and fit a new model for this fold
            fold_model = PLSRegression(n_components=n_components, scale=True)
            fold_model.fit(X_train, y_train)
            
            # Get decision scores
            y_scores = fold_model.predict(X_test).ravel()
            
            # Calculate ROC AUC if we have both classes
            if len(np.unique(y_test)) > 1:
                auc = roc_auc_score(y_test, y_scores)
                cv_scores.append(auc)
            else:
                # If only one class, assign perfect score
                cv_scores.append(1.0)
        
        cv_scores = np.array(cv_scores)
        
        # Calculate model performance metrics
        y_pred = plsda.predict(X)
        r2_score = plsda.score(X, celltype)
        
        plsda_results = {
            'model': plsda,
            'scores': scores,
            'weights': weights,
            'bec_associated': bec_associated,
            'lec_associated': lec_associated,
            'thresholds': {
                'top': pls_top_threshold,
                'bottom': pls_bottom_threshold
            },
            'bootstrap_counts': {
                'top': pls_boot_number_top,
                'bottom': pls_boot_number_bottom,
                'combined': pls_boot_number_combined
            },
            'cross_validation': {
                'cv_scores': cv_scores,
                'mean_auc': cv_scores.mean(),
                'std_auc': cv_scores.std()
            },
            'model_performance': {
                'r2_score': r2_score,
                'cv_score': cv_scores.mean()  # Use our cross-validation AUC as the CV score
            }
        }
        
        # Create visualization (4 subplots as in original)
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # Score plot
        scores_df = scores.copy()
        scores_df = scores_df.merge(self.metadata.set_index('sample_id'), 
                                  left_index=True, right_index=True)
        
        for group in scores_df[target_variable].unique():
            mask = scores_df[target_variable] == group
            axes[0, 0].scatter(scores_df.loc[mask, 'LV1'], scores_df.loc[mask, 'LV2'], 
                             label=group, s=60, alpha=0.7)
        
        axes[0, 0].set_xlabel('Component 1')
        axes[0, 0].set_ylabel('Component 2')
        axes[0, 0].set_title('PLS-DA Score Plot')
        axes[0, 0].legend()
        axes[0, 0].grid(True, alpha=0.3)
        
        # Ordered weights plot with thresholds
        sorted_weights = lv1_weights.sort_values()
        axes[0, 1].scatter(range(len(sorted_weights)), sorted_weights, 
                          marker='.', s=1, alpha=0.5, color='black')
        axes[0, 1].axhline(pls_top_threshold, color='red', linestyle='-', 
                          label=f'{threshold_percentile}% threshold')
        axes[0, 1].axhline(pls_bottom_threshold, color='red', linestyle='-')
        axes[0, 1].set_xlabel('Proteins')
        axes[0, 1].set_ylabel('Weights LV1')
        axes[0, 1].set_title('Ordered PLS-DA Weights with Thresholds')
        axes[0, 1].legend()
        axes[0, 1].grid(True, alpha=0.3)
        
        # Bootstrap distribution plot
        axes[1, 0].hist(pls_boot_number_combined, bins=30, alpha=0.7, 
                       color='lightblue', label='Total', density=True)
        axes[1, 0].hist(pls_boot_number_top, bins=30, alpha=0.7, 
                       color='orange', label='Top', density=True)
        axes[1, 0].hist(pls_boot_number_bottom, bins=30, alpha=0.7, 
                       color='green', label='Bottom', density=True)
        
        # Add observed values as vertical lines
        axes[1, 0].axvline(len(bec_associated) + len(lec_associated), 
                          color='blue', linestyle='--', linewidth=2, label='Observed Total')
        axes[1, 0].axvline(len(lec_associated), color='orange', linestyle='--', linewidth=2)
        axes[1, 0].axvline(len(bec_associated), color='green', linestyle='--', linewidth=2)
        
        axes[1, 0].set_xlabel('Number of proteins')
        axes[1, 0].set_ylabel('Density')
        axes[1, 0].set_title('Permutation Test: Protein Count Distribution')
        axes[1, 0].legend()
        axes[1, 0].grid(True, alpha=0.3)
        
        # Plot cross-validation results (use the cv_scores computed earlier)
        bars = axes[1, 1].bar(range(len(cv_scores)), cv_scores, alpha=0.7, color='skyblue', edgecolor='navy')
        axes[1, 1].axhline(cv_scores.mean(), color='red', linestyle='--', 
                          label=f'Mean CV Score: {cv_scores.mean():.3f}')
        
        # Add value labels on bars
        for i, (bar, score) in enumerate(zip(bars, cv_scores)):
            axes[1, 1].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                           f'{score:.3f}', ha='center', va='bottom', fontsize=10)
        
        axes[1, 1].set_xlabel('CV Fold')
        axes[1, 1].set_ylabel('ROC AUC Score')
        axes[1, 1].set_title('Cross-Validation Performance (PLS-DA)')
        axes[1, 1].set_ylim(0, 1.1)
        axes[1, 1].legend()
        axes[1, 1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        if save_plots:
            plt.savefig(self.output_dir / f'plsda_{target_variable}.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        # Save LV1 weights to CSV with BEC/LEC associations
        weights_output = weights.copy()
        weights_output['Cell_Type_Association'] = 'Neither'
        weights_output.loc[bec_associated.index, 'Cell_Type_Association'] = 'BEC'
        weights_output.loc[lec_associated.index, 'Cell_Type_Association'] = 'LEC'
        
        # Add protein information if available
        if hasattr(self, 'genes_prot_data') and self.genes_prot_data is not None:
            weights_output = weights_output.merge(
                self.genes_prot_data.set_index('Gene names'),
                left_index=True, right_index=True, how='left'
            )
        
        weights_file = self.output_dir / f'plsda_weights_{target_variable}.csv'
        weights_output.to_csv(weights_file)
        logger.info(f"PLS-DA weights saved to {weights_file}")
        
        self.results['plsda'] = plsda_results
        logger.info(f"PLS-DA completed. Found {len(lec_associated)} LEC-associated and "
                   f"{len(bec_associated)} BEC-associated proteins")
        
        return plsda_results
    
    def differential_analysis(self, reference_condition: str, comparison_conditions: List[str],
                            statistical_test: str = 'ttest', multiple_testing_correction: str = 'qvalue',
                            save_results: bool = True) -> Dict:
        """
        Perform differential protein expression analysis.
        
        Parameters:
        -----------
        reference_condition : str
            Reference condition for comparison
        comparison_conditions : list
            List of conditions to compare against reference
        statistical_test : str, default 'ttest'
            Statistical test to use ('ttest', 'mannwhitney')
        multiple_testing_correction : str, default 'qvalue'
            Method for multiple testing correction
        save_results : bool, default True
            Whether to save results to files
            
        Returns:
        --------
        dict : Differential analysis results
        """
        logger.info(f"Performing differential analysis: {comparison_conditions} vs {reference_condition}")
        
        results = {}
        
        for condition in comparison_conditions:
            logger.info(f"Analyzing {condition} vs {reference_condition}")
            
            # Get sample indices
            ref_samples = self.metadata[self.metadata['condition'] == reference_condition]['sample_id']
            comp_samples = self.metadata[self.metadata['condition'] == condition]['sample_id']
            
            if len(ref_samples) == 0 or len(comp_samples) == 0:
                logger.warning(f"No samples found for {condition} or {reference_condition}")
                continue
            
            # Extract data
            ref_data = self.clean_data[ref_samples]
            comp_data = self.clean_data[comp_samples]
            
            # Calculate fold changes and statistics
            ref_means = ref_data.mean(axis=1)
            comp_means = comp_data.mean(axis=1)
            log2_fc = np.log2(comp_means / ref_means)
            
            # Statistical testing
            if statistical_test == 'ttest':
                statistics, p_values = zip(*[
                    stats.ttest_ind(ref_data.iloc[i], comp_data.iloc[i]) 
                    for i in range(len(ref_data))
                ])
            elif statistical_test == 'mannwhitney':
                statistics, p_values = zip(*[
                    stats.mannwhitneyu(ref_data.iloc[i], comp_data.iloc[i], 
                                     alternative='two-sided') 
                    for i in range(len(ref_data))
                ])
            
            # Multiple testing correction
            if multiple_testing_correction == 'qvalue':
                q_values = qvalue(p_values)[1]
            elif multiple_testing_correction == 'bonferroni':
                q_values = np.array(p_values) * len(p_values)
                q_values = np.minimum(q_values, 1.0)
            elif multiple_testing_correction == 'fdr_bh':
                from statsmodels.stats.multitest import multipletests
                _, q_values, _, _ = multipletests(p_values, method='fdr_bh')
            
            # Create results dataframe
            diff_results = pd.DataFrame({
                'gene_name': ref_data.index,
                'log2_fc': log2_fc,
                'p_value': p_values,
                'q_value': q_values,
                'neg_log10_qvalue': -np.log10(q_values),
                'reference_mean': ref_means,
                'comparison_mean': comp_means,
                'significant': (np.abs(log2_fc) >= 0.5) & (q_values <= 0.05)
            })
            
            # Add protein information
            if hasattr(self, 'genes_prot_data') and self.genes_prot_data is not None:
                diff_results = diff_results.merge(
                    self.genes_prot_data, 
                    left_on='gene_name', 
                    right_on='Gene names', 
                    how='left'
                )
            else:
                logger.warning("No protein metadata available for differential analysis")
            
            results[f"{condition}_vs_{reference_condition}"] = diff_results
            
            # Create volcano plot
            self._create_volcano_plot(diff_results, f"{condition} vs {reference_condition}")
            
            if save_results:
                output_file = self.output_dir / f"differential_analysis_{condition}_vs_{reference_condition}.csv"
                diff_results.to_csv(output_file, index=False)
                logger.info(f"Results saved to {output_file}")
        
        self.results['differential_analysis'] = results
        return results
    
    def _create_volcano_plot(self, diff_results: pd.DataFrame, title: str):
        """Create volcano plot for differential analysis results."""
        fig = go.Figure()
        
        # Non-significant proteins
        non_sig = diff_results[~diff_results['significant']]
        fig.add_trace(go.Scatter(
            x=non_sig['log2_fc'],
            y=non_sig['neg_log10_qvalue'],
            mode='markers',
            name='Non-significant',
            marker=dict(color='lightgray', size=4),
            text=non_sig['gene_name'],
            hovertemplate='Gene: %{text}<br>log2FC: %{x:.2f}<br>-log10(q): %{y:.2f}'
        ))
        
        # Significant proteins
        sig = diff_results[diff_results['significant']]
        fig.add_trace(go.Scatter(
            x=sig['log2_fc'],
            y=sig['neg_log10_qvalue'],
            mode='markers',
            name='Significant',
            marker=dict(color='red', size=6),
            text=sig['gene_name'],
            hovertemplate='Gene: %{text}<br>log2FC: %{x:.2f}<br>-log10(q): %{y:.2f}'
        ))
        
        # Add significance thresholds
        fig.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="blue", 
                     annotation_text="q-value = 0.05")
        fig.add_vline(x=0.5, line_dash="dash", line_color="blue")
        fig.add_vline(x=-0.5, line_dash="dash", line_color="blue")
        
        fig.update_layout(
            title=f'Volcano Plot: {title}',
            xaxis_title='log2(Fold Change)',
            yaxis_title='-log10(q-value)',
            hovermode='closest',
            showlegend=True
        )
        
        fig.write_html(self.output_dir / f"volcano_plot_{title.replace(' ', '_')}.html")
        fig.show()
    
    def generate_report(self, save_html: bool = True) -> str:
        """
        Generate a comprehensive analysis report.
        
        Parameters:
        -----------
        save_html : bool, default True
            Whether to save report as HTML file
            
        Returns:
        --------
        str : HTML report content
        """
        logger.info("Generating comprehensive analysis report...")
        
        report_html = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>Endothelial Cell Proteomics Analysis Report</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 40px; }}
                .header {{ background-color: #f0f0f0; padding: 20px; border-radius: 5px; }}
                .section {{ margin: 20px 0; }}
                .metric {{ background-color: #e8f4f8; padding: 10px; margin: 5px 0; border-radius: 3px; }}
                table {{ border-collapse: collapse; width: 100%; }}
                th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
                th {{ background-color: #f2f2f2; }}
            </style>
        </head>
        <body>
            <div class="header">
                <h1>Endothelial Cell Proteomics Analysis Report</h1>
                <p>Generated on: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
                <p>Analysis pipeline version: 1.0</p>
            </div>
            
            <div class="section">
                <h2>Dataset Overview</h2>
                <div class="metric">Total proteins identified: {len(self.raw_data)}</div>
                <div class="metric">Proteins after quality filtering: {len(self.clean_data)}</div>
                <div class="metric">Total samples: {len(self.metadata)}</div>
                <div class="metric">Cell lines: {', '.join(self.metadata['cell_line'].unique())}</div>
                <div class="metric">Time points: {', '.join(self.metadata['time_point'].unique())}</div>
            </div>
        """
        
        # Add PCA results if available
        if 'pca' in self.results:
            pca_results = self.results['pca']
            report_html += f"""
            <div class="section">
                <h2>Principal Component Analysis</h2>
                <div class="metric">PC1 explained variance: {pca_results['explained_variance_ratio'][0]:.2%}</div>
                <div class="metric">PC2 explained variance: {pca_results['explained_variance_ratio'][1]:.2%}</div>
                <div class="metric">Total variance explained (PC1+PC2): {pca_results['explained_variance_ratio'][:2].sum():.2%}</div>
            </div>
            """
        
        # Add PLS-DA results if available
        if 'plsda' in self.results:
            plsda_results = self.results['plsda']
            report_html += f"""
            <div class="section">
                <h2>PLS-DA Analysis</h2>
                <div class="metric">LEC-associated proteins: {len(plsda_results.get('significant_proteins', {}).get('LEC_associated', []))}</div>
                <div class="metric">BEC-associated proteins: {len(plsda_results.get('significant_proteins', {}).get('BEC_associated', []))}</div>
            </div>
            """
        
        # Add differential analysis results if available
        if 'differential_analysis' in self.results:
            diff_results = self.results['differential_analysis']
            report_html += """
            <div class="section">
                <h2>Differential Expression Analysis</h2>
                <table>
                    <tr><th>Comparison</th><th>Significant Proteins</th><th>Upregulated</th><th>Downregulated</th></tr>
            """
            
            for comparison, results in diff_results.items():
                sig_proteins = results[results['significant']]
                upregulated = len(sig_proteins[sig_proteins['log2_fc'] > 0])
                downregulated = len(sig_proteins[sig_proteins['log2_fc'] < 0])
                
                report_html += f"""
                    <tr>
                        <td>{comparison}</td>
                        <td>{len(sig_proteins)}</td>
                        <td>{upregulated}</td>
                        <td>{downregulated}</td>
                    </tr>
                """
            
            report_html += "</table></div>"
        
        report_html += """
        </body>
        </html>
        """
        
        if save_html:
            report_file = self.output_dir / "analysis_report.html"
            with open(report_file, 'w') as f:
                f.write(report_html)
            logger.info(f"Report saved to {report_file}")
        
        return report_html
    
    def run_analysis(self, time_points: List[str] = None, plsda_threshold: float = 90.0) -> Dict:
        """
        Run the complete analysis pipeline.
        
        Parameters:
        -----------
        time_points : list, optional
            List of time points to include in PCA (e.g., ['D2', 'D5'])
        plsda_threshold : float, default 90.0
            Percentile threshold for PLS-DA significant weights
            
        Returns:
        --------
        dict : Dictionary containing all analysis results
        """
        logger.info("Starting complete proteomics analysis pipeline...")
        
        # Ensure data is loaded
        if self.raw_data is None:
            logger.error("No data loaded. Please use load_data() first.")
            return {}
        
        results = {}
        
        try:
            # 1. Preprocess and normalize data
            logger.info("Step 1: Preprocessing and normalizing data...")
            self.preprocess_data()
            self.normalize_data()
            results['preprocessing'] = {'status': 'completed'}
            
            # 2. Quality control
            logger.info("Step 2: Quality control analysis...")
            results['quality_control'] = self.quality_control()
            
            # 3. PCA analysis
            logger.info("Step 3: PCA analysis...")
            if time_points:
                results['pca'] = self.perform_pca(time_points=time_points)
                logger.info(f"PCA completed for time points: {time_points}")
            else:
                results['pca'] = self.perform_pca()
                logger.info("PCA completed for all samples")
            
            # 4. PLS-DA analysis
            logger.info("Step 4: PLS-DA analysis...")
            results['plsda'] = self.perform_plsda(threshold_percentile=plsda_threshold)
            
            if results['plsda'] and 'significant_proteins' in results['plsda']:
                sig_prots = results['plsda']['significant_proteins']
                bec_count = len(sig_prots.get('BEC_associated', []))
                lec_count = len(sig_prots.get('LEC_associated', []))
                logger.info(f"PLS-DA completed: {bec_count} BEC-associated, {lec_count} LEC-associated proteins")
            
            # 5. Generate report
            logger.info("Step 5: Generating analysis report...")
            results['report'] = self.generate_report()
            
            # Store results in the instance
            self.results.update(results)
            
            logger.info("✅ Complete analysis pipeline finished successfully!")
            return results
            
        except Exception as e:
            logger.error(f"Error in analysis pipeline: {e}")
            raise


def main():
    """
    Main analysis pipeline for endothelial cell proteomics data.
    
    This function demonstrates the complete analysis workflow including:
    - Data loading and preprocessing
    - Quality control
    - Normalization
    - PCA analysis
    - PLS-DA for cell type classification
    - Differential expression analysis
    - Report generation
    """
    # Initialize analyzer
    analyzer = ProteomicsAnalyzer(data_path=".", output_dir="results")
    
    try:
        # Load and preprocess data
        analyzer.load_data("Suppl_table_1_Px_data.xlsx")
        analyzer.preprocess_data()
        analyzer.normalize_data(method='zscore')
        
        # Quality control
        qc_results = analyzer.quality_control()
        
        # PCA analysis with D2 and D5 samples only
        pca_results = analyzer.perform_pca(n_components=10, time_points=['D2', 'D5'])
        
        # PLS-DA for cell type discrimination with 90% threshold
        plsda_results = analyzer.perform_plsda(target_variable='cell_type', 
                                             n_permutations=10000, 
                                             threshold_percentile=90.0)
        
        # Differential analysis example
        # Compare different time points within each cell line
        cell_lines = analyzer.metadata['cell_line'].unique()
        time_points = analyzer.metadata['time_point'].unique()
        
        if 'D2' in time_points and 'D5' in time_points:
            for cell_line in cell_lines:
                ref_condition = f"{cell_line}_D2"
                comp_condition = f"{cell_line}_D5"
                
                if (ref_condition in analyzer.metadata['condition'].values and 
                    comp_condition in analyzer.metadata['condition'].values):
                    
                    diff_results = analyzer.differential_analysis(
                        reference_condition=ref_condition,
                        comparison_conditions=[comp_condition]
                    )
        
        # Generate comprehensive report
        report = analyzer.generate_report()
        
        print("Analysis completed successfully!")
        print(f"Results saved to: {analyzer.output_dir}")
        
    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        raise


if __name__ == "__main__":
    main()