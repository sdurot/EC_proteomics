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
            
        except FileNotFoundError:
            logger.error(f"File {excel_file} not found in {self.data_path}")
            raise
        except Exception as e:
            logger.error(f"Error loading data: {e}")
            raise
            
        return self
    
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
    
    def perform_pca(self, n_components: int = 10, save_plots: bool = True) -> Dict:
        """
        Perform Principal Component Analysis.
        
        Parameters:
        -----------
        n_components : int, default 10
            Number of principal components to calculate
        save_plots : bool, default True
            Whether to save plots to output directory
            
        Returns:
        --------
        dict : PCA results including components, explained variance, and scores
        """
        logger.info(f"Performing PCA with {n_components} components...")
        
        # Prepare data for PCA
        pca_data = self.normalized_data.T  # Samples as rows, proteins as columns
        
        # Perform PCA
        pca = PCA(n_components=n_components)
        pca_scores = pca.fit_transform(pca_data)
        
        # Create results dictionary
        pca_results = {
            'pca_model': pca,
            'scores': pd.DataFrame(
                pca_scores, 
                index=self.normalized_data.columns,
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
        pca_df = pca_df.merge(self.metadata.set_index('sample_id'), 
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
                      n_permutations: int = 1000, save_plots: bool = True) -> Dict:
        """
        Perform Partial Least Squares Discriminant Analysis (PLS-DA).
        
        Parameters:
        -----------
        target_variable : str, default 'cell_type'
            Target variable for classification ('cell_type', 'cell_line', etc.)
        n_components : int, default 2
            Number of PLS components
        n_permutations : int, default 1000
            Number of permutations for significance testing
        save_plots : bool, default True
            Whether to save plots to output directory
            
        Returns:
        --------
        dict : PLS-DA results including model, weights, and significance
        """
        logger.info(f"Performing PLS-DA for {target_variable}...")
        
        # Prepare data
        X = self.normalized_data.T.fillna(0)  # Samples x Proteins
        y = pd.get_dummies(self.metadata.set_index('sample_id')[target_variable]).iloc[:, 0]
        
        # Fit PLS-DA model
        plsda = PLSRegression(n_components=n_components, scale=True)
        plsda.fit(X, y)
        
        # Get scores and weights
        scores = plsda.transform(X)
        weights = pd.DataFrame(
            plsda.x_weights_, 
            index=self.normalized_data.index,
            columns=[f'LV{i+1}' for i in range(n_components)]
        )
        
        # Permutation testing
        logger.info(f"Running {n_permutations} permutations for significance testing...")
        null_weights = []
        
        for _ in range(n_permutations):
            y_perm = np.random.permutation(y)
            plsda_perm = PLSRegression(n_components=n_components, scale=True)
            plsda_perm.fit(X, y_perm)
            null_weights.append(plsda_perm.x_weights_[:, 0])
        
        null_weights = np.array(null_weights)
        
        # Calculate significance thresholds
        significance_thresholds = {
            'upper_95': np.percentile(null_weights, 95, axis=0),
            'lower_5': np.percentile(null_weights, 5, axis=0),
            'upper_99': np.percentile(null_weights, 99, axis=0),
            'lower_1': np.percentile(null_weights, 1, axis=0)
        }
        
        # Identify significant proteins
        lv1_weights = weights['LV1']
        significant_proteins = {
            'high_weights': lv1_weights[lv1_weights >= np.mean(significance_thresholds['upper_95'])],
            'low_weights': lv1_weights[lv1_weights <= np.mean(significance_thresholds['lower_5'])]
        }
        
        plsda_results = {
            'model': plsda,
            'scores': pd.DataFrame(scores, index=X.index, columns=[f'LV{i+1}' for i in range(n_components)]),
            'weights': weights,
            'significant_proteins': significant_proteins,
            'significance_thresholds': significance_thresholds,
            'null_distribution': null_weights
        }
        
        # Create visualization
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # Score plot
        scores_df = plsda_results['scores'].copy()
        scores_df = scores_df.merge(self.metadata.set_index('sample_id'), 
                                  left_index=True, right_index=True)
        
        for group in scores_df[target_variable].unique():
            mask = scores_df[target_variable] == group
            axes[0, 0].scatter(scores_df.loc[mask, 'LV1'], scores_df.loc[mask, 'LV2'], 
                             label=group, s=60, alpha=0.7)
        
        axes[0, 0].set_xlabel('LV1')
        axes[0, 0].set_ylabel('LV2')
        axes[0, 0].set_title('PLS-DA Score Plot')
        axes[0, 0].legend()
        axes[0, 0].grid(True, alpha=0.3)
        
        # Weights plot with significance
        sorted_weights = lv1_weights.sort_values()
        axes[0, 1].scatter(range(len(sorted_weights)), sorted_weights, 
                          s=20, alpha=0.6, color='black')
        axes[0, 1].axhline(np.mean(significance_thresholds['upper_95']), 
                          color='red', linestyle='--', label='95% threshold')
        axes[0, 1].axhline(np.mean(significance_thresholds['lower_5']), 
                          color='red', linestyle='--')
        axes[0, 1].set_xlabel('Proteins (ranked)')
        axes[0, 1].set_ylabel('LV1 Weights')
        axes[0, 1].set_title('PLS-DA Weights with Significance Thresholds')
        axes[0, 1].legend()
        axes[0, 1].grid(True, alpha=0.3)
        
        # Null distribution
        axes[1, 0].hist(null_weights.flatten(), bins=50, alpha=0.7, 
                       color='lightblue', label='Null distribution')
        axes[1, 0].axvline(len(significant_proteins['high_weights']) + 
                          len(significant_proteins['low_weights']), 
                          color='red', linestyle='--', 
                          label=f"Observed ({len(significant_proteins['high_weights']) + len(significant_proteins['low_weights'])})")
        axes[1, 0].set_xlabel('Number of significant proteins')
        axes[1, 0].set_ylabel('Frequency')
        axes[1, 0].set_title('Permutation Test Results')
        axes[1, 0].legend()
        axes[1, 0].grid(True, alpha=0.3)
        
        # Cross-validation scores
        cv_scores = cross_val_score(plsda, X, y, cv=5, scoring='roc_auc')
        axes[1, 1].bar(range(len(cv_scores)), cv_scores)
        axes[1, 1].axhline(cv_scores.mean(), color='red', linestyle='--', 
                          label=f'Mean CV Score: {cv_scores.mean():.3f}')
        axes[1, 1].set_xlabel('CV Fold')
        axes[1, 1].set_ylabel('ROC AUC')
        axes[1, 1].set_title('Cross-Validation Performance')
        axes[1, 1].legend()
        axes[1, 1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        if save_plots:
            plt.savefig(self.output_dir / f'plsda_{target_variable}.png')
        plt.show()
        
        self.results['plsda'] = plsda_results
        logger.info(f"PLS-DA completed. Found {len(significant_proteins['high_weights'])} high and "
                   f"{len(significant_proteins['low_weights'])} low weighted significant proteins")
        
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
            protein_info_subset = self.raw_data[['Gene names', 'Protein description', 'Protein.Ids']].drop_duplicates()
            diff_results = diff_results.merge(
                protein_info_subset, 
                left_on='gene_name', 
                right_on='Gene names', 
                how='left'
            )
            
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
                <div class="metric">Significant proteins (high weights): {len(plsda_results['significant_proteins']['high_weights'])}</div>
                <div class="metric">Significant proteins (low weights): {len(plsda_results['significant_proteins']['low_weights'])}</div>
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
        analyzer.load_data("Suppl_table_1.xlsx")
        analyzer.preprocess_data()
        analyzer.normalize_data(method='zscore')
        
        # Quality control
        qc_results = analyzer.quality_control()
        
        # PCA analysis
        pca_results = analyzer.perform_pca(n_components=10)
        
        # PLS-DA for cell type discrimination
        plsda_results = analyzer.perform_plsda(target_variable='cell_type', n_permutations=1000)
        
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