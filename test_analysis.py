#!/usr/bin/env python3
"""
Unit tests for the EC Proteomics Analysis Pipeline.

Run these tests to ensure the analysis pipeline is working correctly.

Usage:
    python test_analysis.py
    or
    pytest test_analysis.py
"""

import unittest
import numpy as np
import pandas as pd
import tempfile
import shutil
from pathlib import Path
import sys

# Add the current directory to Python path
sys.path.append(str(Path(__file__).parent))

from proteomics_analysis import ProteomicsAnalyzer

class TestProteomicsAnalyzer(unittest.TestCase):
    """Test cases for the ProteomicsAnalyzer class."""
    
    def setUp(self):
        """Set up test fixtures before each test method."""
        self.temp_dir = tempfile.mkdtemp()
        self.analyzer = ProteomicsAnalyzer(
            data_path=self.temp_dir,
            output_dir=str(Path(self.temp_dir) / "test_results")
        )
        
        # Create mock data
        self.create_mock_data()
    
    def tearDown(self):
        """Clean up after each test method."""
        shutil.rmtree(self.temp_dir)
    
    def create_mock_data(self):
        """Create mock proteomics data for testing."""
        # Create mock intensity data
        np.random.seed(42)
        n_proteins = 100
        n_samples = 24  # 4 cell lines × 3 time points × 2 replicates
        
        # Generate sample names
        cell_lines = ['HDBEC', 'HDLEC', 'HUVEC', 'iLEC']
        time_points = ['D2', 'D5', 'D7']
        sample_names = []
        
        for cell_line in cell_lines:
            for time_point in time_points:
                for replicate in [1, 2]:
                    sample_names.append(f"{cell_line}_{time_point}_{replicate}")
        
        # Generate mock protein data
        gene_names = [f"GENE{i:03d}" for i in range(n_proteins)]
        protein_ids = [f"P{i:05d}" for i in range(n_proteins)]
        protein_descriptions = [f"Protein {i} description" for i in range(n_proteins)]
        
        # Generate intensity data with some structure
        intensities = np.random.lognormal(mean=10, sigma=1, size=(n_proteins, n_samples))
        
        # Add some missing values
        missing_mask = np.random.random((n_proteins, n_samples)) < 0.1
        intensities[missing_mask] = np.nan
        
        # Create main data DataFrame
        main_data = pd.DataFrame(
            np.column_stack([
                gene_names,
                gene_names,  # duplicate for index
                intensities
            ]),
            columns=['Gene names', 'Gene names'] + sample_names
        )
        
        # Convert intensity columns to numeric
        for col in sample_names:
            main_data[col] = pd.to_numeric(main_data[col], errors='coerce')
        
        # Create protein info DataFrame
        protein_info = pd.DataFrame({
            'Gene names': gene_names,
            'First.Protein.Description': protein_descriptions,
            'Protein.Ids': protein_ids
        })
        
        # Save to Excel file
        excel_path = Path(self.temp_dir) / "test_data.xlsx"
        with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
            main_data.to_excel(writer, sheet_name='Px_LFQ_data', index=False)
            protein_info.to_excel(writer, sheet_name='protein_genes_data', index=False)
        
        self.test_excel_file = "test_data.xlsx"
    
    def test_initialization(self):
        """Test ProteomicsAnalyzer initialization."""
        self.assertIsInstance(self.analyzer, ProteomicsAnalyzer)
        self.assertEqual(str(self.analyzer.data_path), self.temp_dir)
        self.assertTrue(self.analyzer.output_dir.exists())
    
    def test_load_data(self):
        """Test data loading functionality."""
        self.analyzer.load_data(self.test_excel_file)
        
        self.assertIsNotNone(self.analyzer.raw_data)
        self.assertIsNotNone(self.analyzer.protein_info)
        
        # Check data structure
        self.assertIn('Gene names', self.analyzer.raw_data.columns)
        self.assertIn('Protein description', self.analyzer.raw_data.columns)
        self.assertEqual(len(self.analyzer.raw_data), 100)  # 100 proteins
    
    def test_preprocess_data(self):
        """Test data preprocessing."""
        self.analyzer.load_data(self.test_excel_file)
        self.analyzer.preprocess_data()
        
        self.assertIsNotNone(self.analyzer.clean_data)
        self.assertIsNotNone(self.analyzer.metadata)
        
        # Check metadata structure
        expected_columns = ['sample_id', 'cell_line', 'time_point', 'condition', 'cell_type']
        for col in expected_columns:
            self.assertIn(col, self.analyzer.metadata.columns)
        
        # Check cell type assignment
        bec_samples = self.analyzer.metadata[self.analyzer.metadata['cell_type'] == 'BEC']
        lec_samples = self.analyzer.metadata[self.analyzer.metadata['cell_type'] == 'LEC']
        
        self.assertTrue(len(bec_samples) > 0)
        self.assertTrue(len(lec_samples) > 0)
    
    def test_normalize_data(self):
        """Test data normalization."""
        self.analyzer.load_data(self.test_excel_file)
        self.analyzer.preprocess_data()
        self.analyzer.normalize_data(method='zscore')
        
        self.assertIsNotNone(self.analyzer.normalized_data)
        
        # Check that data is normalized (approximately zero mean, unit variance)
        sample_data = self.analyzer.normalized_data.iloc[:, 0].dropna()
        self.assertAlmostEqual(sample_data.mean(), 0, places=10)
        self.assertAlmostEqual(sample_data.std(), 1, places=10)
    
    def test_quality_control(self):
        """Test quality control analysis."""
        self.analyzer.load_data(self.test_excel_file)
        self.analyzer.preprocess_data()
        self.analyzer.normalize_data()
        
        qc_results = self.analyzer.quality_control(save_plots=False)
        
        self.assertIsInstance(qc_results, dict)
        self.assertIn('completeness', qc_results)
        self.assertIn('sample_stats', qc_results)
    
    def test_pca(self):
        """Test PCA analysis."""
        self.analyzer.load_data(self.test_excel_file)
        self.analyzer.preprocess_data()
        self.analyzer.normalize_data()
        
        pca_results = self.analyzer.perform_pca(n_components=5, save_plots=False)
        
        self.assertIsInstance(pca_results, dict)
        self.assertIn('pca_model', pca_results)
        self.assertIn('scores', pca_results)
        self.assertIn('explained_variance_ratio', pca_results)
        
        # Check dimensions
        self.assertEqual(len(pca_results['explained_variance_ratio']), 5)
        self.assertEqual(pca_results['scores'].shape[1], 5)
    
    def test_plsda(self):
        """Test PLS-DA analysis."""
        self.analyzer.load_data(self.test_excel_file)
        self.analyzer.preprocess_data()
        self.analyzer.normalize_data()
        
        plsda_results = self.analyzer.perform_plsda(
            target_variable='cell_type',
            n_components=2,
            n_permutations=10,  # Small number for testing
            save_plots=False
        )
        
        self.assertIsInstance(plsda_results, dict)
        self.assertIn('model', plsda_results)
        self.assertIn('scores', plsda_results)
        self.assertIn('weights', plsda_results)
        self.assertIn('significant_proteins', plsda_results)
    
    def test_differential_analysis(self):
        """Test differential expression analysis."""
        self.analyzer.load_data(self.test_excel_file)
        self.analyzer.preprocess_data()
        self.analyzer.normalize_data()
        
        # Test with available conditions
        available_conditions = self.analyzer.metadata['condition'].unique()
        if len(available_conditions) >= 2:
            ref_condition = available_conditions[0]
            comp_condition = available_conditions[1]
            
            diff_results = self.analyzer.differential_analysis(
                reference_condition=ref_condition,
                comparison_conditions=[comp_condition],
                save_results=False
            )
            
            self.assertIsInstance(diff_results, dict)
            comparison_key = f"{comp_condition}_vs_{ref_condition}"
            self.assertIn(comparison_key, diff_results)
            
            result_df = diff_results[comparison_key]
            expected_columns = ['gene_name', 'log2_fc', 'p_value', 'q_value', 'significant']
            for col in expected_columns:
                self.assertIn(col, result_df.columns)

class TestConfigurationAndSetup(unittest.TestCase):
    """Test configuration and setup functionality."""
    
    def test_imports(self):
        """Test that all required modules can be imported."""
        try:
            import numpy
            import pandas
            import matplotlib.pyplot
            import seaborn
            import scipy.stats
            import sklearn.decomposition
            import sklearn.cross_decomposition
            import plotly.express
            self.assertTrue(True)
        except ImportError as e:
            self.fail(f"Failed to import required module: {e}")
    
    def test_config_import(self):
        """Test that configuration can be imported and used."""
        try:
            from config import get_config, update_config
            
            config = get_config()
            self.assertIsInstance(config, dict)
            self.assertIn('data', config)
            self.assertIn('normalization', config)
            
            # Test config update
            updated_config = update_config({'data': {'excel_file': 'test.xlsx'}})
            self.assertEqual(updated_config['data']['excel_file'], 'test.xlsx')
            
        except ImportError as e:
            self.fail(f"Failed to import config: {e}")

def run_tests():
    """Run all tests."""
    print("🧪 Running EC Proteomics Analysis Pipeline Tests")
    print("=" * 60)
    
    # Create test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add test cases
    suite.addTests(loader.loadTestsFromTestCase(TestProteomicsAnalyzer))
    suite.addTests(loader.loadTestsFromTestCase(TestConfigurationAndSetup))
    
    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    # Print summary
    print("\n" + "=" * 60)
    if result.wasSuccessful():
        print("✅ All tests passed successfully!")
        print("The proteomics analysis pipeline is working correctly.")
    else:
        print("❌ Some tests failed.")
        print(f"Failures: {len(result.failures)}")
        print(f"Errors: {len(result.errors)}")
    
    return result.wasSuccessful()

if __name__ == "__main__":
    success = run_tests()
    sys.exit(0 if success else 1)