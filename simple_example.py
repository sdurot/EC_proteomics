#!/usr/bin/env python3
"""
Simple example script demonstrating the EC Proteomics Analysis Pipeline.

This script shows how to use the proteomics analysis pipeline with minimal setup.
Run this script after installing the requirements and having your data file ready.

Usage:
    python simple_example.py

Make sure you have:
1. Installed all requirements: pip install -r requirements.txt
2. Your data file (Suppl_table_1.xlsx) in the same directory
"""

import sys
from pathlib import Path

# Add the current directory to Python path for imports
sys.path.append(str(Path(__file__).parent))

try:
    from proteomics_analysis import ProteomicsAnalyzer
    print("✅ Successfully imported ProteomicsAnalyzer")
except ImportError as e:
    print(f"❌ Error importing ProteomicsAnalyzer: {e}")
    print("Make sure you have installed all requirements:")
    print("pip install -r requirements.txt")
    sys.exit(1)

def main():
    """Run a simple proteomics analysis example."""
    
    print("🧬 Starting EC Proteomics Analysis Pipeline Example")
    print("=" * 60)
    
    # Check if data file exists
    data_file = "Suppl_table_1.xlsx"
    if not Path(data_file).exists():
        print(f"❌ Data file '{data_file}' not found!")
        print("Please make sure your Excel data file is in the current directory.")
        print("The file should contain:")
        print("  - Sheet 'Px_LFQ_data': Main proteomics data")
        print("  - Sheet 'protein_genes_data': Protein information")
        return
    
    try:
        # Step 1: Initialize the analyzer
        print("🔧 Initializing ProteomicsAnalyzer...")
        analyzer = ProteomicsAnalyzer(
            data_path=".",
            output_dir="results_example"
        )
        
        # Step 2: Load and preprocess data
        print("📊 Loading and preprocessing data...")
        analyzer.load_data(data_file)
        analyzer.preprocess_data()
        
        print(f"   - Raw data: {analyzer.raw_data.shape[0]} proteins, {analyzer.raw_data.shape[1]-3} samples")
        print(f"   - Clean data: {analyzer.clean_data.shape[0]} proteins")
        print(f"   - Cell lines: {', '.join(analyzer.metadata['cell_line'].unique())}")
        print(f"   - Time points: {', '.join(analyzer.metadata['time_point'].unique())}")
        
        # Step 3: Normalize data
        print("🎯 Normalizing data...")
        analyzer.normalize_data(method='zscore')
        
        # Step 4: Quality control
        print("🔍 Performing quality control...")
        qc_results = analyzer.quality_control(save_plots=True)
        
        # Step 5: PCA analysis
        print("📈 Running Principal Component Analysis...")
        pca_results = analyzer.perform_pca(n_components=5, save_plots=True)
        
        pc1_var = pca_results['explained_variance_ratio'][0] * 100
        pc2_var = pca_results['explained_variance_ratio'][1] * 100
        print(f"   - PC1 explains {pc1_var:.1f}% of variance")
        print(f"   - PC2 explains {pc2_var:.1f}% of variance")
        
        # Step 6: PLS-DA analysis
        print("🎯 Running PLS-DA for cell type classification...")
        plsda_results = analyzer.perform_plsda(
            target_variable='cell_type',
            n_components=2,
            n_permutations=100,  # Reduced for faster example
            save_plots=True
        )
        
        n_lec = len(plsda_results['significant_proteins']['high_weights'])
        n_bec = len(plsda_results['significant_proteins']['low_weights'])
        print(f"   - Found {n_lec} LEC-associated proteins")
        print(f"   - Found {n_bec} BEC-associated proteins")
        
        # Step 7: Differential analysis (example)
        print("📊 Running differential expression analysis...")
        
        # Example: Compare D5 vs D2 for HUVEC
        if 'HUVEC_D2' in analyzer.metadata['condition'].values and 'HUVEC_D5' in analyzer.metadata['condition'].values:
            diff_results = analyzer.differential_analysis(
                reference_condition="HUVEC_D2",
                comparison_conditions=["HUVEC_D5"],
                statistical_test='ttest',
                save_results=True
            )
            
            huvec_results = diff_results['HUVEC_D5_vs_HUVEC_D2']
            significant = huvec_results[huvec_results['significant']]
            upregulated = significant[significant['log2_fc'] > 0]
            downregulated = significant[significant['log2_fc'] < 0]
            
            print(f"   - HUVEC D5 vs D2: {len(significant)} significant proteins")
            print(f"     * Upregulated: {len(upregulated)}")
            print(f"     * Downregulated: {len(downregulated)}")
        
        # Step 8: Generate report
        print("📋 Generating comprehensive report...")
        report = analyzer.generate_report(save_html=True)
        
        # Summary
        print("\n🎉 Analysis completed successfully!")
        print("=" * 60)
        print(f"📁 Results saved to: {analyzer.output_dir}")
        print("\n📄 Generated files:")
        
        result_files = [
            "analysis_report.html - Comprehensive analysis report",
            "quality_control.png - Quality control diagnostic plots",
            "pca_analysis.html - Interactive PCA visualization", 
            "plsda_cell_type.png - PLS-DA classification results",
            "*.csv - Data tables for further analysis"
        ]
        
        for file_desc in result_files:
            print(f"   • {file_desc}")
        
        print(f"\n🌐 Open 'analysis_report.html' in your browser to view the complete analysis!")
        
        # Example of accessing results programmatically
        print("\n🔬 Quick results access example:")
        print("# Access PCA scores")
        print("pca_scores = analyzer.results['pca']['scores']")
        print("# Access significant proteins")
        print("lec_proteins = analyzer.results['plsda']['significant_proteins']['high_weights']")
        print("# Access differential analysis")
        print("diff_data = analyzer.results['differential_analysis']")
        
    except FileNotFoundError as e:
        print(f"❌ File not found: {e}")
        print("Make sure your data file is in the correct location.")
    except Exception as e:
        print(f"❌ Analysis failed: {e}")
        print("Check your data format and try again.")
        raise  # Re-raise for debugging

if __name__ == "__main__":
    main()