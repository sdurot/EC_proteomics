# -*- coding: utf-8 -*-
"""
@author: sdurot

"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns
import matplotlib

# Set default matplotlib settings
from matplotlib.pyplot import figure
matplotlib.rc_file_defaults()
matplotlib.rcParams.update({'font.size': 16})
matplotlib.rcParams['font.family'] = 'Arial'

# Change working directory to the project folder
os.chdir('project folder with result file')

# Load and clean proteomics data, clean meaning that an outlier (HDLEC D6 2) was removed. The outlier was 
# detected already using diaNN output 
prot_data = pd.read_excel('Suppl_table_1.xlsx', sheet_name='Px_LFQ_data').rename(columns={'newColName':'Gene names'})

# Information about protein groups (PGs)
prot_genes_data = pd.read_excel('Suppl_table_1.xlsx', sheet_name='protein_genes_data').rename(columns={'Genes':'Gene names'})

# Merging of the two datasets, creation of a 'clean' dataset that has only PGs that have values in all samples
prot_data = prot_data.merge(prot_genes_data, on='Gene names').rename(columns={'First.Protein:Description':'Protein description'})
prot_data_clean = prot_data.dropna()
prot_data_clean.index = prot_data_clean.iloc[:, 1]

# Create DataFrame with gene names, protein descriptions, and protein IDs
genes_prot_data = pd.DataFrame({
    'Gene names': prot_data['Gene names'],
    'Protein description': prot_data['First.Protein.Description'],
    'Protein IDs': prot_data['Protein.Ids']
}).reset_index(drop=True)

# Z-score normalization of the normal and clean datasets
prot_data_clean_zscored = prot_data_clean.iloc[:, 2:-2].T.apply(stats.zscore).T
prot_data_zscored = prot_data.iloc[:, 2:-2].T.apply(stats.zscore).T

# Create metadata DataFrame and metadata lists from column names to facilitate the selection of certain samples later
metadata = pd.DataFrame(prot_data_clean_zscored.columns.str.split('_').tolist())
metadata['samples'] = metadata.iloc[:, 0] + " " + metadata.iloc[:, 1]
celllines = pd.factorize(metadata.iloc[:, 0])[0]
samples = pd.factorize(metadata.iloc[:, 3])[0]
days = pd.factorize(metadata.iloc[:, 1])[0]

# Discriminate into LECs and BECs for later analysis
celltype = np.array([0 if c in ['HDBEC', 'HUVEC'] else 1 for c in metadata.iloc[:, 0]])

# Heatmap of z-scored data to check whether there are unwanted patterns, i.e. batch effects etc.
figure(figsize=(6, 4), dpi=500)
sns.heatmap(prot_data_zscored.set_axis([np.arange(0, 7894, 1)], axis=0), vmin=-3, vmax=3, cmap='RdBu', annot=False)
plt.title('z-scored MaxLFQ normalized data')
plt.show()

# Correlation of day 2 (proliferation) and day 5 (quiescence) samples 
days = ['2D', '5D']
celllines = ['HDBEC', 'HDLEC', 'HUVEC', 'iLEC']
prot_celllines_days = np.concatenate([[f'{day} {cellline}']*3 for day in days for cellline in celllines])

prot_samples_day2 = prot_data_clean_zscored.iloc[:, metadata['samples'].str.contains('D2').tolist()]
prot_samples_day5 = prot_data_clean_zscored.iloc[:, metadata['samples'].str.contains('D5').tolist()]
prot_samples_day2and5 = pd.concat([prot_samples_day2, prot_samples_day5], axis=1)

prot_samples_day2and5_correlation = prot_samples_day2and5.corr(method='spearman')
prot_samples_day2and5_correlation.index = prot_celllines_days
prot_samples_day2and5_correlation.columns = prot_celllines_days

figure(figsize=(6, 6), dpi=500)
sns.heatmap(prot_samples_day2and5_correlation, vmin=-1, vmax=1, cmap='RdBu', annot=False)
plt.title('Spearman correlation of proteomics samples')
plt.show()


#%% Create datasets with average protein expression
# Function to compute average protein expression
def average_protein_expression(data):
    """
    Calculate mean protein expression across replicates.

    Parameters:
    data: proteomics dataset with samples as columns and protein intenities as rows.

    """
    data = data.reset_index(drop=True)
    average_exp_df = pd.DataFrame(index=range(len(data)), columns=samples_pert)
    average_exp_df_std = pd.DataFrame(index=range(len(data)), columns=samples_pert)

    # Calculate mean and standard deviation for each condition
    for i in samples_pert_factor:
        columns = np.where(samples == i)[0]
        means = data.iloc[:, columns].mean(axis=1)
        std = data.iloc[:, columns].std(axis=1)
        average_exp_df.iloc[:, i] = means
        average_exp_df_std.iloc[:, i] = std

    # Z-score the data and plot heatmap
    average_exp_df_zscored = average_exp_df.T.apply(stats.zscore).T
    sns.heatmap(average_exp_df_zscored, vmin=-5, vmax=5, cmap='RdBu', annot=False)
    plt.show()

    # Spearman correlation and heatmap of raw values
    average_exp_correlation = average_exp_df.corr(method='spearman').replace(np.nan, 0)
    sns.heatmap(average_exp_correlation, vmin=-1, vmax=1, cmap='RdBu', annot=False)
    plt.show()

    sns.clustermap(average_exp_correlation, method="complete", cmap='RdBu', annot=False, figsize=(15, 12))
    plt.show()

    return average_exp_df, average_exp_df_std

# Set up metadata
samples_pert = np.unique(metadata['samples']).astype('str')
samples_pert_factor = np.unique(samples)
samples_pert_metadata = pd.DataFrame(samples_pert).iloc[:, 0].str.split(' ').tolist()

# Process and analyze protein data
prot_data_edit = prot_data
average_prot_exp_full_df, average_prot_exp_full_std_df = average_protein_expression(prot_data_edit.iloc[:, 2:-2])
average_prot_exp_NAomit, average_prot_exp_NAomit_std_df = average_protein_expression(prot_data_clean.iloc[:, 2:-2])

average_prot_exp_full_df.index = prot_data['Gene names']

# Reset the index of prot_data and sort columns for average protein expression DataFrames
prot_data_newindex = prot_data.reset_index(drop=True)
column_order = [0, 1, 2, 3, 4, 5, 12, 13, 14, 15, 16, 17, 6, 7, 8, 9, 10, 11, 18, 19, 20, 21, 22, 23]
average_prot_exp_full_df_sorted = average_prot_exp_full_df.iloc[:, column_order]
average_prot_exp_full_std_df_sorted = average_prot_exp_full_std_df.iloc[:, column_order]


## Function to plot protein intensities for days 2 and 5

def plot_prot_intensities_days257(protein_of_interest, columns=None):
    """
    Plot protein intensities for selected columns with error bars.

    Parameters:
    protein_of_interest (str): The name of the protein to plot.
    columns (list): The list of column indices to include in the plot. 
                   Default is [0, 3, 12, 15, 6, 9, 18, 21].

    Raises:
    ValueError: If protein_of_interest is not found in the data.
    """
    # Constants
    DEFAULT_COLUMNS = [0, 3, 12, 15, 6, 9, 18, 21]
    CELL_TYPES = ['HDBEC', 'HDBEC', 'HUVEC', 'HUVEC', 'HDLEC', 'HDLEC', 'iLEC', 'iLEC']
    DAYS = ['D2', 'D5'] * 4
    
    # Set default columns if none provided
    columns = columns or DEFAULT_COLUMNS

    try:
        # Find the index of the protein of interest
        prot_index = prot_data_newindex.index[
            prot_data_newindex['Gene names'].str.startswith(protein_of_interest)
        ].tolist()[0]
    except IndexError:
        raise ValueError(f"Protein '{protein_of_interest}' not found in the dataset")

    # Select intensities and standard deviations
    prot_intensities = average_prot_exp_full_df.iloc[prot_index, columns].replace(np.nan, 0)
    prot_std = average_prot_exp_full_std_df.iloc[prot_index, columns].replace(np.nan, 0)
    
    # Create plot
    fig, ax = plt.subplots(figsize=(4, 6), dpi=500)
    x = np.arange(len(columns))
    
    # Plot bars with error bars
    ax.bar(x, prot_intensities, yerr=prot_std, capsize=5)
    
    # Customize plot
    x_labels = [f"{day} {cell}" for day, cell in zip(DAYS, CELL_TYPES)]
    ax.set_xticks(x)
    ax.set_xticklabels(x_labels, rotation=90)
    ax.ticklabel_format(axis='y', style='sci', scilimits=[-1, 1], useMathText=True)
    ax.set_title(protein_of_interest)
    ax.grid(axis='y')
    ax.set_ylabel('Intensity', fontweight='bold')
    
    # Adjust layout and save
    plt.tight_layout()
    # plt.savefig(f'{protein_of_interest}_bar_D2and5.svg', format='svg')
    plt.show()

if __name__ == "__main__":
    protein = input('Protein: ')
    plot_prot_intensities_days257(protein)


#%% PCA with samples from specific days
from sklearn.decomposition import PCA

# Function to select indices based on specific days
def get_day_indices(metadata, days):
    return np.sort(np.concatenate([metadata.index[metadata.iloc[:, 1].str.startswith(day)].tolist() for day in days]))

# Function to perform PCA and plot cumulative explained variance
def perform_pca_and_plot(data, n_components=10):
    pca = PCA(n_components).fit(data)
    explained_variance = pca.explained_variance_ratio_
    
    plt.plot(np.cumsum(explained_variance))
    plt.xlabel('Number of Components')
    plt.ylabel('Cumulative Explained Variance')
    plt.show()
    
    return pca, explained_variance

# Function to calculate PCA scores
def calculate_pca_scores(data, pca):
    return data.dot(pca.components_.T).reset_index(drop=True)

# Function to plot PCA scatter plot
def plot_pca_scatter(pca_scores, metadata, indices, pc1_expl, pc2_expl):
    pca_df = pd.DataFrame({
        'PC1': pca_scores.iloc[:, 0],
        'PC2': pca_scores.iloc[:, 1],
        'time points (days)': metadata.iloc[indices, 1].reset_index(drop=True),
        'cell line': metadata.iloc[indices, 0].reset_index(drop=True)
    })
    
    plt.figure(figsize=(4, 4), dpi=500)
    sns.scatterplot(data=pca_df, x='PC1', y='PC2', hue='cell line', style='time points (days)', s=100)
    plt.xlabel(f'PC1: {pc1_expl:.2f}%', fontweight='bold')
    plt.ylabel(f'PC2: {pc2_expl:.2f}%', fontweight='bold')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), markerscale=1.5, fontsize=12)
    plt.title('EC proteomics', fontweight='bold')
    # plt.savefig('pca_proteomics.svg', format='svg')
    # plt.savefig('pca_proteomics.png', format='png', dpi=1200)
    plt.show()

# Function to plot PCA loadings
def plot_pca_loadings(pca, data):
    plt.scatter(pca.components_[0, :], pca.components_[1, :])
    plt.xlabel('Loadings PC1')
    plt.ylabel('Loadings PC2')
    plt.show()

    loadings_df = pd.DataFrame({
        'Loading PC1': pca.components_[0, :],
        'Loading PC2': pca.components_[1, :]
    })
    loadings_df.index = data.index
    return loadings_df

# Select indices for specific days
days_to_select = ['D2', 'D5']
days_indices = get_day_indices(metadata, days_to_select)

# Extract and z-score the data for selected days
zscored_df_days = prot_data_clean_zscored.iloc[:, list(days_indices)]

# Perform PCA and plot cumulative explained variance
pca, explained_variance = perform_pca_and_plot(zscored_df_days.T, n_components=10)
pc1_expl, pc2_expl = explained_variance[0] * 100, explained_variance[1] * 100

# Calculate PCA scores
pca_scores = calculate_pca_scores(zscored_df_days.T, pca)

# Plot PCA scatter plot
plot_pca_scatter(pca_scores, metadata, days_indices, pc1_expl, pc2_expl)

# Plot PCA loadings
loadings_df = plot_pca_loadings(pca, prot_data_clean_zscored)


#%% PLS-DA to find proteins that define BEC and LEC identities
from sklearn.cross_decomposition import PLSRegression

celltype = np.where(metadata.iloc[:,0].isin(['HDBEC', 'HUVEC']), 0, 1)

# PLS-DA with the whole dataset
plsr = PLSRegression(n_components=2, scale=True)
data_for_plsr = prot_data.iloc[:,2:-2].replace(np.nan, 0).T
plsr.fit(data_for_plsr, celltype)

# Extract scores and set index
scores = pd.DataFrame(plsr.x_scores_)
scores.index = prot_data.iloc[:,2:-2].columns

# Create DataFrame for Seaborn scatterplot
pls_df = pd.DataFrame({
    'LV1': scores.iloc[:,0].reset_index(drop=True),
    'LV2': scores.iloc[:,1].reset_index(drop=True),
    'time points (days)': metadata.iloc[:,1].reset_index(drop=True),
    'cell line': metadata.iloc[:,0].reset_index(drop=True)
})

# Plotting PLS-DA scatter plot
plt.figure(figsize=(6, 6), dpi=1200)
sns.scatterplot(data=pls_df, x='LV1', y='LV2', hue='cell line', style='time points (days)', s=100)
plt.xlabel('Component 1', fontweight='bold')
plt.ylabel('Component 2', fontweight='bold')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), markerscale=2, fontsize=12)
# plt.savefig('plsda.svg', format='svg')
plt.show()

# Extract weights
pls_weights = pd.DataFrame(plsr.x_weights_, columns=['LV1', 'LV2'], index=prot_data.iloc[:,1])

# Adding IDs to weights
pls_weights['ID'] = prot_data_edit['Protein.Ids']


# Bootstrapping for PLS-DA to find significant cutoffs for weights
prot_data_pls = prot_data.iloc[:,2:-2].replace(np.nan, 0).T
plsr = PLSRegression(n_components=2, scale=True)

pls_boot_means, pls_boot_90_percentile, pls_boot_10_percentile = [], [], []
for _ in range(10000):
    celltype_shuffle = np.random.permutation(celltype)
    plsr.fit(prot_data_pls, celltype_shuffle)
    weights = pd.DataFrame(plsr.x_weights_)[0]
    pls_boot_means.append(weights.mean())
    pls_boot_90_percentile.append(weights.quantile(0.90))
    pls_boot_10_percentile.append(weights.quantile(0.10))


# Plot ordered weights LV1 with bootstrapped cutoffs
plt.figure(figsize=(6, 4), dpi=1200)
plt.scatter(range(len(pls_weights)), pls_weights['LV1'].sort_values(), marker='.', s=1, alpha=0.5, color='k')
plt.axhline(np.mean(pls_boot_90_percentile), color='r')
plt.axhline(np.mean(pls_boot_10_percentile), color='r')
plt.ylabel('Weights LV1', fontweight='bold')
plt.xlabel('Proteins', fontweight='bold')
plt.legend()
plt.tight_layout()
plt.show()

bec_associated = pls_weights.loc[pls_weights['LV1'] <= np.mean(pls_boot_10_percentile), 'LV1']
lec_associated = pls_weights.loc[pls_weights['LV1'] >= np.mean(pls_boot_90_percentile), 'LV1']


## BOOTSTRAPPING TO FIND THE AVERAGE NUMBER OF PROTEINS THAT ARE IN THE TOP AND BOTTOM 10%
# Prepare the data and PLS model
prot_data_pls = prot_data.iloc[:, 2:-2].replace(np.nan, 0).T
celltype_shuffle = celltype.copy()
plsr = PLSRegression(n_components=2, scale=True)

# Define thresholds
pls_top = np.mean(pls_boot_90_percentile)
pls_bottom = np.mean(pls_boot_10_percentile)

# Initialize lists to store results
pls_boot_number_top = []
pls_boot_number_bottom = []
pls_boot_number_combined = []

# Bootstrapping loop
for _ in range(1000):
    np.random.shuffle(celltype_shuffle)
    plsr.fit(prot_data_pls, celltype_shuffle)
    
    weights_lv1 = pd.DataFrame(plsr.x_weights_).iloc[:, 0]
    
    number_top = (weights_lv1 >= pls_top).sum()
    number_bottom = (weights_lv1 <= pls_bottom).sum()
    
    pls_boot_number_top.append(number_top)
    pls_boot_number_bottom.append(number_bottom)
    pls_boot_number_combined.append(number_top + number_bottom)

# Convert results to DataFrame for easier analysis
bootstrapping_results = pd.DataFrame({
    'Top 10% Count': pls_boot_number_top,
    'Bottom 10% Count': pls_boot_number_bottom,
    'Combined Count': pls_boot_number_combined
})

# Display some summary statistics
print(bootstrapping_results.describe())

# Optionally plot the results
figure(figsize=(6, 4), dpi=500)
matplotlib.rcParams.update({'font.size': 12})
sns.kdeplot(pls_boot_number_combined)
sns.kdeplot(pls_boot_number_top)
sns.kdeplot(pls_boot_number_bottom)
plt.axvline((bec_associated.count()+lec_associated.count()), color='tab:blue')
plt.axvline(lec_associated.count(), color='tab:orange')
plt.axvline(bec_associated.count(), color='tab:green')
plt.legend(['Total','Top','Bottom'],loc='center left',bbox_to_anchor=(1,0.5))
#plt.title('Permutation test: number of proteins passing LV1 threshold', fontweight='bold')
plt.ylabel('Density', fontweight='bold')
plt.xlabel('Number of proteins', fontweight='bold')
plt.show()

## CALCULATE P-VALUES WITH BOOTSTRAPPING
def boot_pvalue(distribution, observed_value):
    distribution = np.array(distribution)
    
    # Perform a bootstrapped test
    n_bootstraps = 1000
    bootstrap_means = []
    
    for _ in range(n_bootstraps):
        resampled_data = np.random.choice(distribution, size=len(distribution), replace=True)
        bootstrap_means.append(np.mean(resampled_data))
    
    if observed_value <= np.mean(np.abs(np.array(bootstrap_means))):
        p_value = np.mean(np.abs(np.array(bootstrap_means) <= observed_value))
    
    elif observed_value >= np.mean(np.abs(np.array(bootstrap_means))):
        p_value = np.mean(np.abs(np.array(bootstrap_means) >= observed_value))
    
    print("P-value:", p_value)

boot_pvalue(pls_boot_number_combined, (bec_associated.count()+lec_associated.count()))
boot_pvalue(pls_boot_number_top, lec_associated.count())
boot_pvalue(pls_boot_number_bottom, bec_associated.count())

#%% Calculation of log2(FC) and adj. p-values between days 5/7 and day 2
from multipy.fdr import qvalue

def diff_analysis(data, cellline, day):
    """
    Calculate log2FC, p-values and q-values for each protein between samples.

    Parameters:
    data: proteomics data file with protein intensities as rows and samples as columns.
    cellline: which cellline to perform differential analysis with.
    day: which days to compare.

    """
    # Filter indices for the given cell line and day
    cl_indices = metadata.index[metadata['samples'].str.startswith(cellline)]
    day_indices = cl_indices[metadata.loc[cl_indices, 1].str.contains(str(day))]
    
    # Subset data for Day 2 and the target day
    D2_df = data.iloc[:, cl_indices[0:3]]
    Dx_df = data.iloc[:, day_indices]
    
    # Calculate log2(FC)
    means_D2 = D2_df.mean(axis=1)
    means_Dx = Dx_df.mean(axis=1)
    log2_FC = np.log2(means_Dx / means_D2)
    
    # Calculate p-values
    p_values = [stats.ttest_ind(D2_df.iloc[i, :], Dx_df.iloc[i, :])[1] for i in range(D2_df.shape[0])]
    
    # Adjust p-values using q-value
    q_values = qvalue(pvals=p_values)[1]
    q_values_neg_log10 = -np.log10(q_values)
    
    # Create a DataFrame for differential analysis results
    diff_analysis_df = pd.DataFrame({
        'log2FC': log2_FC,
        'pval': p_values,
        'q-value': q_values,
        'neg log10(q-value)': q_values_neg_log10,
        'Gene names': data.index
    }).reset_index(drop=True)
    
    # Mark significant results based on log2(FC) and q-value thresholds
    diff_analysis_df['significant'] = (np.abs(diff_analysis_df['log2FC']) >= 0.5) & (diff_analysis_df['q-value'] <= 0.05)
    
    # Merge with gene/protein information
    diff_analysis_df = diff_analysis_df.merge(genes_prot_data, on='Gene names')
    
    return diff_analysis_df

# Perform differential analysis for each cell line and day combination
prot_data_edit = prot_data.set_index(prot_data.iloc[:, 1])
diff_results = {}
for cellline in ['HDBEC', 'HDLEC', 'HUVEC', 'iLEC']:
    for day in [5, 7]:
        diff_results[f'{cellline}_{day}D'] = diff_analysis(prot_data_edit.iloc[:, 2:-2], cellline, day)

# Compile log2(FC) and q-value results into a single DataFrame
all_proteins_5D7D_FC = pd.DataFrame({'Gene names': diff_results['HUVEC_7D']['Gene names'],
                                     'Protein description': diff_results['HUVEC_7D']['Protein description'],
                                     'Protein IDs': diff_results['HUVEC_7D']['Protein IDs']})
for key, df in diff_results.items():
    all_proteins_5D7D_FC[f'{key} log2FC'] = df['log2FC']
    all_proteins_5D7D_FC[f'{key} q-value'] = df['q-value']

all_proteins_5D7D_FC.set_index('Gene names', inplace=True)

