import numpy as np
import scipy
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt

def cell_calling(gene_count_matrix,barcodes,genes,U,T,num_sim=1000,seed=42):
    total_transcript_per_cell = np.array(gene_count_matrix.sum(axis=0)).flatten()
    ordered_cell_idx = np.argsort(total_transcript_per_cell)[::-1]
    sorted_gene_count_matrix = gene_count_matrix[:,ordered_cell_idx]
    gene_count = sorted_gene_count_matrix.T
    barcodes = [barcodes[i] for i in ordered_cell_idx]
     # filter genes without count (in ambient barcodes)
    ambient_count = gene_count[T:] 
    total_gene_count_ambient = ambient_count.sum(axis=0)
    nonzero_ind = np.flatnonzero(total_gene_count_ambient)

    # break gene_count matrix into the three count matrices
    real_count = gene_count[:U, nonzero_ind]
    putative_count = gene_count[U:T, nonzero_ind]
    ambient_count = gene_count[T:, nonzero_ind] 
    genes = [genes[i] for i in nonzero_ind]

    # get normalized count matrix 
    ambient_freq = ambient_count / np.array(ambient_count).sum(axis=1, keepdims=True)
    ambient_mean = np.array(ambient_freq.mean(axis=0)).flatten()
    ambient_var = np.array(ambient_freq.var(axis=0)).flatten()

    res = sm.OLS(ambient_var, ambient_mean).fit()
    slope = res.params.item()

    lower_cutoff = int(np.ceil(1/slope))
    '''Both real and putative counts should have per-barcode total count greater than lower_cutoff'''
    ambient_ind_to_keep = np.flatnonzero(ambient_count.sum(axis=1) > lower_cutoff)
    ambient_count = ambient_count[ambient_ind_to_keep]

    putative_pvals = _collect_pvalues(putative_count, ambient_mean, slope, num_sim, seed=seed)

    '''trt setting this error rate, alpha, as low as possible (e.g., 0.002)'''
    putative_reject, putative_pvals_corrected = multipletests(putative_pvals, alpha=0.05, method='fdr_bh')[:2]  
    gene_count_called = np.vstack((real_count, putative_count[putative_reject]))
    putative_idx = [putative_reject[i]*i for i in range(T-U) if putative_reject[i]]
    barcodes = barcodes[:U] + [barcodes[U+i] for i in putative_idx]
    return gene_count_called.T, barcodes, genes

def _collect_pvalues(count_matrix, prior_p, slope, num_sim=1000, seed=42):
    pvals_total = []
    for barcode in count_matrix:
        total_count = barcode.sum()

        mean = prior_p * total_count
        var = (prior_p * slope) * (total_count ** 2) 
        
        phi = np.array([max(i,1e-2) for i in mean**2 / (var - mean)])
        p = mean / var

        rng = np.random.default_rng(seed)

        Lb = scipy.stats.nbinom.logpmf(barcode, phi, p).sum()
        samples = scipy.stats.nbinom.rvs(phi, p, size=(num_sim, len(prior_p)), random_state=rng)
        Lsamples = scipy.stats.nbinom.logpmf(samples, phi, p).sum(axis=1)
        pval = (np.count_nonzero(Lsamples < Lb) + 1) / (num_sim + 1)
        pvals_total.append(pval)   

    return np.array(pvals_total)


def find_UT_by_quantile(gene_count_matrix,U_quantile=0.5,T_quantile=0.1):
    total_transcript_per_cell = np.array(gene_count_matrix.sum(axis=0)).flatten()
    U_threshold = np.quantile(total_transcript_per_cell,U_quantile)
    T_threshold = np.quantile(total_transcript_per_cell,T_quantile)   
    U = sum(total_transcript_per_cell >= U_threshold)
    T = sum(total_transcript_per_cell >= T_threshold)
    return U,T

def plot_loglog(gene_count_matrix,U,T):
    num_genes,num_cells = gene_count_matrix.shape
    total_transcript_per_cell = np.array(gene_count_matrix.sum(axis=0)).flatten()
    
    ordered_cell_idx = np.argsort(total_transcript_per_cell)[::-1]
    sorted_gene_count_matrix = gene_count_matrix[:,ordered_cell_idx]
    sorted_total_transcript_per_cell = total_transcript_per_cell[ordered_cell_idx]
    U_threshold = sorted_total_transcript_per_cell[U-1]
    T_threshold = sorted_total_transcript_per_cell[T-1]
        
    plt.loglog([i for i in range(num_cells)],sorted_total_transcript_per_cell)
    plt.axhline(y = U_threshold, color = 'g', linestyle = '--', label = "U")
    plt.axhline(y = T_threshold, color = 'r', linestyle = '--', label = "T")
    plt.text(1, 10 ** ((np.log10(T_threshold) + np.log10(min(total_transcript_per_cell))) / 2),\
             'Ambient Cells', dict(size=12))
    plt.text(1, 10 ** ((np.log10(T_threshold) + np.log10(U_threshold)) / 2),\
         'Putative Cells', dict(size=12))
    plt.text(1, 10 ** ((np.log10(max(total_transcript_per_cell)) + np.log10(U_threshold)) / 2),\
         'Real Cells', dict(size=12))
    plt.legend(bbox_to_anchor = (1.0, 1), loc = 'upper left')
    plt.xlabel('Rank of cell')
    plt.ylabel('Number of transcripts per cell')
    
    