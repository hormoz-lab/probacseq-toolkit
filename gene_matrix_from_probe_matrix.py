import h5py
import os
import numpy as np
import scipy.sparse as sparse
import scipy.io as sio
from .utils import *


def gene_matrix_from_probe_matrix(data_matrix_h5, probe_to_gene_map, map_method,\
                                  UMI_per_cell_threshold=0,map_function=None):
    if map_method not in ['max', 'median', 'mean', 'custom']:
        raise Exception(f'Unrecorgnized map_method: {map_method}')

    probe_count_matrix, CBs, probe_name = load_10x_h5_matrix(data_matrix_h5)
    total_probe_per_cell = np.array(probe_count_matrix.sum(axis=1)).flatten()
    probe_count_matrix = probe_count_matrix[total_probe_per_cell>UMI_per_cell_threshold,:]
    CBs = CBs[total_probe_per_cell>UMI_per_cell_threshold]
    
    probe_count_matrix = probe_count_matrix.todense()
    N_genes = len(probe_name)
    N_cells = len(CBs)
    print(f'Collapsing probes for {len(CBs)} common cells using method={map_method}')
    
    gene_name = [probe_to_gene_map[probe.decode()] for probe in probe_name]
    gene_name, which_gene = np.unique(gene_name, return_inverse=True)
    gene_idx = [np.argwhere(which_gene == i).flatten()
                for i in range(len(gene_name))]
    
    gene_count_matrix = []
    for idx in gene_idx:
        mat = probe_count_matrix[:, idx]
        if map_method == "max":
            gene_count = np.max(mat,axis=1)
        elif map_method == "median":
            gene_count = np.median(mat,axis=1)
        elif map_method == "mean":
            gene_count = np.mean(mat,axis=1)
        elif map_method == "custom":
            gene_count = map_function(mat)
        gene_count = np.array(gene_count).flatten()
        gene_count_matrix.append(gene_count)
        
    gene_count_matrix = np.array(gene_count_matrix)
    return gene_count_matrix, [b.decode('UTF-8') for b in CBs], gene_name.tolist()
    
          
    
# def gene_matrix_from_probe_matrix(data_matrix_h5, outdir, probe_to_gene_map, map_method, map_function=None):
#     print('\n-----start-------\n')
#     dataset = h5py.File(data_matrix_h5,'r')
#     mtx = dataset['matrix']   
#     genes = list(set(probe_to_gene_map.values()))

#     i = mtx['indices']
#     indices = np.array(i)

#     d = mtx['indptr']
#     indptr = np.array(d)

#     s = mtx['shape']
#     shape = np.array(s)

#     x = mtx['data']
#     y = mtx['features']['id']
#     z = mtx['barcodes']
#     data = np.array(x)
#     features = np.char.decode(np.array(y))
#     barcodes = np.char.decode(np.array(z))

#     start = 0
#     row,col,val = [],[],[]
#     mark=time.time()
#     for c in range(1,len(indptr)-1): #for each cell
#         if c%1000==0:
#             print(c,'cells parsed ......')
#         probes = [features[p] for p in indices[start:indptr[c]]]
#         UMIs = data[start:indptr[c]]

#         cell_dict={}
#         for z in zip(probes,UMIs):
#             if z[0] in probe_to_gene_map:
#                 gene = probe_to_gene_map[z[0]]
#                 if gene in cell_dict:
#                     cell_dict[gene].append(z[1])
#                 else:
#                     cell_dict[gene]=[z[1]]

#         for i in range(len(genes)):
#             gene = genes[i]
#             try:
#                 if len(cell_dict[gene])>0:
#                     row.append(i)
#                     col.append(c-1)
#                     if map_method == "max":
#                         v = max(cell_dict[gene])
#                     elif map_method == "mean":
#                         v = mean(cell_dict[gene])
#                     elif map_method == "median":
#                         v = median(cell_dict[gene])
#                     elif map_method == "custom":
#                         v = map_function(cell_dict[gene])   
#                     val.append(v)
#             except KeyError:
#                 pass
#         start = indptr[c]
    
#     print(c,'cells parsed')
#     print(time.time()-mark,'s')
    
#     # creating sparse matrix
#     gene_count_matrix = sparse.csr_matrix((val, (row, col)), 
#                               shape = (len(genes), len(barcodes)))#.toarray()
    
#     nonzero_barcode = np.flatnonzero(np.sum(gene_count_matrix,axis=0))
#     nonzero_gene = np.flatnonzero(np.sum(gene_count_matrix,axis=1))
#     gene_count_matrix = gene_count_matrix[nonzero_gene,:][:,nonzero_barcode]
#     barcodes = [barcodes[i] for i in nonzero_barcode]
#     genes = [genes[i] for i in nonzero_gene]
#     return gene_count_matrix, barcodes, genes

