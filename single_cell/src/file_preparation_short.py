import urllib.request
import collections
import scipy.sparse as sp_sparse
import tables
import numpy as np
import pandas as pd

file_urls = {"Dress1": "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3892569&format=file&file=GSM3892569%5FSkin%5FDress1%5Ffiltered%5Fgene%5Fbc%5Fmatrices%5Fh5%2Eh5",
             "hv1": "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3892577&format=file&file=GSM3892577%5FSKIN%5FHV1%5FF1%5Ffiltered%5Fgene%5Fbc%5Fmatrices%5Fh5%2Eh5",
             "hv2": "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3892579&format=file&file=GSM3892579%5FSKIN%5FHV2%5Ffiltered%5Fgene%5Fbc%5Fmatrices%5Fh5%2Eh5",
             "hv3": "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3892580&format=file&file=GSM3892580%5FSKIN%5FHV3%5Ffiltered%5Fgene%5Fbc%5Fmatrices%5Fh5%2Eh5",
             "hv4": "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3892581&format=file&file=GSM3892581%5FSKIN%5FHV4%5Ffiltered%5Fgene%5Fbc%5Fmatrices%5Fh5%2Eh5",
             "hv5": "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3892582&format=file&file=GSM3892582%5FSKIN%5FHV5%5Ffiltered%5Fgene%5Fbc%5Fmatrices%5Fh5%2Eh5"}


def download_files(url, filename):
    file = filename + ".h5"
    path = "../data/" + file
    print(f"Downloading {file} to data folder")
    urllib.request.urlretrieve(url, path)


def get_matrix_from_h5(filename):
    with tables.open_file(filename, 'r') as f:
        mat_group = f.get_node(f.root, 'GRCh38')
        barcodes = f.get_node(mat_group, 'barcodes').read()
        data = getattr(mat_group, 'data').read()
        indices = getattr(mat_group, 'indices').read()
        indptr = getattr(mat_group, 'indptr').read()
        shape = getattr(mat_group, 'shape').read()
        matrix = sp_sparse.csc_matrix((data, indices, indptr), shape=shape)
        feature_ref = {}
        genes = f.get_node(mat_group, 'genes')
        gene_names = f.get_node(mat_group, 'gene_names')
        gene_dict = dict(zip(genes, gene_names))
        return gene_dict, matrix

def make_df(file, matrix):
    df_index = [x.decode("UTF-8") for x in gene_dict.keys()]
    df = pd.DataFrame(data=matrix.toarray(),
                      index = df_index,
                      columns = [f"{file}_{x}" for x in range(1, matrix.toarray().shape[1] +1)])
    return df


dataframes = {}
for i in file_urls.keys():
    download_files(file_urls[i], i)
    file = i + ".h5"
    path = "../data/" + file
    gene_dict, matrix = get_matrix_from_h5(path)
    df = make_df(i, matrix)
    dataframes[i] = df



full_df = pd.concat([dress1, hv1.reindex(dress1.index)], axis = 1)
full_df = pd.concat([full_df, dataframes["hv2"].reindex(full_df.index)], axis = 1)
full_df = pd.concat([full_df, dataframes["hv3"].reindex(full_df.index)], axis = 1)
full_df = pd.concat([full_df, dataframes["hv4"].reindex(full_df.index)], axis = 1)
full_df = pd.concat([full_df, dataframes["hv5"].reindex(full_df.index)], axis = 1)
full_df = full_df.fillna(0)
print("Expression data saved to ../data/full_df.csv")
full_df.to_csv("../data/full_df.csv")

phenotype_df = pd.DataFrame(index=full_df.columns)
phenotype_df["phenotype"] = ["DRESS" if x.startswith("Dress1") else "HV" for x in full_df.columns]
phenotype_df_t = phenotype_df.T
print("Phenotype data saved to ../data/pheno_df.csv")
phenotype_df.to_csv("../data/pheno_df.csv")
