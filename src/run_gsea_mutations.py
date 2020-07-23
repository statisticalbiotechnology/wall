import gseapy
import pandas as pd
import time


metabric_data = pd.read_csv('../data/GSEA_expression.csv', index_col=0)
mutation_data = pd.read_csv('../data/mutations_per_patients.csv', index_col=0)
#gene expression dataset and clinical_data

#mutation_patients = mutation_data.index.tolist() # patient names as list
#print(mutation_patients)


mutation_patients = [item for item in mutation_data.index.tolist() if item in metabric_data.columns.tolist()]
#match patient names


metabric_data = metabric_data.T.loc[mutation_patients].T
#select expression data based on which patients are in clinical_data


def GSEA(mutation):
    classes = ['MUT' if x == "MUT" else 'WT' for x in mutation_data[mutation]] #define levels
    destination_path = f'GSEA_/{mutation}/' #out dir for the GSEA file

    gs_res = gseapy.gsea(data=metabric_data,
                         gene_sets='../data/GSEA_reactome.gmt',
                         method="signal_to_noise",
                         permutation_num=1000,
                         max_size=100000,
                         min_size=1,
                         cls=classes,
                         no_plot = True,
                         #processes = 4, #multiprocessing
                         verbose=True,
                         outdir=destination_path)



GSEA('BRCA2')
