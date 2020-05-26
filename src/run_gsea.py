import gseapy
import pandas as pd


x = input("Enter cluster:")


metabric_data = pd.read_csv('../data/GSEA_expression.csv', index_col=0)
clinical_data = pd.read_csv('../data/GSEA_clinical.csv', index_col=0)
#gene expression dataset and clinical_data

clinical_patients = clinical_data.index.tolist() # patient names as list

clinical_patients = [item for item in clinical_data.index.tolist() if item in metabric_data.columns.tolist()]
#match patient names


metabric_data = metabric_data.T.loc[clinical_patients].T
#select expression data based on which patients are in clinical_data


def GSEA(receptor, level):
    classes = ['MUT' if x == level else 'WT' for x in clinical_data[receptor]] #define levels
    destination_path = f'GSEA_/adj_{level}/' #out dir for the GSEA file

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



GSEA('Integrative Cluster', x)
