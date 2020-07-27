import gseapy
import pandas as pd


expression_data = pd.read_csv('../data/GSEA_df.csv', index_col=0)
sample_data = pd.read_csv('../data/GSEA_samples.csv', index_col=0)
#gene expression dataset and clinical_data
print(expression_data.iloc[:10])

samples = [item for item in sample_data.index.tolist() if item in expression_data.columns.tolist()]
#match patient names
print(samples)

expression_data = expression_data.T.loc[samples].T
#select expression data based on which patients are in clinical_data


def GSEA(sample):
    classes = ['MUT' if x == "DRESS" else 'WT' for x in sample_data[sample]] #define levels
    destination_path = 'GSEA_/single_cell/' #out dir for the GSEA file

    gs_res = gseapy.gsea(data=expression_data,
                         gene_sets='Reactome_2016',
                         method="signal_to_noise",
                         permutation_num=10000,
                         max_size=100000,
                         min_size=1,
                         cls=classes,
                         no_plot = True,
                         #processes = 4, #multiprocessing
                         verbose=True,
                         outdir=destination_path)

GSEA('sample')
