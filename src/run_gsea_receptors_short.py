import gseapy
import pandas as pd
import numpy as np


metabric_data = pd.read_csv('../data/expression_data_metabric.csv', index_col=0)
new_clinical_patients = pd.read_csv('../data/GSEA_clinical.csv')


def GSEA(receptor, level):
    classes = ['MUT' if x == level else 'WT' for x in new_clinical_patients[receptor]]
    destination_path = f'GSEA_/{receptor}/'
    print(new_clinical_patients[receptor])
    print(classes)
    gs_res = gseapy.gsea(data=metabric_data,
                         gene_sets='Reactome_2016',
                         method="signal_to_noise",
                         permutation_num=1000,
                         max_size=10000,
                         min_size=2,
                         processes=6,
                         cls=classes,
                         verbose=True,
                         outdir=destination_path)


GSEA('ER Status', 'Positive')
GSEA('PR Status', 'Positive')
GSEA('HER2 Status', 'Positive')
