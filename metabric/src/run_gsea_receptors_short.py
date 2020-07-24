import gseapy
import pandas as pd


metabric_data = pd.read_csv('../data/expression_data.csv', index_col=0)
clinical_data = pd.read_csv('../data/GSEA_clinical.csv', index_col=0)
clinical_patients = clinical_data.index.tolist()

metabric_data = metabric_data.T.loc[clinical_patients].T



def GSEA(receptor, level):
    classes = ['MUT' if x == level else 'WT' for x in clinical_data[receptor]]
    destination_path = f'GSEA_/{receptor}/'

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
GSEA('Triple Neg', True)
GSEA('ER-/PR-/HER2+', True)
