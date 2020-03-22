import gseapy
import pandas as pd
from gprofiler import GProfiler


metabric_data = pd.read_csv('../data/metabric.csv.gz', low_memory=False)
print('Success!')
expression_data = metabric_data.iloc[27:,:]
print(expression_data)

new_clinical_patient = pd.read_csv(
    '../data/data_clinical_patient.txt', sep='\t', index_col=0).iloc[4:]
new_clinical_sample = pd.read_csv('../data/data_clinical_sample.txt',
                                  sep='\t', index_col=0).iloc[4:]
new_clinical = pd.concat(
    [new_clinical_patient, new_clinical_sample.reindex(new_clinical_patient.index)], axis=1)
real_gene_names = pd.read_csv('../data/illumina2symbol.txt', sep="\t", index_col=0)

dtypedict = {}
for i in expression_data.columns[1:]:
    dtypedict[i] = 'str'
expression_data = expression_data.astype(dtypedict)

genes = expression_data['Unnamed: 0'].values.tolist()

gp = GProfiler(return_dataframe=True)
gp = gp.convert(organism='hsapiens',query=genes)

gp = gp.loc[gp['n_converted'] == 1]
gp = gp.loc[gp['name'] != 'None']
gp = gp.set_index('incoming')
gprofiler_names = gp


expression_data.index = expression_data['Unnamed: 0']
expression_data = expression_data.iloc[:, 1:]
expression_data = expression_data.rename(
    index=lambda x: gprofiler_names.loc[x, 'name'] if x in gprofiler_names.index else x)
expression_data = expression_data.rename(
    index=lambda x: real_gene_names.loc[x, 'symbol'] if x in real_gene_names.index else x)


gseapatients = []
for i in expression_data.columns.tolist():
    if i in new_clinical.index:
        gseapatients.append(i)
new_clinical_patients = new_clinical.loc[gseapatients]
new_clinical_patients = new_clinical_patients.iloc[:,9]
gsea_list = new_clinical_patients.tolist()



def GSEA(clust):
    classes = []
    for i in gsea_list:
        if i == clust:
            classes.append("cluster")
        else:
            classes.append("no")
    #classes = ["cluster" if x == clust else "no" for x in new_clinical_patients]
    print(classes)
    destination_path = f'GSEA_/{clust}/'
    print(expression_data)
    gs_res = gseapy.gsea(data=expression_data,
                         gene_sets='Reactome_2016',
                         method="signal_to_noise",
                         permutation_num=1000,
                         max_size=10000,
                         min_size=2,
                         processes=6,
                         cls=classes,
                         verbose=True,
                         outdir=destination_path)

clusterlist = ['1', '2', '3', '4ER+', '4ER-', '5', '6', '7', '8', '9', '10']
for i in clusterlist:
    print(i)
    GSEA(i)
