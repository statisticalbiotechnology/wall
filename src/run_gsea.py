import gseapy
import pandas as pd
from gprofiler import GProfiler


metabric_data = pd.read_csv('../data/metabric.csv.gz')
expression_data = metabric_data.iloc[27:,:]
new_clinical_patient = pd.read_csv('../data/data_clinical_patient.txt', sep='\t', index_col=0).iloc[4:]
new_clinical_sample = pd.read_csv('../data/data_clinical_sample.txt', sep='\t', index_col=0).iloc[4:]
new_clinical = pd.concat([new_clinical_patient, new_clinical_sample.reindex(new_clinical_patient.index)], axis=1)
real_gene_names = pd.read_csv('../data/illumina2symbol.txt', sep="\t", index_col = 0)
'''
new_clinical['Triple Neg'] = new_clinical.apply(lambda row: True if ((row['ER Status'] == 'Negative')
                                                                     and (row['PR Status'] == 'Negative')
                                                                     and (row['HER2 Status'] == 'Negative')) else False, axis = 1)

new_clinical['ER-/PR-/HER2+'] = new_clinical.apply(lambda row: True if ((row['ER Status'] == 'Negative')
                                                                     and (row['PR Status'] == 'Negative')
                                                                     and (row['HER2 Status'] == 'Positive')) else False, axis = 1)

'''
dtypedict = {}
for i in expression_data.columns[1:]:
    dtypedict[i] = 'float32'
expression_data = expression_data.astype(dtypedict)

genes = expression_data['Unnamed: 0'].values.tolist()

gp = GProfiler(return_dataframe = True)
gp = gp.convert(organism='hsapiens',
          query=genes)

gp = gp.loc[gp['n_converted'] == 1]
gp = gp.loc[gp['name'] != 'None']
gp = gp.set_index('incoming')
gprofiler_names = gp

dtypedict = {}
for i in expression_data.columns[1:]:
    dtypedict[i] = 'float32'
expression_data = expression_data.astype(dtypedict)

expression_data.index = expression_data['Unnamed: 0']
expression_data = expression_data.iloc[:,1:]
expression_data = expression_data.rename(index= lambda x: gprofiler_names.loc[x, 'name'] if x in gprofiler_names.index else x)
expression_data = expression_data.rename(index= lambda x: real_gene_names.loc[x, 'symbol'] if x in real_gene_names.index else x)


gseapatients = expression_data.columns.tolist()
new_clinical_patients = new_clinical.loc[gseapatients]
classes = [x for x in new_clinical_patients['Integrative Cluster']]

def GSEA(cluster):
    gseapatients = expression_data.columns.tolist()
    new_clinical_patients = new_clinical.loc[gseapatients]
    classes = ["cluster" if x == cluster else "no" for x in new_clinical_patients['Integrative Cluster']]
    destination_path = 'GSEA_/' + cluster

    gs_res = gseapy.gsea(data = expression_data,
                    gene_sets ='Reactome_2016',
                    method = "t_test",
                    permutation_num = 1000,
                    max_size = 10000,
                    min_size = 2
                    no_plot = True,
                    processes = 4,
                    cls = classes,
                    verbose = True,
                    outdir = destination_path)

for cluster in new_clinical["Integrative Cluster"].unique().tolist():
    GSEA(cluster)
