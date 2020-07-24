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

new_clinical['Triple Neg'] = new_clinical.apply(lambda row: 'Triple Neg' if ((row['ER Status'] == 'Negative')
                                                                     and (row['PR Status'] == 'Negative')
                                                                     and (row['HER2 Status'] == 'Negative')) else 'Not Triple Neg', axis = 1)

new_clinical['ER-/PR-/HER2+'] = new_clinical.apply(lambda row: 'ER-PR-HER2+' if ((row['ER Status'] == 'Negative')
                                                                     and (row['PR Status'] == 'Negative')
                                                                     and (row['HER2 Status'] == 'Positive')) else 'not ER-PR-HER2+', axis = 1)






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

new_clinical.to_csv('../data/gsea_clinical.csv')
expression_data.to_csv('../data/expression_data.csv')
