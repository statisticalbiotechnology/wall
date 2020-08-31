import pandas as pd
from gprofiler import GProfiler
from ast import literal_eval

def read_reactome(file_name, gene_name_start = "ENSG0"):
    df = pd.read_csv(file_name, sep='\t', header=None)

    if gene_name_start == None:
        sub_df = df
    else:
        subset_vec = df[0].str.startswith(gene_name_start)
        sub_df = df.loc[subset_vec]

    genes_df = sub_df.groupby(1)[0].apply(list)
    names_df = sub_df.groupby(1)[3].max()

    out_df = pd.concat([genes_df,names_df], axis=1)
    out_df.columns = ['genes', 'pathway_name']

    return out_df



df = pd.read_csv("../data/lfc_df.csv", index_col=0)
genes = df.index.tolist()

gp = GProfiler(return_dataframe = True)
gp = gp.convert(organism='hsapiens',
          query=genes)

gp = gp.loc[gp['n_converted'] == 1]
gp = gp.loc[gp['name'] != 'None']
gp = gp.set_index('incoming')
gprofiler_names = gp
gprofiler_names = gprofiler_names.set_index("converted")
gprofiler_names

converted_genes = gprofiler_names["name"]
converted_genes = converted_genes.rename("symbol")

react_ome = read_reactome('../../metabric/data/ensembl2reactome_all_Levels.txt.gz')
final_df = pd.DataFrame(index=gprofiler_names["name"])
final_df = pd.DataFrame(index=df.index)

for pathway in react_ome.index:
    genes = react_ome.loc[pathway, "genes"]
    converted_list = []
    for i in genes:
        if i in gprofiler_names.index:
            selected_gene = gprofiler_names.loc[i, "name"]
            if type(selected_gene) == str:
                converted_list.append(selected_gene)
    bin_list = []
    for i in final_df.index:
        if i in converted_list:
            bin_list.append(1)
        else:
            bin_list.append(0)
    final_df[pathway] = bin_list

final_df.to_csv("../data/annotation_reactome.csv")
