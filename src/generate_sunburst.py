import pandas as pd
import numpy as np
import networkx as nx
from networkx.readwrite import json_graph
import json

relation_file = "../data/ReactomePathwaysRelation.txt"
pathway_name = "../data/ReactomePathways.txt"

def generate_tree(relation_file = relation_file):
    rel_df = pd.read_csv(relation_file, sep = "\t", header = None, index_col = 0, names=['id'])
    real_names = pd.read_csv(pathway_name, sep='\t', index_col=0, names=['pathway', 'name', 'organism'])


    cut = rel_df.index.str.contains('MMU') & rel_df['id'].str.contains('MMU')
    rel_df = rel_df.loc[cut]

    namelist = []
    for i in rel_df.index:
        conversion = real_names.loc[i, 'name']
        namelist.append(conversion)

    childlist = []
    for i in rel_df['id']:
        conversion = real_names.loc[i, 'name']
        childlist.append(conversion)
    names = pd.DataFrame()
    names['parent'] = namelist
    names['child'] = childlist


    G = nx.DiGraph()
    G.add_edges_from(names.values)
    roots = [n for n,d in G.in_degree() if d==0]
    print(roots)
    print(type(roots))

    roots_df = pd.DataFrame(columns = ['parentId', 'id'])
    print(roots_df)

    roots_df['id'] = roots
    roots_df['parentId'] = 'Human'

    roots_df = pd.DataFrame(roots_df.values, columns = ['parentId', 'id'])
    rel_df = pd.DataFrame(names.values, columns = ['parentId', 'id'])

    tree = roots_df.append(rel_df)
    return tree

rel_df = generate_tree()

def default(o):
     if isinstance(o, np.integer): return int(o)
     raise TypeError

def sunburst(in_df, outname = 'sun_tree.json'):
    max_val = in_df['value'].max()
    highest_rank = float(len(in_df.index)+1)
    topPaths = rel_df.loc[(rel_df['parentId'] == 'Human'), 'id']
    homoNgenes = np.sum(in_df.loc[[x in topPaths.tolist() for x in in_df.index],'ngenes'])
    homoNode = pd.DataFrame([[0,homoNgenes,"Human", max_val, highest_rank]], columns = ["value", "ngenes", "Organism", "max_val", 'rank']).xs(0)
    homoNode.name = 'Human'

    in_df = in_df.append(homoNode)
    #print(in_df)
    topDict = in_df.to_dict()

    pathways = in_df.index

    n_path = len(pathways)

    subset_vec = [x in pathways for x in rel_df.iloc[:,0]] and [x in pathways for x in rel_df.iloc[:,1]]
    sub_rel_df = rel_df[subset_vec]

    G = nx.DiGraph()

    G.add_nodes_from(pathways)
    G.add_edges_from(sub_rel_df.values)

    tree = nx.algorithms.dag.dag_to_branching(G)

    secondDict = nx.get_node_attributes(tree,'source')
    #print(secondDict)

    thirdDict = {'value':{}, 'ngenes':{}, 'max_val': {}, 'rank': {}}
    for key, value in secondDict.items():
        thirdDict['value'].update({key : topDict['value'][value]})
        thirdDict['ngenes'].update({key : topDict['ngenes'][value]})
        thirdDict['max_val'].update({key : topDict['max_val'][value]})
        thirdDict['rank'].update({key : topDict['rank'][value]})

    nx.set_node_attributes(tree, thirdDict['value'], name = 'value')
    nx.set_node_attributes(tree, thirdDict['ngenes'], name = 'ngenes')
    nx.set_node_attributes(tree, thirdDict['max_val'], name = 'max_val')
    nx.set_node_attributes(tree, thirdDict['rank'], name = 'rank')

    root = [v for v, d in tree.in_degree() if d == 0][0]
    out_json = json_graph.tree_data(tree, root)

    with open(outname, 'w') as outfile:
        print(outname)
        json.dump(out_json, outfile, default=default)


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
    out_df.index = out_df.pathway_name

    return out_df




def make_the_json_files():
    cluster_df = pd.read_csv("../exp/GSEA_qvalues.csv", index_col = 0)
    #clusterindex = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21]
    #cluster_df = cluster_df.iloc[:,[1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21]]
    reactome_ngenes = read_reactome("../data/Ensembl2Reactome_All_Levels.txt.gz")

    length_dict = {}
    for i in cluster_df.index:
        if i in reactome_ngenes.index:
            nr_genes = len(reactome_ngenes.loc[i, "genes"])
        else:
            print('shit')
        length_dict[i] = nr_genes

    cluster_df['ngenes'] = cluster_df.index.map(length_dict)


    df_dict = {}

    for i in cluster_df.iloc[:,:-1]:
        df = pd.DataFrame(index = cluster_df.index)
        df['value'] = cluster_df.loc[:,i]
        df['ngenes'] = cluster_df.loc[:,'ngenes']
        df['Organism'] = 'Human'
        df['max_val'] = round(df['value'].max(),3)
        sorted_in_df = df.sort_values(by='value', ascending = False)
        sorted_in_df['rank'] = df.reset_index().index +1
        df = sorted_in_df
        print(df)
        df_dict[i] = df

    for i in df_dict:
        clust = i.strip('cluster qva')
        sunburst(df_dict[i], outname = f'sunburst/GSEA_clust_{clust}.json')

make_the_json_files()
