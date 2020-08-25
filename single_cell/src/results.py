import pandas as pd
import numpy as np

#importing the dataset and the pathway name information
df = pd.read_csv("../data/idea_pathways.csv", index_col = 0)
pathway_names = pd.read_csv("../../metabric/data/ReactomePathwaysHuman.txt", sep="\t", index_col = 0, names=["name", "species"])

#format pathway names to match
new_index = [i.replace(".", "-") for i in df.index.tolist()]

#reindex the dataframe
df["index"] = new_index
df = df.set_index("index")

#change df index to the real names instead of reactome notation
df["full_name"] = [pathway_names.loc[i, "name"] for i in new_index]
df = df.set_index("full_name")

#new column with -log10(pvalue) to use for sunburst
df["-log10(pvalue)"] = -np.log10(df["pvalue_louis"])

#select the -log10(pvalues) and write them to a csv
sunburst_df = df["-log10(pvalue)"]
sunburst_df = sunburst_df.apply(lambda x: 100 if x == np.Inf else x)
sunburst_df.to_csv("../exp/iDEA.csv", sep = "\t")
