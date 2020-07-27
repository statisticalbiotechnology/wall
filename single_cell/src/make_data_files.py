import pandas as pd

files = ["Skin_Dress1", "SKIN_HV1_F1", "SKIN_HV1_F2", "SKIN_HV3", "SKIN_HV4", "SKIN_HV5"]

dress1 = pd.read_csv(f"../data/{files[0]}/{files[0]}.csv", index_col = 0)
hv1 = pd.read_csv(f"../data/{files[1]}/{files[1]}.csv", index_col = 0)
hv2 = pd.read_csv(f"../data/{files[2]}/{files[2]}.csv", index_col = 0)
hv3 = pd.read_csv(f"../data/{files[3]}/{files[3]}.csv", index_col = 0)
hv4 = pd.read_csv(f"../data/{files[4]}/{files[4]}.csv", index_col = 0)
hv5 = pd.read_csv(f"../data/{files[5]}/{files[5]}.csv", index_col = 0)

print("finished importing files")

def make_genecount_df(input_df):
    df = pd.DataFrame()
    for i in input_df["gene_name"].unique().tolist():
        df[i] = pd.Series(input_df.loc[input_df["gene_name"] == i]["data"].values)
    return df


dress1_df = make_genecount_df(dress1)
dress1_df.to_csv("../data/finished_data_files/dress1.csv")
print("finished with dress1 file")

hv1_df = make_genecount_df(hv1)
hv1_df.to_csv("../data/finished_data_files/hv1.csv")
print("finished with hv1 file")


hv2_df = make_genecount_df(hv2)
hv2_df.to_csv("../data/finished_data_files/hv2.csv")
print("finished with hv2 file")


hv3_df = make_genecount_df(hv3)
hv3_df.to_csv("../data/finished_data_files/hv3.csv")
print("finished with hv3 file")


hv4_df = make_genecount_df(hv4)
hv4_df.to_csv("../data/finished_data_files/hv4.csv")
print("finished with hv4 file")


hv5_df = make_genecount_df(hv5)
hv5_df.to_csv("../data/finished_data_files/hv5.csv")
print("finished with hv5 file")
