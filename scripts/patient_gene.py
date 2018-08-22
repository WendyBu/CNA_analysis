"""
for patient_gene heatmap
x: individual patient
y: genelist in each chrom
"""

import pandas as pd
import numpy as np
import glob
import os
import os.path


def order_table(df, chr):
    chr_df = pd.read_csv("../chrom_gene_list/"+chr+"_gene.csv", sep="\t", index_col=4)
    merged_df = chr_df.join(df, how="right")
    merged_df.sort_values("start", inplace=True)
    return merged_df


def clean_df(df):
    vlist = [-1, -2, 0, 1, 2]
    test = df.iloc[3, :]  # randomly picked the third row
    test_result = test[~test.isin(vlist)].index
    df.drop(test_result, axis=1, inplace=True)
    if ("chrom" in df.columns):
        df.drop("chrom", axis=1, inplace=True)
    return df


def generate_heatmap_table(keyword="Melanoma", chr="chr8"):
    main_df = pd.DataFrame()
    input_path = os.path.join("..", keyword+"_temp", chr, "part_CNA")
    for file in glob.glob(input_path+"/*"):
        df_cna = pd.read_csv(file, sep="\t",index_col=0)
        if "Cytoband" in df_cna.columns:
            df_cna.drop("Cytoband", axis=1, inplace=True)
        if not main_df.empty:
            main_df = main_df.join(df_cna, lsuffix='', rsuffix='_y')
        else:
            main_df = df_cna
        # print main_df.shape
    df_ordered = order_table(main_df, chr)
    df_ordered = df_ordered.replace(r'\s+', np.nan, regex=True)
    df_ordered.fillna(0, inplace=True)
    df_clean = clean_df(df_ordered)
    output_dir = os.path.join("..", keyword+"_temp", "reports", "heatmap_"+keyword+"_"+chr+".xls")
    df_clean.to_csv(output_dir,sep="\t")


def main(keyword="Melanoma"):
    chr_list = list(range(1,23))
    chr_list.extend(["X", "Y"])
    chr_list = ["chr"+str(chr) for chr in chr_list]
    for chr in chr_list:
        generate_heatmap_table(keyword, chr)
    pass


# keywordlist = ["Prostate","Breast", "Bladder", "Colon", "GBM"]
keywordlist=["Prostate"]
if __name__ == "__main__":
    for key in keywordlist:
        main(key)

