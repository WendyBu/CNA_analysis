"""
organize all the TCGA gistic results
copy number variable in each arms of chromosome
"""


import pandas as pd
import glob, re
pd.set_option("display.max_columns", 100)


df_ampFreq = pd.DataFrame()
df_delFreq = pd.DataFrame()
df_ampQ = pd.DataFrame()
df_delQ = pd.DataFrame()
df_ampZ = pd.DataFrame()
df_delZ = pd.DataFrame()
cnv_path = "/Users/yiwenbu/Desktop/CNV_TCGA_GISTIC_RESULTS/*/broad_significance_results.txt"

for file in glob.glob(cnv_path):
    try:
        tumorName = re.search('org_(.+?)-T', file).group(1)
    except AttributeError:
        print file
    df = pd.read_csv(file, sep="\t", index_col=0)
    df_ampFreq[tumorName] = df["Amp frequency"]
    df_delFreq[tumorName] = df["Del Frequency"]
    df_ampQ[tumorName] = df["Amp q-value"]
    df_delQ[tumorName] = df["Del q-value"]
    df_ampZ[tumorName] = df["Amp z-score"]
    df_delZ[tumorName] = df["Del z-score"]


df_ampFreq.to_csv("forms/amp_frequency.xls", sep="\t")
df_delFreq.to_csv("forms/del_frequency.xls", sep="\t")
df_ampQ.to_csv("forms/amp_q_value.xls", sep="\t")
df_delQ.to_csv("forms/del_q_value.xls", sep="\t")
df_ampZ.to_csv("forms/amp_z_score.xls", sep="\t")
df_delZ.to_csv("forms/del_z_score.xls", sep="\t")














