import pandas as pd
import glob
import os
import os.path
import numpy as np


vlist = [-1, -2, 0, 1, 2]

df = pd.read_csv("../Prostate_temp/reports/heatmap_Prostate_chr8.xls", sep="\t", index_col=0)
print df.shape
test = df.iloc[3,:]
test_result = test[~test.isin(vlist)].index
df.drop(test_result,axis=1,inplace=True)
print df.shape


# print np.sum(df.columns.str.contains("-"))
# print np.sum(~df.columns.str.contains("-"))
# df.drop("Entrez_Gene_Id", axis=1, inplace=True)

# vlist = [-1, -2, 0, 1, 2]
# for i in range(0, 5):
#     test = df.iloc[i, :]
#     drop_col = test[~test.isin(vlist)].index
#
#     df.drop(labels=drop_col, axis=1, inplace=True)
# Entrez_Gene_Id

# print df.shape

# print file.head()
# print file.shape


# print file.shape
# for col in file:
#     print file[file[col].isin(vlist)]
# print file.index.values




# file.to_csv("../Prostate_temp/reports/chr8_test.xls", sep="\t")