import pandas as pd
import sys

df_all = pd.read_table("all_bioportal_list.txt", sep="\t", index_col=0)
df = df_all[df_all["CNA"]>0]


def group_project(keylist):
    for key in keylist:
        project_list = df[df.iloc[:,0].str.lower().str.contains(key)]
        project_size = project_list.shape[0]
        patients_num = project_list.sum(axis=0)[1]
        yield key
        yield project_list[["Name", "Reference", "CNA"]]
        yield project_size
        yield patients_num
        yield


keylist = ["breast", "prostate", "bladder", "lung", "GBM", "glio", "colo", "ovarian", "kidney", "uterian", "head",
            "thyroid", "stomach", "cervical", "sarcoma", "hepat", "leukemia", "esophageal", "pancrea", "rectum",
            "melanoma", "renal", "lymphoma"]
print len(keylist)

with open("project_info_group.txt", "w") as f:
    for i in group_project(keylist):
        print >>f, i

