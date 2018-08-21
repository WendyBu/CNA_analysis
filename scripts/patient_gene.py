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




def generate_heatmap_table(keyword="Melanoma", chr="chr8"):
    main_df = pd.DataFrame()
    input_path = os.path.join("..", keyword+"_temp", chr, "part_CNA")
    for file in glob.glob(input_path+"/*"):
        df_cna = pd.read_csv(file, sep="\t",index_col=0)
        if not main_df.empty:
            main_df = main_df.join(df_cna, lsuffix='', rsuffix='_y')
        else:
            main_df = df_cna
        # print main_df.shape
    main_df =  main_df.replace(r'\s+', np.nan, regex=True)
    main_df.fillna(0, inplace=True)
    output_path = os.path.join("..", keyword+"_temp", chr, "heatmap_"+keyword+"_"+chr+".xls")
    main_df.to_csv(output_path, sep="\t")



def main(keyword="Melanoma"):
    chr_list = list(range(1,23))
    chr_list.extend(["X", "Y"])
    chr_list = ["chr"+str(chr) for chr in chr_list]
    for chr in chr_list:
        generate_heatmap_table(keyword, chr)
    pass


if __name__ == "__main__":
    main()

