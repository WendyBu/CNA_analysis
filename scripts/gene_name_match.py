import pandas as pd
import glob
import os
import os.path

pd.set_option('display.max_columns', 100)


full_table = pd.read_csv("../genename_conversion_table.xls", sep="\t", index_col='Ensembl_ID1')
# print full_table.head(), full_table.columns
input_genelist = pd.read_csv("../chrom_gene_list/chr1_gene.csv", sep="\t", index_col=0)
# print input_genelist.shape  #5279 genes, include lincRNAs
df_gene_merge = input_genelist.merge(full_table, how="left", left_index=True, right_index=True)
# print df_gene_merge.head(), df_gene_merge.shape
df_gene_merge.drop_duplicates(inplace=True)
df_gene_ids = df_gene_merge[["Symbol", "Entrez_ID_x"]]
df_gene_ids.drop_duplicates(inplace=True)  # drop all lincRNAs, because no symbol or entrez_ID name. will not find match from cBioportal
# print df_gene_ids.shape  #3484
# print df_gene_ids.columns  #[u'Symbol', u'Entrez_ID_x']
input_list = df_gene_ids["Symbol"].tolist()
# print len(set(input_list))

#
# for file in glob.glob("/Users/yiwenbu/PycharmProjects/chrom/cBioportal_data/Breast/*/data_CNA.txt"):
df_cna = pd.read_csv("/Users/yiwenbu/PycharmProjects/chrom/cBioportal_data/Breast/brca_igr_2015/data_CNA.txt", sep="\t", index_col=0)
df_cna_input = df_cna.loc[input_list, :]
not_matched = set(input_list).difference(df_cna.index.values)

