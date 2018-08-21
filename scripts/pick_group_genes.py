## giving the tumor_type and a gene list, get a excel table on their amp, deep_amp, del, deep_del
## total amp, total del

import pandas as pd
import glob
import os
import os.path


def get_simple_name(filename):
    """
    split chromosome filename. input "chr1_gene.csv", output "chr1"
    """
    n1 = os.path.basename(filename)
    n2 = n1.split("_")[0]
    return n2


def trim_cna_df(tumor_type, genelist, file_save_name):
    """
    generate the rbp related CNA table, 1078 rows,
    Each gene one row.
    :return:  "Breast_CNA" folder
    """
    input_dir = os.path.join("..", "cBioportal_data", tumor_type)
    subdir_list = os.listdir(input_dir)
    if '.DS_Store' in subdir_list:
        subdir_list.remove('.DS_Store')
    for subfolder in subdir_list:
        f = os.path.join(input_dir, subfolder, "data_CNA.txt")
        # print f
        df = pd.read_table(f, sep="\t", index_col=0, low_memory=False)
        df_input = df[df.index.isin(genelist)]
        df_list = pd.Series(genelist)
        nomatch = df_list[~df_list.isin(df.index)]
        if df_input.columns[0] == "Entrez_Gene_Id":
            df_input.drop(columns="Entrez_Gene_Id", axis=0, inplace=True)
        sample_num = df_input.shape[1]
        new_folder = os.path.join("..", tumor_type + "_temp", file_save_name, "part_CNA/")
        if not os.path.exists(new_folder):
            os.makedirs(new_folder)
        outf = os.path.join(new_folder, subfolder + "_CNA.txt")
        # outf_nomatch = os.path.join(new_folder, "nomatch", subfolder + "no_match.txt")
        df_input.to_csv(outf, sep="\t")
        # nomatch.to_csv(outf_nomatch, sep="\t")
        yield sample_num


def generate_report(df_cna):
    """
    create report table for each study
    :param df_cna:
    :return:
    """
    df_cnaC = pd.DataFrame(index=df_cna.index)
    df_cnaC['Hetero_amp'] = df_cna[df_cna.iloc[:,:] == 1].count(axis=1)
    df_cnaC['deep_amp'] = df_cna[df_cna.iloc[:,:] == 2].count(axis=1)
    df_cnaC['hetero_del'] = df_cna[df_cna.iloc[:,:] == -1].count(axis=1)
    df_cnaC['deep_del'] = df_cna[df_cna.iloc[:,:] == -2].count(axis=1)
    df_cnaC['amp_total'] = df_cnaC['sum_1'] + df_cnaC['deep_amp']
    df_cnaC['del_total'] = df_cnaC['sum_-1'] + df_cnaC['deep_del']
    return df_cnaC


def save_report(keyword, simple_name):
    """
    Save CNA information report in
    :return: "Breast_report" folder for each study
    """
    input_folder = os.path.join("..", keyword + "_temp", simple_name, "part_CNA")
    print input_folder
    for file in glob.glob(input_folder+"/*"):
        df_cna = pd.read_table(file, sep="\t", index_col=0)
        df_cna_report = generate_report(df_cna)

        out_folder = os.path.join("..", keyword + "_temp", simple_name, "sum_CNA")
        if not os.path.exists(out_folder):
            os.makedirs(out_folder)
        filename = os.path.split(file)[1]
        output_name = os.path.join(out_folder, filename)
        df_cna_report.to_csv(output_name, sep="\t")
        yield df_cna_report


def generate_genelist(chr_range):  # give a range in chrom, find the genelist
    chr8 = pd.read_table("GRch38_chr8_genes.txt", sep="\t", index_col=0)
    KIAA_pos = chr8.loc["KIAA1429", "chromStart"]
    up = KIAA_pos - chr_range
    down = KIAA_pos + chr_range
    gene_range = chr8[chr8["chromStart"].between(up, down, inclusive=True)]
    genelist = gene_range.index.tolist()
    return genelist


def create_genelist(chrom_file):
    input_file = pd.read_csv(chrom_file, sep="\t", index_col=0)
    return input_file["gene_name"].tolist()


def main(keyword="Breast"):
    # step1: generate rbp_related CNA table
    # genelist1 = generate_genelist(chr_range)   # give a range, find the genelist
    filename = "../chrom_gene_list/chr8_gene.csv"
    simple_name = get_simple_name(filename)
    genelist1 = create_genelist(filename)

    total_sample_number = 0
    for n in trim_cna_df(keyword, genelist1, simple_name):
        total_sample_number += n  ## total sample number in all Breast patients

    # step2: Count all the breast samples CNA
    df_total = pd.DataFrame()
    for i in save_report(keyword, simple_name):
        if df_total.empty:
            df_total = i
        else:
            df_total = df_total.add(i, fill_value=0)
    df_total.drop_duplicates(inplace=True)
    #
    # step3: calculate frequency
    print "total sample number:", total_sample_number
    total_sample_number /= 100.00
    df_freq = df_total.divide(total_sample_number)
    #
    # ranking the genename by original sequence
    pos = pd.read_table(filename, sep="\t", index_col=0)
    df_ordered = df_freq.merge(pos, how="left", left_index = True, right_on = "gene_name")
    df_ordered["gene_middle"] = (df_ordered["start"] + df_ordered["end"])/2
    # virma_center = df_ordered.loc["KIAA1429","Middle"]
    # df_ordered["centered"] = df_ordered["Middle"] - virma_center

    # save the file
    filepath_dir = os.path.join("..", keyword + "_temp", "reports")
    filepath2 = os.path.join(filepath_dir, simple_name + ".xls" )
    if not os.path.exists(filepath_dir):
        os.makedirs(filepath_dir)
    df_ordered.to_csv(filepath2, sep="\t")
    pass


## keyword : "Breast"  6290
## keyword : "Prostate"  3836
# "Ovarian
# Melanoma  896
# Bladder 1696
# Colon 1310
# GBM 5022
# Head&Neck 1470
# Kidney  2303
# Liver 968
# Lung 3989
# Stomach 1900
# Thyroid 1610
# Uterine 1426




if __name__ == "__main__":
    main()


### myc location NC_000008.11 (127735434..127742951)
