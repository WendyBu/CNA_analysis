"""
calculate the coefficiency among arms across all the tumor types.
"""



import pandas as pd
import numpy as np
import scipy.stats


def cal_corr(df):
    ### correlation : both methods works.
    # print scipy.stats.pearsonr(a, b)
    # print np.corrcoef(a, b)[0,1]
    s = []
    for i in df.index:
        for j in df.index:
            results = scipy.stats.pearsonr(df.ix[i, :], df.ix[j,:])
            row_data = [i, j]
            row_data.extend(list(results))
            s.append(row_data)
    df_s = pd.DataFrame(s,columns=["arm_1", "arm_2", "correlation", "p_value"])
    return df_s


def np_corr():
    pass






def main():
### generate both correlation and P value
    amp = "forms/amp_frequency.xls"
    amp_df = pd.read_csv(amp, sep="\t", index_col=0)
    # amp_df_corr = cal_corr(amp_df)
    # amp_df_corr.to_csv("secondForms/amp_correlations.xls", sep="\t", index=False)

    del_file = "forms/del_frequency.xls"
    del_df = pd.read_csv(del_file, sep="\t", index_col=0)
    # del_df_corr = cal_corr(del_df)
    # del_df_corr.to_csv("secondForms/del_correlations.xls", sep="\t", index=False)

### the third way: dataframe.corr(), only generate correlation, without p value
### triangle shape
    amp_correlation = amp_df.T.corr(method='pearson')
    amp_correlation.loc[:, :] = np.tril(amp_correlation, k=-1)
    amp_correlation.to_csv("secondForms/amp_corr2.xls", sep="\t")
    del_correlation = del_df.T.corr(method='pearson')
    del_correlation.loc[:, :] = np.tril(del_correlation, k=-1)
    del_correlation.to_csv("secondForms/del_corr2.xls", sep="\t")




if __name__ == "__main__":
    main()