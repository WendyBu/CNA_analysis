import pandas as pd
import glob
import os
import os.path


def pick20(sort_key, n):
    df_amp20 = pd.DataFrame()
    for file in glob.glob("reports/*"):
        df = pd.read_csv(file, sep="\t", index_col=0)
        tumor_type = os.path.basename(file).split("_")[0]
        df1 = df.copy()
        df1.sort_values(sort_key, ascending=False, inplace=True)
        df1_top20 = pd.Series(df1.ix[0:n,sort_key], name=tumor_type)
        if df_amp20.empty:
            df_amp20 = pd.concat([df_amp20, df1_top20], axis=1, sort=False)
        else:
            df_amp20 = df_amp20.join(df1_top20.to_frame(), how="outer")
    df_amp20.fillna('0', inplace=True)
    df_amp20.drop_duplicates(inplace=True)
    return df_amp20



df_top_amp = pick20(sort_key='deep_amp', n=20)
df_top_amp.to_csv("analysis/top_amp_chart.xls", sep='\t')

df_top_del = pick20(sort_key='deep_del', n=20)
df_top_del.to_csv("analysis/top_del_chart.xls", sep='\t')












