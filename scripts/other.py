"""
combine amplification and deletion together
results stored in secondForms/
"""

import pandas as pd

## 1.combine amp_freq and del_freq
amp_freq = pd.read_csv("forms/amp_frequency.xls", sep="\t")
del_freq = pd.read_csv("forms/del_frequency.xls", sep="\t")

amp_freq["Arm"] = amp_freq["Arm"].apply(lambda x: x + str("_amp"))
amp_freq.set_index("Arm", inplace=True)

del_freq["Arm"] = del_freq["Arm"].apply(lambda x: x + str("_del"))
del_freq.set_index("Arm", inplace=True)
# del_freq = del_freq.applymap(lambda x: -x)
total_freq = amp_freq.append(del_freq)

total_freq.to_csv("secondForms/total_frequency2.xls", sep="\t")




