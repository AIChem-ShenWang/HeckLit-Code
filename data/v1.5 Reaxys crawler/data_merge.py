import pandas as pd
from utils.rxn_class import *
import os

# 1.打开文件夹
file_list = os.listdir("./data")
key = input("请输入反应名称")
data = list()
for file in file_list:
    if key in file and "raw data" in file:
        data.append(pd.read_excel("./data/" + file))

# 2.将文件转为rxn_class类型
rxn_list = list()
for df in data:
    for i in df_to_rxn_list(df):
        rxn_list.append(i)

# 3.清理reference
for rxn in rxn_list:
    rxn.ref = str(rxn.ref)

    if "Full Text" in rxn.ref:
        rxn.ref = rxn.ref.replace("Full Text", "")
    if "Details" in rxn.ref:
        rxn.ref = rxn.ref.replace("Details", "")
    if "Abstract" in rxn.ref:
        rxn.ref = rxn.ref.replace("Abstract", "")
    if "Cited" in rxn.ref:
        pos = rxn.ref.find("Cited")
        rxn.ref = rxn.ref[:(pos-1)]

# 4.将rxn_list转为df
df = rxn_list_to_df(rxn_list)
df.to_excel("./data/%s merged data.xlsx" % key)