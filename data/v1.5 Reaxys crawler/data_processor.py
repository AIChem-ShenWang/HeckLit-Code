from utils.molecule_manager import *
import pandas as pd
from selenium import webdriver
import re

# 1.导入数据
data = pd.read_excel("./data/raw data.xlsx").copy()
len_rows = data.shape[0]
len_cols = data.shape[1] - 1

# 2.处理分子smi式
# 询问是否使用cactus
pd = input("是否同时使用cactus, 若需要请输入 y")
if pd == "y":
    pd = True
else:
    pd = False
# PubChem模拟浏览器初始化
bro = webdriver.Chrome("chromedriver.exe")
mol_manager = Mol_Manager(bro=bro)
# 遍历所有非smi的分子式
for col in data.columns:
    if "reagents" in col or "catalysts" in col or "solvents" in col:
        for i in range(len_rows):
            name = data[col][i]
            # 若能找到smi
            if mol_manager.get_smi(name) != None:
                data.loc[i, col] = mol_manager.get_smi(name, cactus=pd)
bro.quit()

# 3.处理产率栏
for i in range(len_rows):
    y = data["Yield"][i]
    y = re.findall("-?[0-9]+\.[0-9]*|-?[0-9]+",y)
    if len(y) != 0:
        data.loc[i, "Yield"] = y[0]

# 4.生成处理后excel
data.to_excel("./data/processed data.xlsx", sheet_name="processed data")