from utils.molecule import *
import pandas as pd
from selenium import webdriver


# 0. set start row
st = int(input("Start cleaning from which row? (start with 0)"))
cnt = st
# 1.import data
if st == 0:
    data = pd.read_excel("../../data/Heck/Heck preprocessed data1.xlsx").copy()
else:
    data = pd.read_excel("../../data/Heck/Heck preprocessed data2.xlsx").copy()
len_rows = data.shape[0]
len_cols = data.shape[1] - 1

# 2.handle smi
# PubChem mode
path = "../../utils/geckodriver.exe"
bro = webdriver.Firefox(executable_path=path)
mol_manager = Mol_Manager(bro=bro)

# search all molecule
try:
    for i in range(st, len_rows):
        for col in data.columns:
            if "reagents" in col or "catalysts" in col or "solvents" in col:
                name = data[col][i]
                # if there is smi
                if mol_manager.get_smi(name) != None:
                    data.loc[i, col] = mol_manager.get_smi(name)
        cnt = i

except:
    print("Cleaned up to row %d" % (cnt))

bro.quit()

# 3.generate excel
data.to_excel("../../data/Heck/Heck preprocessed data2.xlsx")