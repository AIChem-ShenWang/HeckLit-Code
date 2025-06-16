import random
import numpy as np
from utils.rxn import *
from utils.molecule import *
import re

# 1. import data
data = pd.read_excel("../data/Heck/Heck preprocessed data2.xlsx")
raw_rxn_list = df_to_rxn_list(data)
smi_rxn_list = list()

# 2.check whether SMILES can be converted into mol, and then clean reagents & catalysts & solvents
for rxn in raw_rxn_list:
    # cleaning
    try:
        for i in range(len(rxn.reagents)):
            rxn.reagents[i] = Mol_Clean(rxn.reagents[i], ChooseLargeFrag=False, Uncharge=False)

        for i in range(len(rxn.cats)):
            rxn.cats[i] = Mol_Clean(rxn.cats[i], ChooseLargeFrag=False, Uncharge=False)

        for i in range(len(rxn.solvents)):
            rxn.solvents[i] = Mol_Clean(rxn.solvents[i], ChooseLargeFrag=False, Uncharge=False)

        if rxn.time != "/":
            if float(rxn.time) > 100:
                continue
    except:
        continue

    smi_rxn_list.append(rxn)

print("There are %d set(s) of data which contains smi can not convert to mol will be deleted" %(len(raw_rxn_list) - len(smi_rxn_list)))

# 3.check whether there is metal Pd in the rxn
rxn_list = list()
for rxn in smi_rxn_list:
    if Heck_filter2(rxn):
        rxn_list.append(rxn)

print("There are %d set(s) of data which do not contain Metal Pd will be deleted" %(len(smi_rxn_list) - len(rxn_list)))

# 4.divide into intra & inter
rxn_intra = list() # intra
rxn_inter = list() # inter
for rxn in rxn_list:
    if len(rxn.reactants) == 1:
        rxn_intra.append(rxn)
    if len(rxn.reactants) == 2:
        rxn_inter.append(rxn)
mol_intra = [Chem.MolFromSmiles(rxn.products[0]) for rxn in rxn_intra]
mol_inter = [Chem.MolFromSmiles(rxn.products[0]) for rxn in rxn_inter]

# 5.check smiles again and then create a new database
checked_rxn_list = list()
for rxn in rxn_list:
    if smi_checker(rxn):
        checked_rxn_list.append(rxn)

rxn_list = checked_rxn_list
print("There are %d case(s) in the dataset" % len(checked_rxn_list))

# 5.Calculate the average yield for the same reaction
rxn_dict = {}
for i in tqdm(range(len(rxn_list))):
    rxn = rxn_list[i]
    # Chemical dict
    chem_list = rxn.reactants + rxn.reagents + rxn.cats + rxn.solvents + rxn.products
    chem_dict = {}
    for chem in chem_list:
        if chem not in chem_dict.keys():
            chem_dict[chem] = 0
    chem_list = []
    for chem in chem_dict.keys():
        chem_list.append(chem)
    chem_names = str(chem_list)

    # if 2 reaction have the same chemicals
    for j in rxn_dict.keys():
        if j == chem_names:
            if isinstance(rxn_dict[j].rxn_id, list):
                if rxn.rxn_id not in rxn_dict[j].rxn_id:
                    # Yield
                    rxn_dict[j].rxn_yield.append(rxn.rxn_yield)
                    rxn.rxn_yield = rxn_dict[j].rxn_yield
                    # time
                    rxn_dict[j].time.append(rxn.time)
                    rxn.time = rxn_dict[j].time
                    # temp
                    rxn_dict[j].temp.append(rxn.temp)
                    rxn.temp = rxn_dict[j].temp
                    # rxn_id
                    rxn_dict[j].rxn_id.append(rxn.rxn_id)
                    rxn.rxn_id = rxn_dict[j].rxn_id
                    # ref
                    rxn_dict[j].ref.append(rxn.ref)
                    rxn.ref = rxn_dict[j].ref
            else:
                if rxn.rxn_id != rxn_dict[j].rxn_id:
                    # Yield
                    rxn.rxn_yield = [rxn.rxn_yield, rxn_dict[j].rxn_yield]
                    # time
                    rxn.time = [rxn.time, rxn_dict[j].time]
                    # temp
                    rxn.temp = [rxn.temp, rxn_dict[j].temp]
                    # rxn_id
                    rxn.rxn_id = [rxn.rxn_id, rxn_dict[j].rxn_id]
                    # ref
                    rxn.ref = [rxn.ref, rxn_dict[j].ref]

    rxn_dict[chem_names] = rxn

rxn_list = []

for key in rxn_dict.keys():
    rxn = rxn_dict[key]
    if isinstance(rxn.rxn_yield, list):
        rxn.rxn_yield = np.array([float(item) for item in rxn.rxn_yield]).mean()
        # temp
        try:
            if all(item == '/' for item in rxn.temp):
                rxn.temp = "/"
            else:
                rxn.temp = np.array([float(item) for item in rxn.temp if item != '/' and item is not None]).mean()
        except:
            rxn.temp = "/"
        # time
        try:
            if all(item == '/' for item in rxn.time):
                rxn.time = "/"
            else:
                rxn.time = np.array([float(item) for item in rxn.time if item != '/' and item is not None]).mean()
        except:
            rxn.time = "/"
    rxn_list.append(rxn_dict[key])

df = rxn_list_to_df(rxn_list)
df.to_excel("../data/Heck/Heck processed data.xlsx")