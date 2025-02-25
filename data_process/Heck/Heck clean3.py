import random
from utils.rxn import *
from utils.molecule import *
import re

# 1. import data
data = pd.read_excel("../../data/Heck/Heck preprocessed data2.xlsx")
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

# # 5.Filling Missing value by interval temperature / time median
# # intra-molecular
# intra_yield_region_temp = list()
# intra_yield_region_time = list()
# # 10 yield region 0-10, 10-20 ... 90-100
# for i in range(10):
#     intra_yield_region_temp.append([])
#     intra_yield_region_time.append([])
# for rxn in rxn_intra:
#     for i in range(0, 100, 10):
#         if rxn.rxn_yield >= i and rxn.rxn_yield <= i + 10: # in the yield region
#             if rxn.temp != "/":
#                 intra_yield_region_temp[int(i / 10)].append(rxn.temp)
#             if rxn.time != "/":
#                 intra_yield_region_time[int(i / 10)].append(rxn.time)
# # Filling Missing Value by median value
# for rxn in rxn_intra:
#     if rxn.temp == "/":
#         for i in range(0, 100, 10):
#             if rxn.rxn_yield >= i and rxn.rxn_yield <= i + 10:  # in the yield region
#                 rxn.temp = np.median(np.array(intra_yield_region_temp[int(i / 10)]).astype(float), axis=0)
#                 break
#     if rxn.time == "/":
#         for i in range(0, 100, 10):
#             if rxn.rxn_yield >= i and rxn.rxn_yield <= i + 10:  # in the yield region
#                 rxn.time = np.median(np.array(intra_yield_region_time[int(i / 10)]).astype(float), axis=0)
#                 break
#
# # inter-molecular
# inter_yield_region_temp = list()
# inter_yield_region_time = list()
#
# # 10 yield region 0-10, 10-20 ... 90-100
# for i in range(10):
#     inter_yield_region_temp.append([])
#     inter_yield_region_time.append([])
# for rxn in rxn_inter:
#     for i in range(0, 100, 10):
#         if rxn.rxn_yield >= i and rxn.rxn_yield <= i + 10: # in the yield region
#             if rxn.temp != "/":
#                 inter_yield_region_temp[int(i / 10)].append(rxn.temp)
#             if rxn.time != "/":
#                 inter_yield_region_time[int(i / 10)].append(rxn.time)
#
# # Filling Missing Value by median value
# for rxn in rxn_inter:
#     if rxn.temp == "/":
#         for i in range(0, 100, 10):
#             if rxn.rxn_yield >= i and rxn.rxn_yield <= i + 10:  # in the yield region
#                 rxn.temp = np.median(np.array(inter_yield_region_temp[int(i / 10)]).astype(float), axis=0)
#                 break
#     if rxn.time == "/":
#         for i in range(0, 100, 10):
#             if rxn.rxn_yield >= i and rxn.rxn_yield <= i + 10:  # in the yield region
#                 rxn.time = np.median(np.array(inter_yield_region_time[int(i / 10)]).astype(float), axis=0)
#                 break
#
# rxn_list = rxn_intra + rxn_inter
# random.shuffle(rxn_list)

# 6.check smiles again and then create a new database
checked_rxn_list = list()
for rxn in rxn_list:
    if smi_checker(rxn):
        checked_rxn_list.append(rxn)

print("There are %d case(s) in the dataset" % len(checked_rxn_list))
df = rxn_list_to_df(checked_rxn_list)
df.to_excel("../../data/Heck/Heck processed data.xlsx")