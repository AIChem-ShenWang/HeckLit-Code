import pandas as pd
import re
from utils.rxn import *
from utils.molecule import *
import rdkit.Chem
import rdkit.Chem.AllChem

# import data
data = pd.read_excel("../data/Heck/Heck raw data.xlsx")
rxn_list = df_to_rxn_list(data)
rxn_selected = list()

cnt = 0
# select rxn
for rxn in rxn_list:
    # processing counter
    cnt += 1

    # at most two reactants
    if len(rxn.reactants) > 2:
       continue

    # only one product
    if len(rxn.products) > 1:
        continue

    # rxn without cat & reagent / sol
    if len(rxn.reagents) + len(rxn.cats) == 0 or len(rxn.solvents) == 0:
        continue

    # Normalization
    # reactants
    for i in range(len(rxn.reactants)):
        if Mol_Clean(rxn.reactants[i]) == False:
            continue
        rxn.reactants[i] = Mol_Clean(rxn.reactants[i])

    # products
    if Mol_Clean(rxn.products[0]) == False:
        continue
    else:
        rxn.products[0] = Mol_Clean(rxn.products[0])

    # yield
    if isinstance(rxn.rxn_yield, float) \
            or "%" not in rxn.rxn_yield \
            or rxn.rxn_yield == "/" \
            or rxn.rxn_yield == "" \
            or "A" in rxn.rxn_yield \
            or "g" in rxn.rxn_yield:
        continue

    rxn.rxn_yield = re.findall("-?[0-9]+\.[0-9]*|-?[0-9]+", rxn.rxn_yield)[0]

    # check it is Heck rxn

    if Heck_filter1(rxn) == False:
        continue

    # it is a Heck reaction
    rxn_selected.append(rxn)

    # process
    if cnt % 1000 == 0:
        print("Processing: %.2f" % (cnt / len(rxn_list) * 100), "%")

df = rxn_list_to_df(rxn_selected)
df.to_excel("../data/Heck/Heck preprocessed data1.xlsx")