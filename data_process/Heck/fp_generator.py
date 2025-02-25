import pandas as pd
from utils.rxn import *
from tqdm import tqdm

# 1.import data
data_Heck = pd.read_excel("../../data/Heck/Heck processed data.xlsx")
Heck_rxn = df_to_rxn_list(data_Heck)
data_BH = pd.read_excel("../../data/BH_HTE/BH_HTE_data.xlsx")
data_Suzuki = pd.read_excel("../../data/Suzuki_HTE/Suzuki_HTE_data.xlsx")
np.set_printoptions(threshold=np.inf)

# 2.generate rxnfp
# Heck
Heck_rxnfp = list()
Heck_React_rxnfp = list()
Heck_Reagent_rxnfp = list()
Heck_drfp = list()
Heck_React_drfp = list()
Heck_Reagent_drfp = list()

for i in tqdm(range(len(Heck_rxn))):
    rxn = Heck_rxn[i]

    Heck_rxnfp.append(str(get_Heck_rxnfp(rxn)))
    Heck_React_rxnfp.append(str(get_Heck_React_rxnfp(rxn)))
    Heck_Reagent_rxnfp.append(str(get_Heck_Reagent_rxnfp(rxn)))

    Heck_drfp.append(str(get_Heck_drfp(rxn)))
    Heck_React_drfp.append(str(get_Heck_React_drfp(rxn)))
    Heck_Reagent_drfp.append(str(get_Heck_Reagent_drfp(rxn)))

Heck_rxnfp = pd.DataFrame(Heck_rxnfp, columns=["rxnfp"])
Heck_React_rxnfp = pd.DataFrame(Heck_React_rxnfp, columns=["rxnfp React"])
Heck_Reagent_rxnfp = pd.DataFrame(Heck_Reagent_rxnfp, columns=["rxnfp Reagent"])

Heck_drfp = pd.DataFrame(Heck_drfp, columns=["drfp"])
Heck_React_drfp = pd.DataFrame(Heck_React_drfp, columns=["drfp React"])
Heck_Reagent_drfp = pd.DataFrame(Heck_Reagent_drfp, columns=["drfp Reagent"])

Hcek_df = pd.concat([data_Heck,
                     Heck_rxnfp, Heck_React_rxnfp, Heck_Reagent_rxnfp,
                     Heck_drfp, Heck_React_drfp, Heck_Reagent_drfp], axis=1)
Hcek_df.to_excel("../../data/Heck/Heck_fp.xlsx")

# BH
BH_rxnfp = list()
BH_drfp = list()
for i in tqdm(range(data_BH.shape[0])):
    rxn = data_BH.loc[i]
    BH_rxnfp.append(str(get_Buchwald_rxnfp(rxn)))
    BH_drfp.append(str(get_Buchwald_drfp(rxn)))
BH_rxnfp = pd.DataFrame(BH_rxnfp, columns=["rxnfp"])
BH_drfp = pd.DataFrame(BH_drfp, columns=["drfp"])
BH_df = pd.concat([data_BH, BH_rxnfp, BH_drfp], axis=1)
BH_df.to_excel("../../data/BH_HTE/BH_HTE_fp.xlsx")

# Suzuki
Suzuki_rxnfp = list()
Suzuki_drfp = list()
for i in tqdm(range(data_Suzuki.shape[0])):
    rxn = data_Suzuki.loc[i]
    Suzuki_rxnfp.append(str(get_Suzuki_rxnfp(rxn)))
    Suzuki_drfp.append(str(get_Suzuki_drfp(rxn)))

Suzuki_rxnfp = pd.DataFrame(Suzuki_rxnfp, columns=["rxnfp"])
Suzuki_drfp = pd.DataFrame(Suzuki_drfp, columns=["drfp"])
Suzuki_df = pd.concat([data_Suzuki, Suzuki_rxnfp, Suzuki_drfp], axis=1)
Suzuki_df.to_excel("../../data/Suzuki_HTE/Suzuki_HTE_fp.xlsx")