import pandas
import pandas as pd
import matplotlib.pyplot as plt
from utils.rxn import *

# Color
data_color={"Intra":[114/255, 188/255, 213/255],
            "Inter":[55/255, 103/255, 149/255],
            "Total":'dodgerblue',
            "BH":[255/255, 208/255, 111/255],
            "Suzuki":[231/255, 98/255, 84/255]}

# 1. import data
# Intra
Heck_intra = pd.read_excel("../../data/Heck/Intramolecular data DR.xlsx")
rxn_intra = df_to_rxn_list(Heck_intra)

# MDS
# rxnfp
intra_rxnfp_MDS1 = Heck_intra["Intra_rxnfp_MDS1"]
intra_rxnfp_MDS2 = Heck_intra["Intra_rxnfp_MDS2"]

# drfp
intra_drfp_MDS1 = Heck_intra["Intra_drfp_MDS1"]
intra_drfp_MDS2 = Heck_intra["Intra_drfp_MDS2"]

# t-SNE
# reaction
# rxnfp
intra_rxnfp_TSNE1 = Heck_intra["Intra_rxnfp_TSNE1"]
intra_rxnfp_TSNE2 = Heck_intra["Intra_rxnfp_TSNE2"]
# drfp
intra_drfp_TSNE1 = Heck_intra["Intra_drfp_TSNE1"]
intra_drfp_TSNE2 = Heck_intra["Intra_drfp_TSNE2"]

# reactant and product
# rxnfp
intra_react_rxnfp_TSNE1 = Heck_intra["Intra_react_rxnfp_TSNE1"]
intra_react_rxnfp_TSNE2 = Heck_intra["Intra_react_rxnfp_TSNE2"]
# drfp
intra_react_drfp_TSNE1 = Heck_intra["Intra_react_drfp_TSNE1"]
intra_react_drfp_TSNE2 = Heck_intra["Intra_react_drfp_TSNE2"]

# reagents
# rxnfp
intra_reagent_rxnfp_TSNE1 = Heck_intra["Intra_reagent_rxnfp_TSNE1"]
intra_reagent_rxnfp_TSNE2 = Heck_intra["Intra_reagent_rxnfp_TSNE2"]
# drfp
intra_reagent_drfp_TSNE1 = Heck_intra["Intra_reagent_drfp_TSNE1"]
intra_reagent_drfp_TSNE2 = Heck_intra["Intra_reagent_drfp_TSNE2"]


# Inter
Heck_inter = pd.read_excel("../../data/Heck/Intermolecular data DR.xlsx")
rxn_inter = df_to_rxn_list(Heck_inter)

# MDS
# rxnfp
inter_rxnfp_MDS1 = Heck_inter["Inter_rxnfp_MDS1"]
inter_rxnfp_MDS2 = Heck_inter["Inter_rxnfp_MDS2"]
# drfp
inter_drfp_MDS1 = Heck_inter["Inter_drfp_MDS1"]
inter_drfp_MDS2 = Heck_inter["Inter_drfp_MDS2"]

# t-SNE
# reaction
# rxnfp
inter_rxnfp_TSNE1 = Heck_inter["Inter_rxnfp_TSNE1"]
inter_rxnfp_TSNE2 = Heck_inter["Inter_rxnfp_TSNE2"]
# drfp
inter_drfp_TSNE1 = Heck_inter["Inter_drfp_TSNE1"]
inter_drfp_TSNE2 = Heck_inter["Inter_drfp_TSNE2"]

# reactant and product
# rxnfp
inter_react_rxnfp_TSNE1 = Heck_inter["Inter_react_rxnfp_TSNE1"]
inter_react_rxnfp_TSNE2 = Heck_inter["Inter_react_rxnfp_TSNE2"]
# drfp
inter_react_drfp_TSNE1 = Heck_inter["Inter_react_drfp_TSNE1"]
inter_react_drfp_TSNE2 = Heck_inter["Inter_react_drfp_TSNE2"]

# reagents
# rxnfp
inter_reagent_rxnfp_TSNE1 = Heck_inter["Inter_reagent_rxnfp_TSNE1"]
inter_reagent_rxnfp_TSNE2 = Heck_inter["Inter_reagent_rxnfp_TSNE2"]
# drfp
inter_reagent_drfp_TSNE1 = Heck_inter["Inter_reagent_drfp_TSNE1"]
inter_reagent_drfp_TSNE2 = Heck_inter["Inter_reagent_drfp_TSNE2"]


# BH
BH_data = pd.read_excel("../../data/BH_HTE/BH_HTE DR.xlsx")
# MDS
# rxnfp
BH_rxnfp_MDS1 = BH_data["BH_rxnfp_MDS1"]
BH_rxnfp_MDS2 = BH_data["BH_rxnfp_MDS2"]

# drfp
BH_drfp_MDS1 = BH_data["BH_drfp_MDS1"]
BH_drfp_MDS2 = BH_data["BH_drfp_MDS2"]


# Suzuki
Suzuki_data = pd.read_excel("../../data/Suzuki_HTE/Suzuki_HTE DR.xlsx")
# MDS
# rxnfp
Suzuki_rxnfp_MDS1 = Suzuki_data["Suzuki_rxnfp_MDS1"]
Suzuki_rxnfp_MDS2 = Suzuki_data["Suzuki_rxnfp_MDS2"]

# drfp
Suzuki_drfp_MDS1 = Suzuki_data["Suzuki_drfp_MDS1"]
Suzuki_drfp_MDS2 = Suzuki_data["Suzuki_drfp_MDS2"]


# 2.Figures
# MDS
# rxnfp
plt.figure(dpi=500, figsize=(5, 5))

# inter
plt.scatter(inter_rxnfp_MDS1,
            inter_rxnfp_MDS2,
            color=data_color["Inter"],
            marker=".", s=75, alpha=0.8)

# intra
plt.scatter(intra_rxnfp_MDS1,
            intra_rxnfp_MDS2,
            color=data_color["Intra"],
            marker=".", s=75, alpha=0.8)

# BH
plt.scatter(BH_rxnfp_MDS1,
            BH_rxnfp_MDS2,
            color=data_color["BH"],
            marker=".", s=75, alpha=0.8)

# Suzuki
plt.scatter(Suzuki_rxnfp_MDS1,
            Suzuki_rxnfp_MDS2,
            color=data_color["Suzuki"],
            marker=".", s=75, alpha=0.8)

plt.xlabel("MDS1", fontsize=11)
plt.ylabel("MDS2", fontsize=11)
plt.xticks([])
plt.yticks([])
plt.title("MDS analysis(RXNFP)", fontsize=14)
plt.legend(["Intermolecular Heck", "Intramolecular Heck", "Buchwald HTE", "Suzuki HTE"], loc="best", prop={'size': 11})
plt.tight_layout()
plt.savefig("../../figures/RXNFP MDS Analysis.png")
# plt.show()

# drfp
plt.figure(dpi=500, figsize=(5, 5))
# inter
plt.scatter(inter_drfp_MDS1,
            inter_drfp_MDS2,
            color=data_color["Inter"],
            marker=".", s=75, alpha=0.8)

# intra
plt.scatter(intra_drfp_MDS1,
            intra_drfp_MDS2,
            color=data_color["Intra"],
            marker=".", s=75, alpha=0.8)

# BH
plt.scatter(BH_drfp_MDS1,
            BH_drfp_MDS2,
            color=data_color["BH"],
            marker=".", s=75, alpha=0.8)

# Suzuki
plt.scatter(Suzuki_drfp_MDS1,
            Suzuki_drfp_MDS2,
            color=data_color["Suzuki"],
            marker=".", s=75, alpha=0.8)

plt.xlabel("MDS1", fontsize=11)
plt.ylabel("MDS2", fontsize=11)
plt.xticks([])
plt.yticks([])
plt.title("MDS analysis(DRFP)", fontsize=14)
plt.legend(["Intermolecular Heck", "Intramolecular Heck", "Buchwald HTE", "Suzuki HTE"], loc="best", prop={'size': 11})
plt.tight_layout()
plt.savefig("../../figures/DRFP MDS Analysis.png")
# plt.show()

# t-SNE
# RXNFP
fig, ax = plt.subplots(ncols=3, dpi=500, figsize=(12, 5))
# Reaction
ax[0].scatter(inter_rxnfp_TSNE1,
             inter_rxnfp_TSNE2,
             color=data_color["Inter"],
             marker=".", s=75, alpha=0.8)
ax[0].scatter(intra_rxnfp_TSNE1,
             intra_rxnfp_TSNE2,
             color=data_color["Intra"],
             marker=".", s=75, alpha=0.8)
ax[0].set_ylabel("t-SNE2", fontsize=11)
ax[0].set_title("Reaction", fontsize=13)

# Only reactants and products
ax[1].scatter(inter_react_rxnfp_TSNE1,
             inter_react_rxnfp_TSNE2,
             color=data_color["Inter"],
             marker=".", s=75, alpha=0.8)
ax[1].scatter(intra_react_rxnfp_TSNE1,
             intra_react_rxnfp_TSNE2,
             color=data_color["Intra"],
             marker=".", s=75, alpha=0.8)
ax[1].set_title("Only Reactants and Products", fontsize=13)

# Only Reagents
ax[2].scatter(inter_reagent_rxnfp_TSNE1,
             inter_reagent_rxnfp_TSNE2,
             color=data_color["Inter"],
             marker=".", s=75, alpha=0.8)
ax[2].scatter(intra_reagent_rxnfp_TSNE1,
             intra_reagent_rxnfp_TSNE2,
             color=data_color["Intra"],
             marker=".", s=75, alpha=0.8)
ax[2].set_title("Only Reagents", fontsize=13)
ax[2].legend(["Intermolecular", "Intramolecular"], loc="upper right", prop={'size': 11})

# Beautify
for i in range(3):
    ax[i].set_xlabel("t-SNE1", fontsize=11)
    ax[i].set_xticks([])
    ax[i].set_yticks([])

plt.suptitle("Reaction Space described by t-SNE(RXNFP)", fontsize=16)
plt.tight_layout()
plt.savefig("../../figures/RXNFP t-SNE Analysis.png")
# plt.show()

# DRFP
fig, ax = plt.subplots(ncols=3, dpi=500, figsize=(12, 5))
# Reaction
ax[0].scatter(inter_drfp_TSNE1,
             inter_drfp_TSNE2,
             color=data_color["Inter"],
             marker=".", s=75, alpha=0.8)
ax[0].scatter(intra_drfp_TSNE1,
             intra_drfp_TSNE2,
             color=data_color["Intra"],
             marker=".", s=75, alpha=0.8)
ax[0].set_ylabel("t-SNE2", fontsize=11)
ax[0].set_title("Reaction", fontsize=13)

# Only reactants and products
ax[1].scatter(inter_react_drfp_TSNE1,
             inter_react_drfp_TSNE2,
             color=data_color["Inter"],
             marker=".", s=75, alpha=0.8)
ax[1].scatter(intra_react_drfp_TSNE1,
             intra_react_drfp_TSNE2,
             color=data_color["Intra"],
             marker=".", s=75, alpha=0.8)
ax[1].set_title("Only Reactants and Products", fontsize=13)

# Only Reagents
ax[2].scatter(inter_reagent_drfp_TSNE1,
             inter_reagent_drfp_TSNE2,
             color=data_color["Inter"],
             marker=".", s=75, alpha=0.8)
ax[2].scatter(intra_reagent_drfp_TSNE1,
             intra_reagent_drfp_TSNE2,
             color=data_color["Intra"],
             marker=".", s=75, alpha=0.8)
ax[2].set_title("Only Reagents", fontsize=13)
ax[2].legend(["Intermolecular", "Intramolecular"], loc="upper right", prop={'size': 11})

# Beautify
for i in range(3):
    ax[i].set_xlabel("t-SNE1", fontsize=11)
    ax[i].set_xticks([])
    ax[i].set_yticks([])

plt.suptitle("Reaction Space described by t-SNE(DRFP)", fontsize=16)
plt.tight_layout()
plt.savefig("../../figures/DRFP t-SNE Analysis.png")
# plt.show()