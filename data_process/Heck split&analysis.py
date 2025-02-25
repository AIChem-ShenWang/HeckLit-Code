import pandas as pd
from sklearn.manifold import MDS
from utils.rxn import *
from utils.dataset_analysis import *
import matplotlib.pyplot as plt
from tqdm import tqdm
from sklearn.manifold import TSNE
import os
import umap
np.set_printoptions(threshold=np.inf)

# 1.import data
data_Heck = pd.read_excel("../../data/Heck/Heck_fp.xlsx")
data_BH = pd.read_excel("../../data/BH_HTE/BH_HTE_fp.xlsx")
data_Suzuki = pd.read_excel("../../data/Suzuki_HTE/Suzuki_HTE_fp.xlsx")


# 2.convert df into rxn_class list
rxn_list = df_to_rxn_list(data_Heck)


# 3.divide into intra & inter and other subsets
rxn_intra = list() # intra
rxn_inter = list() # inter
# The RXNFP & DRFP are in str() type and need to be converted into numpy using read_rxnfp or read_drfp
# rxnfp
rxnfp_intra = list()
rxnfp_inter = list()
# drfp
drfp_intra = list()
drfp_inter = list()
# Split into Intramolecular & Intermolecular
for i in range(len(rxn_list)):
    rxn = rxn_list[i]
    if len(rxn.reactants) == 1:
        rxn_intra.append(rxn)
        rxnfp_intra.append(data_Heck.loc[i]["rxnfp"])
        drfp_intra.append(data_Heck.loc[i]["drfp"])
    if len(rxn.reactants) == 2:
        rxn_inter.append(rxn)
        rxnfp_inter.append(data_Heck.loc[i]["rxnfp"])
        drfp_inter.append(data_Heck.loc[i]["drfp"])

data_color={"Intra":[114/255, 188/255, 213/255],
            "Inter":[55/255, 103/255, 149/255],
            "Total":'dodgerblue',
            "BH":[255/255, 208/255, 111/255],
            "Suzuki":[231/255, 98/255, 84/255]}

df_intra = rxn_list_to_df(rxn_intra)
df_inter = rxn_list_to_df(rxn_inter)
# Split datasets by reaction diversity
if os.path.exists("../../data/Heck/Subsets") == False:
    os.mkdir("../../data/Heck/Subsets")

# Intra
Intra_Insertion = [[], []] # Alpha, Beta
Intra_LG = [[], [], [], [], []] # F, Cl, Br, I, *
Intra_RS = [[], [], [], [], []] # 5, 6, 7, 8, 9
for i in tqdm(range(len(rxn_intra))):
    rxn = rxn_intra[i]
    Type_Add, Type_LG, MemberSize = HeckIntra_Classify(rxn)
    # Insertion Type
    if Type_Add == "Alpha":
        Intra_Insertion[0].append([rxn, rxnfp_intra[i], drfp_intra[i]])
    if Type_Add == "Beta":
        Intra_Insertion[1].append([rxn, rxnfp_intra[i], drfp_intra[i]])
    # LG
    if Type_LG == "F":
        Intra_LG[0].append([rxn, rxnfp_intra[i], drfp_intra[i]])
    if Type_LG == "Cl":
        Intra_LG[1].append([rxn, rxnfp_intra[i], drfp_intra[i]])
    if Type_LG == "Br":
        Intra_LG[2].append([rxn, rxnfp_intra[i], drfp_intra[i]])
    if Type_LG == "I":
        Intra_LG[3].append([rxn, rxnfp_intra[i], drfp_intra[i]])
    if Type_LG == "*":
        Intra_LG[4].append([rxn, rxnfp_intra[i], drfp_intra[i]])
    # Ring Size
    RS = MemberSize[0]
    if RS == 5:
        Intra_RS[0].append([rxn, rxnfp_intra[i], drfp_intra[i]])
    if RS == 6:
        Intra_RS[1].append([rxn, rxnfp_intra[i], drfp_intra[i]])
    if RS == 7:
        Intra_RS[2].append([rxn, rxnfp_intra[i], drfp_intra[i]])
    if RS == 8:
        Intra_RS[3].append([rxn, rxnfp_intra[i], drfp_intra[i]])
    if RS == 9:
        Intra_RS[4].append([rxn, rxnfp_intra[i], drfp_intra[i]])
# Generate excel
# Insertion
InsertionType = ["Alpha", "Beta"]
for i in range(len(Intra_Insertion)):
    pd.concat(
        [rxn_list_to_df([case[0] for case in Intra_Insertion[i]]),
           pd.DataFrame([case[1] for case in Intra_Insertion[i]], columns=["rxnfp"]),
           pd.DataFrame([case[2] for case in Intra_Insertion[i]], columns=["drfp"])], axis=1).to_excel(
        "../../data/Heck/Subsets/Intramolecular %s-Insertion.xlsx" % InsertionType[i])
# LG
LGType = ["F", "Cl", "Br", "I", "other"]
for i in range(len(Intra_LG)):
    pd.concat(
        [rxn_list_to_df([case[0] for case in Intra_LG[i]]),
           pd.DataFrame([case[1] for case in Intra_LG[i]], columns=["rxnfp"]),
           pd.DataFrame([case[2] for case in Intra_LG[i]], columns=["drfp"])], axis=1).to_excel(
        "../../data/Heck/Subsets/Intramolecular LG = %s.xlsx" % LGType[i])
# Ring Size
RSType = [5, 6, 7, 8, 9]
for i in range(len(Intra_RS)):
    pd.concat(
        [rxn_list_to_df([case[0] for case in Intra_RS[i]]),
           pd.DataFrame([case[1] for case in Intra_RS[i]], columns=["rxnfp"]),
           pd.DataFrame([case[2] for case in Intra_RS[i]], columns=["drfp"])], axis=1).to_excel(
        "../../data/Heck/Subsets/Intramolecular %s-membered Ring.xlsx" % RSType[i])

# Inter
Inter_Insertion = [[], []] # Alpha, Beta
Inter_LG = [[], [], [], [], []] # F, Cl, Br, I, *
for i in tqdm(range(len(rxn_inter))):
    rxn = rxn_inter[i]
    Type_Add, Type_LG = HeckInter_Classify(rxn)
    # Insertion Type
    if Type_Add == "Alpha":
        Inter_Insertion[0].append([rxn, rxnfp_inter[i], drfp_inter[i]])
    if Type_Add == "Beta":
        Inter_Insertion[1].append([rxn, rxnfp_inter[i], drfp_inter[i]])
    # LG
    if Type_LG == "F":
        Inter_LG[0].append([rxn, rxnfp_inter[i], drfp_inter[i]])
    if Type_LG == "Cl":
        Inter_LG[1].append([rxn, rxnfp_inter[i], drfp_inter[i]])
    if Type_LG == "Br":
        Inter_LG[2].append([rxn, rxnfp_inter[i], drfp_inter[i]])
    if Type_LG == "I":
        Inter_LG[3].append([rxn, rxnfp_inter[i], drfp_inter[i]])
    if Type_LG == "*":
        Inter_LG[4].append([rxn, rxnfp_inter[i], drfp_inter[i]])
# Generate excel
# Insertion
InsertionType = ["Alpha", "Beta"]
for i in range(len(Inter_Insertion)):
    pd.concat(
        [rxn_list_to_df([case[0] for case in Inter_Insertion[i]]),
           pd.DataFrame([case[1] for case in Inter_Insertion[i]], columns=["rxnfp"]),
           pd.DataFrame([case[2] for case in Inter_Insertion[i]], columns=["drfp"])], axis=1).to_excel(
        "../../data/Heck/Subsets/Intermolecular %s-Insertion.xlsx" % InsertionType[i])
# LG
LGType = ["F", "Cl", "Br", "I", "other"]
for i in range(len(Inter_LG)):
    pd.concat(
        [rxn_list_to_df([case[0] for case in Inter_LG[i]]),
           pd.DataFrame([case[1] for case in Inter_LG[i]], columns=["rxnfp"]),
           pd.DataFrame([case[2] for case in Inter_LG[i]], columns=["drfp"])], axis=1).to_excel(
        "../../data/Heck/Subsets/Intermolecular LG = %s.xlsx" % LGType[i])

# 4.counting for rxn
report = open("Dataset Report.txt", mode="w+")

# intra
intra_res = IntraCount(rxn_intra, file=report)

# inter
inter_res = InterCount(rxn_inter, file=report)

# Total
total_res = HeckCount(rxn_list, file=report)

report.close()

# Counting Figure
fig, ax = plt.subplots(ncols=2, dpi=500, figsize=(13, 5))

# Intra
label = list()
num = list()
for key in intra_res.keys():
    label.append(key)
    num.append(intra_res[key] / len(rxn_intra) * 100)
# ax[0].bar(x=[i for i in range(len(label))], height=num, color=data_color["Intra"])
for i in range(len(num)):
    # Insertion Type
    if i < 2:
        ax[0].bar(x=i, height=num[i], color="#fce8dd", edgecolor="#e9212c")
    # LG Type
    if i >= 2 and i <= 5:
        ax[0].bar(x=i, height=num[i], color="#d2e7d4", edgecolor="#01844f")
    # RingSize
    if i > 5:
        ax[0].bar(x=i, height=num[i], color="#d0d8eb", edgecolor="#7195c5")
for i in range(len(label)):
    ax[0].text(i, num[i] + 0.6, '%.1f' % num[i] + "%", ha='center', va='center', fontsize=9)

ax[0].set_title("Intramolecular Dataset", fontsize=16)
ax[0].set_ylabel("% of Intramolecular Dataset", fontsize=16)
ax[0].set_xticks([i for i in range(len(label))], label, rotation=20, fontsize=9)

# Inter
label = list()
num = list()
for key in inter_res.keys():
    label.append(key)
    num.append(inter_res[key] / len(rxn_inter) * 100)
for i in range(len(num)):
    # Insertion Type
    if i < 2:
        ax[1].bar(x=i, height=num[i], color="#fce8dd", edgecolor="#e9212c")
    # LG Type
    if i >= 2:
        ax[1].bar(x=i, height=num[i], color="#d2e7d4", edgecolor="#01844f")

# label
for i in range(len(label)):
    ax[1].text(i, num[i] + 1.2, '%.1f' % num[i] + "%", ha='center', va='center', fontsize=9)

ax[1].set_title("Intermolecular Dataset", fontsize=16)
ax[1].set_ylabel("% of Intermolecular Dataset", fontsize=16)
ax[1].set_xticks([i for i in range(len(label))], label, rotation=20, fontsize=9)

plt.suptitle("Heck Reaction Diversity", fontsize=18)
plt.tight_layout()
plt.savefig("../../figures/Reaction Diversity.png")
# plt.show()


# 5.analysis of temp, time & yield distribution
# intra
intra_time_dis = list()
intra_temp_dis = list()
intra_yield_dis = list()
for rxn in rxn_intra:
    if rxn.time != "/":
        intra_time_dis.append(rxn.time)
    if rxn.temp != "/":
        intra_temp_dis.append(rxn.temp)
    intra_yield_dis.append(rxn.rxn_yield)

# inter
inter_time_dis = list()
inter_temp_dis = list()
inter_yield_dis = list()
for rxn in rxn_inter:
    if rxn.time != "/":
        inter_time_dis.append(rxn.time)
    if rxn.temp != "/":
        inter_temp_dis.append(rxn.temp)
    inter_yield_dis.append(rxn.rxn_yield)

total_yield_dis = list()
for rxn in rxn_list:
    total_yield_dis.append(rxn.rxn_yield)

attr = ["Time", "Temperature", "Yield"]
name = ["Intramolecular", "Intermolecular"]

fig, ax = plt.subplots(nrows=3, dpi=500, figsize=(8, 12))
ax[0].boxplot(np.array(intra_time_dis).astype(float), positions=[1])
ax[0].boxplot(np.array(inter_time_dis).astype(float), positions=[1.8])

ax[1].boxplot(np.array(intra_temp_dis).astype(float), positions=[1])
ax[1].boxplot(np.array(inter_temp_dis).astype(float), positions=[1.8])

vio_part = ['cbars', 'cmins', 'cmaxes', 'cmedians']
vio = ax[2].violinplot(intra_yield_dis, positions=[1], showmedians=True)
for vp in vio['bodies']:
    vp.set_facecolor(data_color["Intra"])
    vp.set_edgecolor('#000')
    vp.set_alpha(.9)
for p in vio_part:
    vp = vio[p]
    vp.set_edgecolor('#000')
    vp.set_alpha(.9)

vio = ax[2].violinplot(inter_yield_dis, positions=[1.8], showmedians=True)
for vp in vio['bodies']:
    vp.set_facecolor(data_color["Inter"])
    vp.set_edgecolor('#000')
    vp.set_alpha(.9)
for p in vio_part:
    vp = vio[p]
    vp.set_edgecolor('#000')
    vp.set_alpha(.9)

# beautify
for i in range(3):
    ax[i].set_ylabel("%s" % attr[i], fontsize=14)
    ax[i].set_xticks([1, 1.8], name, fontsize=14)
    ax[i].set_title("%s Distribution" % attr[i], fontsize=16)
fig.suptitle("Time/Temperature/Yield Distribution for Heck Reaction", fontsize=18)
plt.tight_layout()
plt.savefig("../../figures/Time Temperature Yield Distribution.png")
# plt.show()


# 6.Yield Distribution Compare with HTE
BH_yield_dis = list()
for i in range(data_BH.shape[0]):
    BH_yield_dis.append(float(data_BH.loc[i]["yield"]))

Suzuki_yield_dis = list()
for i in range(data_Suzuki.shape[0]):
    Suzuki_yield_dis.append(float(data_Suzuki.loc[i]["Product_Yield_PCT_Area_UV"]))

plt.figure(dpi=500, figsize=(9, 5))
vio = plt.violinplot(total_yield_dis, positions=[1.0], showmedians=True)
for vp in vio['bodies']:
    vp.set_facecolor(data_color["Total"])
    vp.set_edgecolor('#000')
    vp.set_alpha(.9)
for p in vio_part:
    vp = vio[p]
    vp.set_edgecolor('#000')
    vp.set_alpha(.9)

vio = plt.violinplot(BH_yield_dis, positions=[1.6], showmedians=True)
for vp in vio['bodies']:
    vp.set_facecolor(data_color["BH"])
    vp.set_edgecolor('#000')
    vp.set_alpha(.9)
for p in vio_part:
    vp = vio[p]
    vp.set_edgecolor('#000')
    vp.set_alpha(.9)

vio = plt.violinplot(Suzuki_yield_dis, positions=[2.2], showmedians=True)
for vp in vio['bodies']:
    vp.set_facecolor(data_color["Suzuki"])
    vp.set_edgecolor('#000')
    vp.set_alpha(.9)
for p in vio_part:
    vp = vio[p]
    vp.set_edgecolor('#000')
    vp.set_alpha(.9)

plt.xticks([1, 1.6, 2.2], ["Heck", "Buchwald HTE", "Suzuki HTE"], fontsize=10)
plt.ylabel("Yield(%)", fontsize=10)
plt.title("Yield for Different Datasets", fontsize=12)
plt.tight_layout()
plt.savefig("../../figures/Yield Distribution Comparison.png")
plt.show()


# 7. Yield distribution with reaction number
shot_color = {
    "Few-shot":[210/255, 34/255, 37/255],
    "Medium-shot":[16/255, 118/255, 63/255],
    "Many-shot":[37/255, 48/255, 122/255]
}
split_num = 100 # 0-100% is split into split_num regions
BH_yield_dict = shot_classifier(BH_yield_dis, split_num, upper=50, lower=25)
Suzuki_yield_dist = shot_classifier(Suzuki_yield_dis, split_num, upper=65, lower=20)
intra_yield_dist = shot_classifier(intra_yield_dis, split_num, upper=20, lower=5)
inter_yield_dist = shot_classifier(inter_yield_dis, split_num, upper=125, lower=30)
total_yield_dist = shot_classifier(total_yield_dis, split_num, upper=150, lower=50)
yield_dict_list = [BH_yield_dict,
                   Suzuki_yield_dist,
                   intra_yield_dist,
                   inter_yield_dist,
                   total_yield_dist]
yield_name = ["Buchwald HTE", "Suzuki HTE", "Heck Intramolecular", "Heck Intermolecular", "Heck"]

# Figure of Yield Distribution
# BH
for i in range(5):
    yield_dict = yield_dict_list[i]
    type = [[], [], []]
    plt.figure(dpi=500)
    for key in yield_dict.keys():
        color = shot_color[yield_dict[key][1]]
        b = plt.bar(x=key, height=yield_dict[key][0], width=100/split_num, color=color)

        if yield_dict[key][1] == "Few-shot":
            type[0] = b
        if yield_dict[key][1] == "Medium-shot":
            type[1] = b
        if yield_dict[key][1] == "Many-shot":
            type[2] = b

    plt.title("Yield for %s Dataset" % yield_name[i], fontsize=12)
    plt.xlabel("Reaction Number", fontsize=8)
    plt.ylabel("Yield", fontsize=8)
    plt.legend(type, ["Few-shot", "Medium-shot", "Many-shot"])
    plt.tight_layout()
    plt.savefig("../../figures/%s Yield Distribution.png" % yield_name[i])


# 8.MDS analysis for RXNFP & DRFP
# RXNFP
# Heck
HeckInter_rxnfp = list()
HeckIntra_rxnfp = list()
HeckInter_react_rxnfp = list()
HeckIntra_react_rxnfp = list()
HeckInter_reagent_rxnfp = list()
HeckIntra_reagent_rxnfp = list()
for i in range(data_Heck.shape[0]):
    if data_Heck.loc[i]["reactants 2"] == "/":
        HeckIntra_rxnfp.append(read_rxnfp(data_Heck.loc[i]["rxnfp"]))
        HeckIntra_react_rxnfp.append(read_rxnfp(data_Heck.loc[i]["rxnfp React"]))
        HeckIntra_reagent_rxnfp.append(read_rxnfp(data_Heck.loc[i]["rxnfp Reagent"]))
    else:
        HeckInter_rxnfp.append(read_rxnfp(data_Heck.loc[i]["rxnfp"]))
        HeckInter_react_rxnfp.append(read_rxnfp(data_Heck.loc[i]["rxnfp React"]))
        HeckInter_reagent_rxnfp.append(read_rxnfp(data_Heck.loc[i]["rxnfp Reagent"]))

# BH
BH_rxnfp = list()
for i in range(data_BH.shape[0]):
    BH_rxnfp.append(read_rxnfp(data_BH.loc[i]["rxnfp"]))


# Suzuki
Suzuki_rxnfp = list()
for i in range(data_Suzuki.shape[0]):
    Suzuki_rxnfp.append(read_rxnfp(data_Suzuki.loc[i]["rxnfp"]))

rxnfp_list = HeckInter_rxnfp + HeckIntra_rxnfp + BH_rxnfp + Suzuki_rxnfp

xy = fp_mds(rxnfp_list, 2)


# inter
accum = len(HeckInter_rxnfp)
inter_rxnfp_MDS1 = pd.DataFrame(xy[:accum, 0], columns=["Inter_rxnfp_MDS1"])
inter_rxnfp_MDS2 = pd.DataFrame(xy[:accum, 1], columns=["Inter_rxnfp_MDS2"])
df_inter = pd.concat([df_inter, inter_rxnfp_MDS1, inter_rxnfp_MDS2], axis=1)

# intra
accum += len(HeckIntra_rxnfp)
intra_rxnfp_MDS1 = pd.DataFrame(xy[accum - len(HeckIntra_rxnfp):accum, 0], columns=["Intra_rxnfp_MDS1"])
intra_rxnfp_MDS2 = pd.DataFrame(xy[accum - len(HeckIntra_rxnfp):accum, 1], columns=["Intra_rxnfp_MDS2"])
df_intra = pd.concat([df_intra, intra_rxnfp_MDS1, intra_rxnfp_MDS2], axis=1)

# BH
accum += len(BH_rxnfp)
BH_rxnfp_MDS1 = pd.DataFrame(xy[accum - len(BH_rxnfp):accum, 0], columns=["BH_rxnfp_MDS1"])
BH_rxnfp_MDS2 = pd.DataFrame(xy[accum - len(BH_rxnfp):accum, 1], columns=["BH_rxnfp_MDS2"])
data_BH = pd.concat([data_BH, BH_rxnfp_MDS1, BH_rxnfp_MDS2], axis=1)

# Suzuki
Suzuki_rxnfp_MDS1 = pd.DataFrame(xy[accum:, 0], columns=["Suzuki_rxnfp_MDS1"])
Suzuki_rxnfp_MDS2 = pd.DataFrame(xy[accum:, 1], columns=["Suzuki_rxnfp_MDS2"])
data_Suzuki = pd.concat([data_Suzuki, Suzuki_rxnfp_MDS1, Suzuki_rxnfp_MDS2], axis=1)

# DRFP
# Heck
HeckInter_drfp = list()
HeckIntra_drfp = list()
HeckInter_react_drfp = list()
HeckIntra_react_drfp = list()
HeckInter_reagent_drfp = list()
HeckIntra_reagent_drfp = list()
for i in range(data_Heck.shape[0]):
    if data_Heck.loc[i]["reactants 2"] == "/":
        HeckIntra_drfp.append(read_drfp(data_Heck.loc[i]["drfp"]))
        HeckIntra_react_drfp.append(read_drfp(data_Heck.loc[i]["drfp React"]))
        HeckIntra_reagent_drfp.append(read_drfp(data_Heck.loc[i]["drfp Reagent"]))
    else:
        HeckInter_drfp.append(read_drfp(data_Heck.loc[i]["drfp"]))
        HeckInter_react_drfp.append(read_drfp(data_Heck.loc[i]["drfp React"]))
        HeckInter_reagent_drfp.append(read_drfp(data_Heck.loc[i]["drfp Reagent"]))

# BH
BH_drfp = list()
for i in range(data_BH.shape[0]):
    BH_drfp.append(read_drfp(data_BH.loc[i]["drfp"]))


# Suzuki
Suzuki_drfp = list()
for i in range(data_Suzuki.shape[0]):
    Suzuki_drfp.append(read_drfp(data_Suzuki.loc[i]["drfp"]))

drfp_list = HeckInter_drfp + HeckIntra_drfp + BH_drfp + Suzuki_drfp

xy = fp_mds(drfp_list, 2)


# Heck
# inter
accum = len(HeckInter_drfp)
inter_drfp_MDS1 = pd.DataFrame(xy[:accum, 0], columns=["Inter_drfp_MDS1"])
inter_drfp_MDS2 = pd.DataFrame(xy[:accum, 1], columns=["Inter_drfp_MDS2"])
df_inter = pd.concat([df_inter, inter_drfp_MDS1, inter_drfp_MDS2], axis=1)

# intra
accum += len(HeckIntra_drfp)
intra_drfp_MDS1 = pd.DataFrame(xy[accum - len(HeckIntra_drfp):accum, 0], columns=["Intra_drfp_MDS1"])
intra_drfp_MDS2 = pd.DataFrame(xy[accum - len(HeckIntra_drfp):accum, 1], columns=["Intra_drfp_MDS2"])
df_intra = pd.concat([df_intra, intra_drfp_MDS1, intra_drfp_MDS2], axis=1)

# BH
accum += len(BH_drfp)
BH_drfp_MDS1 = pd.DataFrame(xy[accum - len(BH_drfp):accum, 0], columns=["BH_drfp_MDS1"])
BH_drfp_MDS2 = pd.DataFrame(xy[accum - len(BH_drfp):accum, 1], columns=["BH_drfp_MDS2"])
data_BH = pd.concat([data_BH, BH_drfp_MDS1, BH_drfp_MDS2], axis=1)

# Suzuki
Suzuki_drfp_MDS1 = pd.DataFrame(xy[accum:, 0], columns=["Suzuki_drfp_MDS1"])
Suzuki_drfp_MDS2 = pd.DataFrame(xy[accum:, 1], columns=["Suzuki_drfp_MDS2"])
data_Suzuki = pd.concat([data_Suzuki, Suzuki_drfp_MDS1, Suzuki_drfp_MDS2], axis=1)


# 9.Cosine Similarity density
import seaborn as sns
# RXNFP
plt.figure(dpi=500)

sns.kdeplot(fp_CosDensity(HeckInter_rxnfp), color=[55/255, 103/255, 149/255], label="HecK Intermolecular", fill=True)
sns.kdeplot(fp_CosDensity(HeckIntra_rxnfp), color=[114/255, 188/255, 213/255], label="Heck Intramolecular", fill=True)
sns.kdeplot(fp_CosDensity(HeckInter_rxnfp + HeckIntra_rxnfp), color="dodgerblue", label="Heck", fill=True)
sns.kdeplot(fp_CosDensity(BH_rxnfp), color=[255/255, 208/255, 111/255], label="Buchwald HTE", fill=True)
sns.kdeplot(fp_CosDensity(Suzuki_rxnfp), color=[231/255, 98/255, 84/255], label="Suzuki HTE", fill=True)

plt.legend()
plt.xlim([0, 1])
plt.xlabel("Cosine Similarity", fontsize=11)
plt.ylabel("Density", fontsize=11)
plt.yticks([])
plt.title("Cosine Similarity Distribution(RXNFP)", fontsize=14)
plt.tight_layout()
plt.savefig("../../figures/RXNFP Cosine Similarity density.png")
# plt.show()

# DRFP
plt.figure(dpi=500)

sns.kdeplot(fp_CosDensity(HeckInter_drfp), color=[55/255, 103/255, 149/255], label="Heck Intermolecular", fill=True)
sns.kdeplot(fp_CosDensity(HeckIntra_drfp), color=[114/255, 188/255, 213/255], label="Heck Intramolecular", fill=True)
sns.kdeplot(fp_CosDensity(HeckInter_drfp + HeckIntra_drfp), color="dodgerblue", label="Heck", fill=True)
sns.kdeplot(fp_CosDensity(BH_drfp), color=[255/255, 208/255, 111/255], label="Buchwald HTE", fill=True)
sns.kdeplot(fp_CosDensity(Suzuki_drfp), color=[231/255, 98/255, 84/255], label="Suzuki HTE", fill=True)

plt.legend()
plt.xlim([0, 1])
plt.xlabel("Cosine Similarity", fontsize=11)
plt.ylabel("Density", fontsize=11)
plt.yticks([])
plt.title("Cosine Similarity Distribution(DRFP)", fontsize=14)
plt.tight_layout()
plt.savefig("../../figures/DRFP Cosine Similarity density.png")
# plt.show()


# 10.analysis of database quality
# Attention: reagents & catalysts are all called catalysts
top_n = 15
# intra
intra_cat, intra_sol = top_coverage(rxn_intra, n=top_n)
# inter
inter_cat, inter_sol = top_coverage(rxn_inter, n=top_n)
# total
total_cat, total_sol = top_coverage(rxn_list, n=top_n)
# Figure
x_label = np.linspace(1, top_n, top_n)
fig, ax = plt.subplots(ncols=3, dpi=500, figsize=(10, 5))
ax[0].plot(x_label, intra_cat, color=[236/255, 164/255, 124/255], marker="^")
ax[0].plot(x_label, intra_sol, color=[117/255, 157/255, 219/255], marker="s")
ax[1].plot(x_label, inter_cat, color=[236/255, 164/255, 124/255], marker="^")
ax[1].plot(x_label, inter_sol, color=[117/255, 157/255, 219/255], marker="s")
ax[2].plot(x_label, total_cat, color=[236/255, 164/255, 124/255], marker="^")
ax[2].plot(x_label, total_sol, color=[117/255, 157/255, 219/255], marker="s")
# beautify
fig.suptitle("Reagents Diversity", fontsize=16)
ax[0].set_title("Intramolecular", fontsize=14)
ax[1].set_title("Intermolecular", fontsize=14)
ax[2].set_title("Total", fontsize=14)
for i in range(3):
    ax[i].legend(["catalysts", "solvents"], loc="upper left", prop={'size': 10})
    ax[i].set_xlabel("Top N", fontsize=12)
    ax[i].set_ylabel("Reaction Coverage", fontsize=12)
    ax[i].grid(True, alpha=0.5, linestyle="--")
plt.tight_layout()
plt.savefig("../../figures/Reagents Diversity.png")
# plt.show()


# 11.t-SNE for Reactants, Reagents, Reactants + Reagents
# rxnfp
# Intra
# Reaction
tsne = TSNE(n_components=2, n_iter=250, min_grad_norm=1e-6).fit(np.array(HeckIntra_rxnfp))
embeded_tsne = tsne.fit_transform(np.array(HeckIntra_rxnfp))
tsne1 = pd.DataFrame(embeded_tsne[:, 0], columns=["Intra_rxnfp_TSNE1"])
tsne2 = pd.DataFrame(embeded_tsne[:, 1], columns=["Intra_rxnfp_TSNE2"])
df_intra = pd.concat([df_intra, tsne1, tsne2], axis=1)

# Only Reactants and products
tsne = TSNE(n_components=2, n_iter=250, min_grad_norm=1e-6).fit(np.array(HeckIntra_react_rxnfp))
embeded_tsne = tsne.fit_transform(np.array(HeckIntra_react_rxnfp))
tsne1 = pd.DataFrame(embeded_tsne[:, 0], columns=["Intra_react_rxnfp_TSNE1"])
tsne2 = pd.DataFrame(embeded_tsne[:, 1], columns=["Intra_react_rxnfp_TSNE2"])
df_intra = pd.concat([df_intra, tsne1, tsne2], axis=1)

# Only Reagent
tsne = TSNE(n_components=2, n_iter=250, min_grad_norm=1e-6).fit(np.array(HeckIntra_reagent_rxnfp))
embeded_tsne = tsne.fit_transform(np.array(HeckIntra_reagent_rxnfp))
tsne1 = pd.DataFrame(embeded_tsne[:, 0], columns=["Intra_reagent_rxnfp_TSNE1"])
tsne2 = pd.DataFrame(embeded_tsne[:, 1], columns=["Intra_reagent_rxnfp_TSNE2"])
df_intra = pd.concat([df_intra, tsne1, tsne2], axis=1)

# Inter
# Reaction
tsne = TSNE(n_components=2, n_iter=250, min_grad_norm=1e-6).fit(np.array(HeckInter_rxnfp))
embeded_tsne = tsne.fit_transform(np.array(HeckInter_rxnfp))
tsne1 = pd.DataFrame(embeded_tsne[:, 0], columns=["Inter_rxnfp_TSNE1"])
tsne2 = pd.DataFrame(embeded_tsne[:, 1], columns=["Inter_rxnfp_TSNE2"])
df_inter = pd.concat([df_inter, tsne1, tsne2], axis=1)

# Only Reactants and products
tsne = TSNE(n_components=2, n_iter=250, min_grad_norm=1e-6).fit(np.array(HeckInter_react_rxnfp))
embeded_tsne = tsne.fit_transform(np.array(HeckInter_react_rxnfp))
tsne1 = pd.DataFrame(embeded_tsne[:, 0], columns=["Inter_react_rxnfp_TSNE1"])
tsne2 = pd.DataFrame(embeded_tsne[:, 1], columns=["Inter_react_rxnfp_TSNE2"])
df_inter = pd.concat([df_inter, tsne1, tsne2], axis=1)

# Only Reagent
tsne = TSNE(n_components=2, n_iter=250, min_grad_norm=1e-6).fit(np.array(HeckInter_reagent_rxnfp))
embeded_tsne = tsne.fit_transform(np.array(HeckInter_reagent_rxnfp))
tsne1 = pd.DataFrame(embeded_tsne[:, 0], columns=["Inter_reagent_rxnfp_TSNE1"])
tsne2 = pd.DataFrame(embeded_tsne[:, 1], columns=["Inter_reagent_rxnfp_TSNE2"])
df_inter = pd.concat([df_inter, tsne1, tsne2], axis=1)

# drfp
# Intra
# Reaction
tsne = TSNE(n_components=2, n_iter=250, min_grad_norm=1e-6).fit(np.array(HeckIntra_drfp))
embeded_tsne = tsne.fit_transform(np.array(HeckIntra_drfp))
tsne1 = pd.DataFrame(embeded_tsne[:, 0], columns=["Intra_drfp_TSNE1"])
tsne2 = pd.DataFrame(embeded_tsne[:, 1], columns=["Intra_drfp_TSNE2"])
df_intra = pd.concat([df_intra, tsne1, tsne2], axis=1)

# Only Reactants and products
tsne = TSNE(n_components=2, n_iter=250, min_grad_norm=1e-6).fit(np.array(HeckIntra_react_drfp))
embeded_tsne = tsne.fit_transform(np.array(HeckIntra_react_drfp))
tsne1 = pd.DataFrame(embeded_tsne[:, 0], columns=["Intra_react_drfp_TSNE1"])
tsne2 = pd.DataFrame(embeded_tsne[:, 1], columns=["Intra_react_drfp_TSNE2"])
df_intra = pd.concat([df_intra, tsne1, tsne2], axis=1)

# Only Reagent
tsne = TSNE(n_components=2, n_iter=250, min_grad_norm=1e-6).fit(np.array(HeckIntra_reagent_drfp))
embeded_tsne = tsne.fit_transform(np.array(HeckIntra_reagent_drfp))
tsne1 = pd.DataFrame(embeded_tsne[:, 0], columns=["Intra_reagent_drfp_TSNE1"])
tsne2 = pd.DataFrame(embeded_tsne[:, 1], columns=["Intra_reagent_drfp_TSNE2"])
df_intra = pd.concat([df_intra, tsne1, tsne2], axis=1)

# Inter
# Reaction
tsne = TSNE(n_components=2, n_iter=250, min_grad_norm=1e-6).fit(np.array(HeckInter_drfp))
embeded_tsne = tsne.fit_transform(np.array(HeckInter_drfp))
tsne1 = pd.DataFrame(embeded_tsne[:, 0], columns=["Inter_drfp_TSNE1"])
tsne2 = pd.DataFrame(embeded_tsne[:, 1], columns=["Inter_drfp_TSNE2"])
df_inter = pd.concat([df_inter, tsne1, tsne2], axis=1)

# Only Reactants and products
tsne = TSNE(n_components=2, n_iter=250, min_grad_norm=1e-6).fit(np.array(HeckInter_react_drfp))
embeded_tsne = tsne.fit_transform(np.array(HeckInter_react_drfp))
tsne1 = pd.DataFrame(embeded_tsne[:, 0], columns=["Inter_react_drfp_TSNE1"])
tsne2 = pd.DataFrame(embeded_tsne[:, 1], columns=["Inter_react_drfp_TSNE2"])
df_inter = pd.concat([df_inter, tsne1, tsne2], axis=1)

# Only Reagent
tsne = TSNE(n_components=2, n_iter=250, min_grad_norm=1e-6).fit(np.array(HeckInter_reagent_drfp))
embeded_tsne = tsne.fit_transform(np.array(HeckInter_reagent_drfp))
tsne1 = pd.DataFrame(embeded_tsne[:, 0], columns=["Inter_reagent_drfp_TSNE1"])
tsne2 = pd.DataFrame(embeded_tsne[:, 1], columns=["Inter_reagent_drfp_TSNE2"])
df_inter = pd.concat([df_inter, tsne1, tsne2], axis=1)

# Generate DR data
df_intra.to_excel("../../data/Heck/Intramolecular data DR.xlsx")
df_inter.to_excel("../../data/Heck/Intermolecular data DR.xlsx")
data_BH.to_excel("../../data/BH_HTE/BH_HTE DR.xlsx")
data_Suzuki.to_excel("../../data/Suzuki_HTE/Suzuki_HTE DR.xlsx")