import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from utils.rxn import *
from utils.dataset_analysis import *

# Import data
data_JCP = pd.read_excel("../data/Heck/JCP processed data.xlsx")
DRFP_JCP = [read_drfp(data_JCP.loc[i]["drfp"]) for i in range(data_JCP.shape[0])]
data_HeckLit = pd.read_excel("../data/Heck/Heck_fp.xlsx")
DRFP_HeckLit = [read_drfp(data_HeckLit.loc[i]["drfp"]) for i in range(data_HeckLit.shape[0])]

print("Calculating Start")
plt.figure(dpi=500)

sns.kdeplot(fp_CosDensity(DRFP_JCP), color="#24a645", label="Das et al.", fill=True)
sns.kdeplot(fp_CosDensity(DRFP_HeckLit), color="dodgerblue", label="HeckLit", fill=True)

plt.legend(prop={'size': 16})
plt.xlim([0, 1])
plt.xlabel("Cosine Similarity", fontsize=18)
plt.ylabel("Density", fontsize=18)
plt.yticks([])
plt.title("Cosine Similarity Distribution(DRFP)", fontsize=20)
plt.tight_layout()
plt.savefig("../figures/Cosine Similarity density Das&HeckLit.png")
# plt.show()
