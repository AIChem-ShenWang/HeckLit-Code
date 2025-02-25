import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

n_basic = 16
# Color
data_color={"Intra":[114/255, 188/255, 213/255],
            "Inter":[55/255, 103/255, 149/255],
            "Total":'dodgerblue',
            "BH":[255/255, 208/255, 111/255],
            "Suzuki":[231/255, 98/255, 84/255]}
color_num = ["Intra", "Inter", "BH", "Suzuki"]

# Exp Addtiion
exp_name = ["Heck Intramolecular", "Heck Intermolecular", "Buchwald HTE", "Suzuki HTE"]
# Import data
df = pd.read_excel("../../data/Heck/Heck Exp Addition n=%d.xlsx" % n_basic, sheet_name=None)
df_list = [
    df["Intra"],
    df["Inter"],
    df["BH"],
    df["Suzuki"]
]
data_list = list() # meta:[[avg], [max], [min]]
x = list() # x-axis
for i in range(len(df_list)):
    df = df_list[i]
    meta = list()

    x = np.array(df.columns[1:]).astype(np.int64)
    avg = np.array(df.iloc[5, 1:])
    maxx = list()
    minn = list()
    for j in range(1, 6):
        maxx.append(np.array(df.iloc[:5, j]).max())
        minn.append(np.array(df.iloc[:5, j]).min())
    data_list.append([avg, maxx, minn])


# Figure
plt.figure(dpi=500)
x_axis = np.linspace(1, len(x), len(x))
for i in range(len(data_list)):
    plt.fill_between(x_axis, data_list[i][1], data_list[i][2], color=data_color["%s" % color_num[i]], alpha=0.2)
    plt.plot(x_axis, data_list[i][0], color=data_color["%s" % color_num[i]])
    plt.scatter(x_axis, data_list[i][0], color=data_color["%s" % color_num[i]], label=exp_name[i])
plt.plot([1, len(x)], [data_list[0][0][0], data_list[0][0][0]], linestyle="--", color="k", label="Reference Line")


plt.grid(True, linestyle="--", alpha=0.5)
plt.legend(loc="upper left", prop={'size': 11})
plt.xlabel("Number of Data", fontsize=11)
plt.ylabel("RMSE", fontsize=11)
plt.xticks(x_axis, x)
plt.title("number of Heck Intramolecular=%d" % n_basic, fontsize=16)
plt.tight_layout()
plt.savefig("../../figures/Addition Test n=%d.png" % n_basic)
# plt.show()