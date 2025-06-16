import re
import numpy as np
import pandas as pd
from utils.rxn import *
from utils.molecule import *
from gensim.models import Word2Vec
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns

# mission_list = ["HeckLit", "JCP"]
mission_list = ["JCP"]

# size of word vector
size = 128

for mission in mission_list:
    train_set = list()
    vocab_dict = dict()
    seq_len = list()
    if mission == "HeckLit":
        Heck_df = pd.read_excel("../data/Heck/Heck processed data.xlsx")
        rxn_list = df_to_rxn_list(Heck_df)
        for rxn in rxn_list:
            text = get_Heck_RxnSmi(rxn)
            train_set.append(smi_tokenizer(text))
            seq_len.append(len(smi_tokenizer(text)))
        print("Maximum sequence length is:%d" % max(np.array(seq_len)))

    if mission == "JCP":
        df = pd.read_excel("../data/Heck/JCP processed data.xlsx")
        for i in range(df.shape[0]):
            text = get_JCP_Heck_RxnSmi(df.iloc[i, :])
            train_set.append(smi_tokenizer(text))
            seq_len.append(len(smi_tokenizer(text)))
        print("Maximum sequence length is:%d" % max(np.array(seq_len)))

    # sequence length distribution
    plt.figure(dpi=500)
    sns.kdeplot(seq_len, fill=True)
    plt.xlabel("Sequence Length")
    plt.ylabel("Counting")
    plt.yticks([])
    plt.title("Sequence Distribution")
    plt.tight_layout()
    plt.savefig("../figures/%s Sequence Distribution.png" % mission)

    # generate vocab file
    print("The length of train set is: %d" % len(train_set))
    word_id = dict()
    word_vec = list()
    if mission == "HeckLit":
        model = Word2Vec(train_set, vector_size=size, window=50, min_count=5, epochs=30, sg=1)
    if mission == "JCP":
        model = Word2Vec(train_set, vector_size=size, window=30, min_count=3, epochs=50, sg=1)

    for i, w in enumerate(model.wv.index_to_key):
        vocab_dict[w] = model.wv[w]
        # record for evaluation
        word_id[w] = i
        word_vec.append(model.wv[w])
    vocab_dict_to_txt(vocab_dict, rxn_name="%s" % mission)
    print("There are %d atoms in the dict" % len(vocab_dict))

    # Evaluation
    X_reduced = PCA(n_components=2).fit_transform(np.array(word_vec))
    plt.figure(dpi=500)
    plt.scatter(X_reduced[:, 0], X_reduced[:, 1], color="black")

    for w in word_id.keys():
        xy = X_reduced[word_id[w], :]
        plt.scatter(xy[0], xy[1], color="r")
        plt.text(xy[0], xy[1], w, color="b")

    plt.title("%s Vocab Distribution" % mission)
    plt.tight_layout()
    plt.savefig("../figures/%s Vocab Distribution.png" % mission)

