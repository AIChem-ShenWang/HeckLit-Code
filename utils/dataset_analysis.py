"""
This Module contains some tools to analyse the data:
    1.Reaction feature analysis
    2.Reaction label analysis
"""
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from utils.rxn import *
from tqdm import tqdm

# 1.Reaction feature analysis
from sklearn.metrics.pairwise import euclidean_distances
def fp_mds(rxnfp_list, n):
    """
    :param morgan_list: morgan fp list
    :param n: Dimension after dimensionality reduction
    :return: result of mds
    """
    # 1.build dist matrix
    num = len(rxnfp_list)
    D = euclidean_distances(np.array(rxnfp_list))
    print("[MDS]: Finished generating distance matrix")
    # 2.centralization dist matrix
    I = np.eye(num)
    i = np.ones([num, num])
    H = I - (1 / num) * i
    # 3.calculate inner product matrix
    B = -1 / 2 * np.dot(np.dot(H, D ** 2), H)
    # 4.eigenvalue decomposition
    Be, Bv = np.linalg.eigh(B)  # eigenvalues of Be matrix B, normalized eigenvectors of Bv
    simple_vector = Bv[:, np.argsort(-Be)[:n]]
    simple_val = Be[np.argsort(-Be)[:n]]
    Z = simple_vector * simple_val ** 0.5
    return Z

from sklearn.metrics.pairwise import cosine_similarity
def fp_CosDensity(rxn_list):
    D = cosine_similarity(np.array(rxn_list))
    S = list()
    for i in range(D.shape[0]):
        for j in range(D.shape[1]):
            if i < j:
                S.append(D[i][j])
    return np.array(S)

def top_coverage(rxn_list, n):
    """
    :param rxn_list: rxn_class list
    :param n: top n in number
    :return: 3 dicts, shows the coverage rate of top n substances
    """
    rxn_cat_dict = dict()
    rxn_sol_dict = dict()
    num_cat = 0
    num_sol = 0
    for rxn in rxn_list:
        # reagents
        for reag in rxn.reagents:
            num_cat += 1
            if reag not in rxn_cat_dict:
                rxn_cat_dict["%s" % reag] = 1
            else:
                rxn_cat_dict["%s" % reag] += 1
        # catalysts
        for cat in rxn.cats:
            num_cat += 1
            if cat not in rxn_cat_dict:
                rxn_cat_dict["%s" % cat] = 1
            else:
                rxn_cat_dict["%s" % cat] += 1
        # solvents
        for sol in rxn.solvents:
            num_sol += 1
            if sol not in rxn_sol_dict:
                rxn_sol_dict["%s" % sol] = 1
            else:
                rxn_sol_dict["%s" % sol] += 1
    # number of top n
    rxn_cat_dict = list(zip(rxn_cat_dict.values(), rxn_cat_dict.keys()))
    rxn_cat_dict = sorted(rxn_cat_dict, reverse=True)
    rxn_sol_dict = list(zip(rxn_sol_dict.values(), rxn_sol_dict.keys()))
    rxn_sol_dict = sorted(rxn_sol_dict, reverse=True)
    # return coverage rate list
    rxn_cat_list = list()
    rxn_sol_list = list()
    for i in range(n):
        if i == 0:
            rxn_cat_list.append(rxn_cat_dict[0:n][i][0] / num_cat)
            rxn_sol_list.append(rxn_sol_dict[0:n][i][0] / num_sol)
        else:
            rxn_cat_list.append(rxn_cat_list[i-1] + rxn_cat_dict[0:n][i][0] / num_cat)
            rxn_sol_list.append(rxn_sol_list[i-1] + rxn_sol_dict[0:n][i][0] / num_sol)

    return rxn_cat_list, rxn_sol_list


# 2.Reaction label analysis
# classify the reaction by yield distribution
def shot_classifier(yield_dis, split_num, upper, lower):
    yield_dict = dict()

    # Initialize
    interval = 100 / split_num
    for i in range(split_num):
        # interval * i~interval interval * (i+1)
        yield_dict[interval * i] = 0

    for y in yield_dis:
        for key in yield_dict.keys():
            if y >= key and y < key + interval:
                yield_dict[key] += 1

    # give the shot-class
    for key in yield_dict.keys():
        num = yield_dict[key]
        # Few-shot
        if num < lower:
            yield_dict[key] = [num, "Few-shot"]
            continue

        # Many-shot
        if num > upper:
            yield_dict[key] = [num, "Many-shot"]
            continue

        # Medium-shot
        if num <= upper and num >= lower:
            yield_dict[key] = [num, "Medium-shot"]

    return yield_dict
