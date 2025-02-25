"""
This module includes some operations towards rxns:
    1.Extraction of rxn data
    2.Counting and Filter the rxns
    3.Convert RXN SMILES into RXNFP or DRFP
"""
from utils.molecule import *

# 1.Extraction of rxn data
import re
import pandas as pd
import torch

class RXN():
    # Prameters of RXN
    def __init__(self):
        self.index = str() # Index

        self.reactants = list() # reactants
        self.products = list() # products

        self.reagents = list() # regents
        self.cats = list() # catalysts
        self.solvents = list() # solvents
        self.temp = str() # temperature
        self.time = str() # time

        self.rxn_yield = str() # yield

        self.ref = str() # reference
        self.rxn_id = str() # Reaxys ID

    # get the information of a rxn using Python Spider
    def get_info(self, html_parser, i):
        """
        :param html_parser: html parser
        :param i: the i-th condition in the rxn
        """
        # rxn index
        self.index = html_parser.xpath("//span[@class='rx-element-index']/text()")[0]

        rxn_index = html_parser.xpath("//div[@class='rx-reactions-table__conditions__steps']")[i]
        # if there is no stages-row
        rxn_index = rxn_index.xpath('./div[@class="stages-row"]')
        if len(rxn_index) == 0:
            return
        else:
            rxn_index = rxn_index[0]
        # reagent
        rxn_r = rxn_index.xpath('./span[@class="stage-reagents"]')
        if len(rxn_r) == 0:
            pass
        else:
            for rea_html in rxn_r:
                rea = str(rea_html.xpath("string()"))
                if len(re.findall(".+[^;\xa0]", rea)) == 0:
                    self.reagents.append(rea)
                else:
                    rea = re.findall(".+[^;\xa0]", rea)[0] # html文本中有奇怪字符
                    self.reagents.append(rea)
        # catalyst
        rxn_c = rxn_index.xpath("./span[@class='stage-catalyst']")
        if len(rxn_c) == 0:
            pass
        else:
            for cat_html in rxn_c:
                cat = str(cat_html.xpath("string()"))
                if len(re.findall(".+[^;\xa0]", cat)) == 0:
                    self.cats.append(cat)
                else:
                    cat = re.findall(".+[^;\xa0]", cat)[0] # if there is strange notes in html
                    self.cats.append(cat)
        # solvent
        rxn_s = rxn_index.xpath("./span[@class='stage-solvents']")
        if len(rxn_s) == 0:
            pass
        else:
            for sol_html in rxn_s:
                sol = str(sol_html.xpath("string()"))
                if len(re.findall(".+[^;\xa0]", sol)) == 0:
                    self.solvents.append(sol)
                else:
                    sol = re.findall(".+[^;\xa0]", sol)[0] # if there is strange notes in html
                    self.solvents.append(sol)
        # time & temperature
        rxn_cond = rxn_index.xpath("string(./span[@class='conditions'])")
        temp = re.findall("(?<=at\s)-?[0-9]+\.[0-9]*|-?[0-9]+(?=℃)", rxn_cond)
        time = re.findall("(?<=for\s)-?[0-9]+\.[0-9]*|-?[0-9]+(?=h)", rxn_cond)
        if len(temp) == 0:
            self.temp = "/"
        else:
            self.temp = temp[0]

        if len(time) == 0:
            self.time = "/"
        else:
            self.time = time[0]

        # yield
        rxn_y = html_parser.xpath("//td[@class='rx-reactions-table__yield display-table-cell']")[i]
        if len(rxn_y) == 0:
            self.rxn_yield = "/"
        else:
            self.rxn_yield = rxn_y.xpath("string()")

        # reference
        rxn_ref = html_parser.xpath("//div[@class='citation clear']")[i]
        if len(rxn_ref) == 0:
            self.ref = "/"
        else:
            self.ref = rxn_ref.xpath("string()")

    def show_info(self):
        print(self.index, self.reactants, self.products, self.reagents, self.cats, self.solvents, self.temp, self.time, self.rxn_yield, self.rxn_id, self.ref)

# ConvertDataframe into list
def df_to_rxn_list(df):
    """
    :param df: pandas的dataframe
    :return: the rxn_class for all rxns, list type
    """

    rxn_list = list()  # list
    data_size = df.shape[0]  # the counting of the rxn

    for num in range(data_size):
        rxn = RXN() # rxn class
        # get index
        rxn.index = df.loc[num]["rxn_index"]
        # get temperature /C
        rxn.temp = df.loc[num]["temperature /C"]
        # get time /h
        rxn.time = df.loc[num]["time /h"]
        # get Yield
        rxn.rxn_yield = df.loc[num]["Yield"]
        # get Reaction ID
        rxn.rxn_id = df.loc[num]["Reaction ID"]
        # get reference
        rxn.ref = df.loc[num]["Reference"]

        # get reactant
        for col in df.columns:
            if "reactants" in col:
                if df.loc[num][col] != "/": # remove /
                    rxn.reactants.append(df.loc[num][col])
        # get product
        for col in df.columns:
            if "products" in col:
                if df.loc[num][col] != "/":  # remove /
                    rxn.products.append(df.loc[num][col])
        # get reagent
        for col in df.columns:
            if "reagents" in col:
                if df.loc[num][col] != "/":  # remove /
                    rxn.reagents.append(df.loc[num][col])
        # get catalyst
        for col in df.columns:
            if "catalysts" in col:
                if df.loc[num][col] != "/":  # remove /
                    rxn.cats.append(df.loc[num][col])
        # get solvents
        for col in df.columns:
            if "solvents" in col:
                if df.loc[num][col] != "/":  # remove /
                    rxn.solvents.append(df.loc[num][col])

        rxn_list.append(rxn) #  put into list

    return rxn_list

# Convert rxn_list into Dataframe
def rxn_list_to_df(rxn_list):
# search rxn_list, find the maximum number of reactants, products to determine the size of Dataframe
    num_reactants = 0
    num_products = 0
    num_reagents = 0
    num_cats = 0
    num_sol = 0

    for rxn in rxn_list:
        num_reactants = max(len(rxn.reactants), num_reactants)
        num_products = max(len(rxn.products), num_products)
        num_reagents = max(len(rxn.reagents), num_reagents)
        num_cats = max(len(rxn.cats), num_cats)
        num_sol = max(len(rxn.solvents), num_sol)

    # set the columns of dataframe
    cols = list()
    # rxnn index
    cols.append("rxn_index")
    # reactants
    for i in range(num_reactants):
        cols.append("reactants %d" % (i + 1))
    # products
    for i in range(num_products):
        cols.append("products %d" % (i + 1))
    # reagents
    for i in range(num_reagents):
        cols.append("reagents %d" % (i + 1))
    # catalysts
    for i in range(num_cats):
        cols.append("catalysts %d" % (i + 1))
    # solvents
    for i in range(num_sol):
        cols.append("solvents %d" % (i + 1))
    # temperature
    cols.append("temperature /C")
    # time
    cols.append("time /h")
    # yield
    cols.append("Yield")
    # Reaction ID
    cols.append("Reaction ID")
    # Referenc
    cols.append("Reference")

    # input data to dataframe
    data = list()
    for rxn in rxn_list:
        meta_data = list()
        # index
        meta_data.append(rxn.index)
        # reactants
        for reactant in rxn.reactants:
            meta_data.append(reactant)
        while len(rxn.reactants) < num_reactants:
            meta_data.append("/")
            rxn.reactants.append("/")
        # products
        for product in rxn.products:
            meta_data.append(product)
        while len(rxn.products) < num_products:
            meta_data.append("/")
            rxn.products.append("/")
        # reagents
        for reagent in rxn.reagents:
            meta_data.append(reagent)
        while len(rxn.reagents) < num_reagents:
            meta_data.append("/")
            rxn.reagents.append("/")
        # catalysts
        for cat in rxn.cats:
            meta_data.append(cat)
        while len(rxn.cats) < num_cats:
            meta_data.append("/")
            rxn.cats.append("/")
        # solvents
        for sol in rxn.solvents:
            meta_data.append(sol)
        while len(rxn.solvents) < num_sol:
            meta_data.append("/")
            rxn.solvents.append("/")
        # temperture
        meta_data.append(rxn.temp)
        # time
        meta_data.append(rxn.time)
        # yield
        meta_data.append(rxn.rxn_yield)
        # Reaction ID
        meta_data.append(rxn.rxn_id)
        # Reference
        meta_data.append(rxn.ref)

        # put the rxn to rxn_list
        data.append(meta_data)

    # produce dataframe
    df = pd.DataFrame(data, columns=cols)
    return df



# 2.Counting and Filter the rxns
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
from tqdm import tqdm
from collections import Counter

# Find the insertion(addition) type for C=C
def IntraInsertType(reactant_mol, product_mol):

    # Match CHR=CH2 to check wether the reaction is regio
    rea_EneNum = len(reactant_mol.GetSubstructMatches(Chem.MolFromSmarts("[CH1]=[CH2]")))
    prod_EneNum = len(product_mol.GetSubstructMatches(Chem.MolFromSmarts("[CH1]=[CH2]")))

    # Neither Alpha nor Beta
    if rea_EneNum == 0 or rea_EneNum == prod_EneNum:
        return "*"

    # Match CR2=CH2 & CHR=CHR
    rea_CH1Num = len(reactant_mol.GetSubstructMatches(Chem.MolFromSmarts("[CH1]")))
    prod_CH1Num = len(product_mol.GetSubstructMatches(Chem.MolFromSmarts("[CH1]")))

    # Alpha
    if rea_CH1Num - prod_CH1Num == 1:
        return "Alpha"

    # Beta
    if rea_CH1Num - prod_CH1Num == -1:
        return "Beta"

# Find LG for Heck Intramolecular rxn
def CntAtomNum(mol):
    Atom_Num = list()
    if mol is not None:
        # 使用列表推导来获取所有原子的符号
        atom_symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]

        # 使用Counter来统计每种原子的数量
        atom_counts = Counter(atom_symbols)

        # 输出结果
        for atom, count in atom_counts.items():
            Atom_Num.append([atom, count])
        return Atom_Num

def IntraLGType(reactant_mol, product_mol):
    reactantAtoms= CntAtomNum(reactant_mol)
    productAtoms = CntAtomNum(product_mol)
    difference = [item for item in reactantAtoms if item not in productAtoms]
    if len(difference) == 0:
        return "*"
    return difference[0][0] # list type, list inside, i.e. [[], [], ...]

# The member of the ring forms in the rxn
def CntChangeRing(reactant_mol, product_mol):
    Change_list = list()  # The member size of changed ring

    ReaSSSR = Chem.GetSymmSSSR(reactant_mol)
    ReaRingHash = [0] * 10
    ProdSSSR = Chem.GetSymmSSSR(product_mol)
    ProdRingHash = [0] * 10
    for ring in ReaSSSR:
        ReaRingHash[len(list(ring))] += 1
    for ring in ProdSSSR:
        ProdRingHash[len(list(ring))] += 1
    # Get Change
    for i in range(10):
        diff = ProdRingHash[i] - ReaRingHash[i]
        for j in range(diff):
            Change_list.append(i)

    return len(Change_list), Change_list  # return the Number of changed ring and its member size

def HeckIntra_Classify(rxn):
    """
        This function will recognize:
            1.The addition Type of Alkene in Heck rxn - "Alpha" or" Beta" or "*"(Other)
            2.The type of Leaving group in ArylHalide - "F", "Cl", "Br", "I", "*"(Other)
            3.The member of the ring forms in the rxn
    """

    reactant_mol = Chem.MolFromSmiles(rxn.reactants[0])
    product_mol = Chem.MolFromSmiles(rxn.products[0])

    # 1.Addition Type
    Type_Add = IntraInsertType(reactant_mol, product_mol)

    # 2.Leaving-Group
    Type_LG = IntraLGType(reactant_mol, product_mol)
    X = ["F", "Cl", "Br", "I"]
    if Type_LG not in X:
        Type_LG = "*"

    # 3.The member of the ring forms in the rxn
    ChangeNum, MemberSize = CntChangeRing(reactant_mol, product_mol)

    return Type_Add, Type_LG, MemberSize

def IntraCount(rxn_list, file):
    Type = dict()
    reactant_dict = dict()
    AlphaInsert = 0
    BetaInsert = 0
    X = ["F", "Cl", "Br", "I", "*"]
    X_Num = [0] * 5
    RingType = [0] * 12

    for i in tqdm(range(len(rxn_list))):
        rxn = rxn_list[i]

        # Count Types
        if rxn.rxn_id not in Type.keys():
            Type[rxn.rxn_id] = 1

        # Count Reactants
        reactant_dict[rxn.reactants[0]] = 1

        # Reaction Type
        Type_Add, Type_LG, MemberSize = HeckIntra_Classify(rxn)
        # Insertion Type
        if Type_Add == "Alpha":
            AlphaInsert += 1
        if Type_Add == "Beta":
            BetaInsert += 1

        # LG type
        for i in range(len(X[:-1])):
            if X[i] == Type_LG:
                X_Num[i] += 1

        # RingSize
        RingType[MemberSize[0]] += 1


    print("[Counting of Heck Intramoleculer Reaction]")
    print("types: %d cases: %d" % (len(Type), len(rxn_list)))
    print("Reactants Num:%d" % len(reactant_dict))
    print("AlphaInsertion Num:%d" % AlphaInsert)
    print("BetaInsertion Num:%d" % BetaInsert)
    X_Num[-1] = len(rxn_list) - np.array(X_Num[:-1]).sum()
    for i in range(len(X)):
        print("LG=%s Num:%d" % (X[i], X_Num[i]))
    for i in range(10):
        if RingType[i] != 0:
            print("RingType:%d-membered Num = %d" % (i, RingType[i]))

    file.write("[Counting of Heck Intramoleculer Reaction]\n")
    file.write("types: %d cases: %d\n" % (len(Type), len(rxn_list)))
    file.write("Reactants Num:%d\n" % len(reactant_dict))
    file.write("AlphaInsertion Num:%d\n" % AlphaInsert)
    file.write("BetaInsertion Num:%d\n" % BetaInsert)
    for i in range(len(X)):
        file.write("LG=%s Num:%d\n" % (X[i], X_Num[i]))
    file.write("\n")
    for i in range(12):
        if RingType[i] != 0:
            file.write("RingType:%d-membered Num = %d\n" % (i, RingType[i]))
    file.write("\n")

    return {"Alpha-Insertion":AlphaInsert,
            "Beta-Insertion":BetaInsert,
            "LG = F":X_Num[0],
            "LG = Cl":X_Num[1],
            "LG = Br":X_Num[2],
            "LG = I":X_Num[3],
            "5-membered ring":RingType[5],
            "6-membered ring":RingType[6],
            "7-membered ring":RingType[7],
            "8-membered ring":RingType[8],
            "9-membered ring":RingType[9],
            }

def HeckInter_Classify(rxn):
    """
        This function will recognize:
            1.The addition Type of Alkene in Heck rxn - "Alpha" or" Beta" or "*"(Other)
            2.The type of Leaving group in ArylHalide - "F", "Cl", "Br", "I", "*"(Other)
    """

    reactants_mol = [Chem.MolFromSmiles(i) for i in rxn.reactants]
    product_mol = Chem.MolFromSmiles(rxn.products[0])

    # 1.Addition Type
    Type_Add = "*"

    # Alpha
    rxnTemp = AllChem.ReactionFromSmarts(
        '[C,c;H1:1](=,:[C,c;H2:2]).[C,c:3]-[*]>>[C,c;H0:1](=,:[C,c;H2:2])[C,c:3]')
    ps = rxnTemp.RunReactants(reactants_mol)
    for p in ps:
        if Chem.MolToSmiles(p[0]) == Chem.MolToSmiles(product_mol) or \
                Chem.MolToSmiles(p[0]).replace("-", "") == Chem.MolToSmiles(product_mol):
            Type_Add = "Alpha"

    # Beta
    rxnTemp = AllChem.ReactionFromSmarts(
        '[C,c;H1:1](=,:[C,c;H2:2]).[C,c:3]-[*]>>[C,c;H1:1]=,:[C,c;H1:2]-[C,c:3]')
    ps = rxnTemp.RunReactants(reactants_mol)
    for p in ps:
        if Chem.MolToSmiles(p[0]) == Chem.MolToSmiles(product_mol) or \
                Chem.MolToSmiles(p[0]).replace("-", "") == Chem.MolToSmiles(product_mol):
            Type_Add = "Beta"

    # 2.Leaving-Group
    X = ["F", "Cl", "Br", "I", "*"]
    Type_LG = "*"
    for i in range(len(X[:-1])):
        rxnTemp = AllChem.ReactionFromSmarts('[C,c:1](=,:[C,c:2]).[C,c:3]-[%s]>>[C,c:1]=,:[C,c:2]-[C,c:3]' % X[i])
        ps = rxnTemp.RunReactants(reactants_mol)
        for p in ps:
            if Chem.MolToSmiles(p[0]) == Chem.MolToSmiles(product_mol) or \
                    Chem.MolToSmiles(p[0]).replace("-", "") == Chem.MolToSmiles(product_mol):
                Type_LG = X[i]

    return Type_Add, Type_LG

def InterCount(rxn_list, file):
    Type = dict()
    Alkene_dict = dict()
    AlphaInsert = 0
    BetaInsert = 0
    ArylHalide_dict = dict()
    X = ["F", "Cl", "Br", "I", "*"] # LG = F/Cl/Br/I/others
    X_Num = [0] * 5

    for i in tqdm(range(len(rxn_list))):
        rxn = rxn_list[i]

        # Count Types
        if rxn.rxn_id not in Type.keys():
            Type[rxn.rxn_id] = 1

        # Count Reactants
        Alkene_dict[rxn.reactants[0]] = 1
        ArylHalide_dict[rxn.reactants[1]] = 1

        # Reaction Type
        Type_Add, Type_LG = HeckInter_Classify(rxn)
        # Addition Type
        if Type_Add == "Alpha":
            AlphaInsert += 1
        if Type_Add == "Beta":
            BetaInsert += 1

        # Leaving-Group Type
        for i in range(len(X[:-1])):
            if X[i] == Type_LG:
                X_Num[i] += 1

    print("[Counting of Heck Intermolecule Reaction]")
    print("types: %d cases: %d" % (len(Type), len(rxn_list)))
    print("Alkene Num:%d" % len(Alkene_dict))
    print("ArylHalide Num:%d" % len(ArylHalide_dict))
    print("AlphaInsertion Num:%d" % AlphaInsert)
    print("BetaInsertion Num:%d" % BetaInsert)
    X_Num[-1] = len(rxn_list) - np.array(X_Num[:-1]).sum()
    for i in range(len(X)):
        print("LG=%s Num:%d" % (X[i], X_Num[i]))

    file.write("[Counting of Heck Intermolecule Reaction]\n")
    file.write("types: %d cases: %d\n" % (len(Type), len(rxn_list)))
    file.write("Alkene Num:%d\n" % len(Alkene_dict))
    file.write("ArylHalide Num:%d\n" % len(ArylHalide_dict))
    file.write("AlphaInsertion Num:%d\n" % AlphaInsert)
    file.write("BetaInsertion Num:%d\n" % BetaInsert)
    for i in range(len(X)):
        file.write("LG=%s Num:%d\n" % (X[i], X_Num[i]))
    file.write("\n")

    return {"Alph-Insertion":AlphaInsert,
            "Beta-Insertion":BetaInsert,
            "LG = F":X_Num[0],
            "LG = Cl":X_Num[1],
            "LG = Br":X_Num[2],
            "LG = I":X_Num[3]}

def HeckCount(rxn_list, file):
    reactants_dict = dict()
    products_dict = dict()
    Pdcats_dict = dict()
    cats_dict = dict()
    solvents_dict = dict()

    for rxn in rxn_list:

        for reactant in rxn.reactants:
            reactants_dict[reactant] = 1
        for product in rxn.products:
            products_dict[product] = 1

        for cat in rxn.reagents:
            if "Pd" in cat:
                Pdcats_dict[cat] = 1
            else:
                cats_dict[cat] = 1
        for cat in rxn.cats:
            if "Pd" in cat:
                Pdcats_dict[cat] = 1
            else:
                cats_dict[cat] = 1
        for solvent in rxn.solvents:
            solvents_dict[solvent] = 1
    print("[Counting of Heck Reaction]")
    print("Reactants Num:%d" % len(reactants_dict))
    print("Products Num:%d" % len(products_dict))
    print("Catalysts with Pd Num:%d" % len(Pdcats_dict))
    print("Catalysts without Pd Num:%d" % len(cats_dict))
    print("Solvents Num:%d" % len(solvents_dict))

    file.write("[Counting of Heck Reaction]\n")
    file.write("Reactants Num:%d\n" % len(reactants_dict))
    file.write("Products Num:%d\n" % len(products_dict))
    file.write("Catalysts with Pd Num:%d\n" % len(Pdcats_dict))
    file.write("Catalysts without Pd Num:%d\n" % len(cats_dict))
    file.write("Solvents Num:%d\n" % len(solvents_dict))
    file.write("\n")

# Heck Filter
def Heck_filter1(rxn):
    # Convert into mol
    try:
        reactants_mol = [Chem.AddHs(Chem.MolFromSmiles(i)) for i in rxn.reactants]
        product_mol = Chem.AddHs(Chem.MolFromSmiles(rxn.products[0]))

        # check the change of the number of the atom
        num_reactants = np.array([i.GetNumAtoms() for i in reactants_mol]).sum()
        num_product = product_mol.GetNumAtoms()
        if num_reactants - num_product != 2: # HX
            if num_reactants - num_product != 9: # TfOH
                if num_reactants - num_product != 19: # TsOH
                    if num_reactants - num_product != 8: # AcOH
                        return False

        # In the reactants, there must be a halogen atom & a C=C group
        LG = ["F", "Cl", "Br", "I", "OS(=O)(=O)C(F)(F)F", "CC1=CC=C(C=C1)S(=O)(=O)O", "CC(=O)O"]
        patt_DoubleBond = Chem.MolFromSmarts('C=C')

        # Intramolecular rxn
        if len(rxn.reactants) == 1:
            pd = False

            # check C=C
            if reactants_mol[0].HasSubstructMatch(patt_DoubleBond) == False:
                return False

            # check LG
            for lg in LG:
                patt_lg = Chem.MolFromSmarts(lg)
                if reactants_mol[0].HasSubstructMatch(patt_lg):
                    pd = True

            # check the change number of the ring
            rea_mol = Chem.MolFromSmiles(rxn.reactants[0])
            prod_mol = Chem.MolFromSmiles(rxn.products[0])
            ChangeNum, ChangeList = CntChangeRing(rea_mol, prod_mol)
            if ChangeNum != 1:
                pd = False

            return pd

        # Intermolecular rxn
        if len(rxn.reactants) == 2:
            pd1_lg = False
            pd2_lg = False
            # C=C
            pd1_db = reactants_mol[0].HasSubstructMatch(patt_DoubleBond)
            pd2_db = reactants_mol[1].HasSubstructMatch(patt_DoubleBond)

            # X
            for lg in LG:
                patt_lg = Chem.MolFromSmarts(lg)
                if reactants_mol[0].HasSubstructMatch(patt_lg):
                    pd1_lg = True
                if reactants_mol[1].HasSubstructMatch(patt_lg):
                    pd2_lg = True

            if int(pd1_lg) + int(pd2_lg) == 0 or int(pd1_db) + int(pd2_db) == 0:
                return False

            # check template
            reactants_mol = [Chem.MolFromSmiles(i) for i in rxn.reactants]
            product_mol = Chem.MolFromSmiles(rxn.products[0])

            rxnTemp = AllChem.ReactionFromSmarts(
                '[C,c:1](=,:[C,c:2]).[C,c:3]-[*]>>[C,c:1](=,:[C,c:2])[C,c:3]')
            # situ1
            ps = rxnTemp.RunReactants(reactants_mol)
            for p in ps:
                if Chem.MolToSmiles(p[0]) == Chem.MolToSmiles(product_mol) or \
                   Chem.MolToSmiles(p[0]).replace("-", "") == Chem.MolToSmiles(product_mol):
                    return True

            # situ2
            reactants_mol.reverse()
            rxn.reactants.reverse()
            ps = rxnTemp.RunReactants(reactants_mol)
            for p in ps:
                if Chem.MolToSmiles(p[0]) == Chem.MolToSmiles(product_mol) or \
                   Chem.MolToSmiles(p[0]).replace("-", "") == Chem.MolToSmiles(product_mol):
                    return True

    except:
        return False

def Heck_filter2(rxn):
    # Check whether there is Pd in the catalyst
    pd = False
    try:
        for reagent in rxn.reagents:
            if "Pd" in str(reagent):
                pd = True
        for cat in rxn.cats:
            if "Pd" in str(cat):
                pd = True

    except:
        pd = False

    return pd


# 3. Convert RXN SMILES into RXNFP or DRFP
from rxnfp.transformer_fingerprints import *
from drfp import DrfpEncoder

# read rxnfp and drfp
def read_rxnfp(arr_str):
    arr_str = arr_str[1:]
    arr_str = arr_str[:-1]
    arr_str = arr_str.split(", ")
    for i in range(len(arr_str)):
        arr_str[i] = float(arr_str[i])
    return np.array(arr_str).astype(np.float16)

def read_drfp(arr_str):
    arr_str = arr_str[1:]
    arr_str = arr_str[:-1]
    arr_str = arr_str.split(" ")
    for i in range(len(arr_str)):
        if "\n" in arr_str[i]:
            arr_str[i] = arr_str[i][0]
        arr_str[i] = float(arr_str[i])
    return np.array(arr_str).astype(np.float16)

def get_Buchwald_RxnSmi(BH_HTE_df):
    base = str(BH_HTE_df.loc["base_smiles"])
    ligand = str(BH_HTE_df.loc["ligand_smiles"])
    aryl_halide = str(BH_HTE_df.loc["aryl_halide_smiles"])
    additive = str(BH_HTE_df.loc["additive_smiles"])
    product = str(BH_HTE_df.loc["product_smiles"])

    text = [str("CC1=CC=C(N)C=C1")] + ["."] + [aryl_halide] + [">"] + [additive] + ["."] + [base] + ["."] + [ligand] + [
        ">"] + [product]

    return "".join(text)

def get_Buchwald_rxnfp(BH_HTE_df):
    base = str(BH_HTE_df.loc["base_smiles"])
    ligand = str(BH_HTE_df.loc["ligand_smiles"])
    aryl_halide = str(BH_HTE_df.loc["aryl_halide_smiles"])
    additive = str(BH_HTE_df.loc["additive_smiles"])
    product = str(BH_HTE_df.loc["product_smiles"])

    text = [str("CC1=CC=C(N)C=C1")] + ["."] + [aryl_halide] + ["."] + [additive] + ["."] + [base] + ["."] + [ligand] + [
        ">>"] + [product]
    text = "".join(text)

    model, tokenizer = get_default_model_and_tokenizer()
    rxnfp_generator = RXNBERTFingerprintGenerator(model, tokenizer)
    rxnfp = rxnfp_generator.convert(text)

    return rxnfp

def get_Buchwald_drfp(BH_HTE_df):
    base = str(BH_HTE_df.loc["base_smiles"])
    ligand = str(BH_HTE_df.loc["ligand_smiles"])
    aryl_halide = str(BH_HTE_df.loc["aryl_halide_smiles"])
    additive = str(BH_HTE_df.loc["additive_smiles"])
    product = str(BH_HTE_df.loc["product_smiles"])

    text = [str("CC1=CC=C(N)C=C1")] + ["."] + [aryl_halide] + ["."] + [additive] + ["."] + [base] + ["."] + [ligand] + [
        ">>"] + [product]
    text = "".join(text)

    drfp = DrfpEncoder.encode(text)[0]

    return drfp

def get_Heck_RxnSmi(Heck_rxn):
    text = list()
    for reactant in Heck_rxn.reactants:
        text = text + [str(reactant)] + ["."]

    text = text[:-1] + [">"]

    for reagent in Heck_rxn.reagents:
        text = text + [str(reagent)] + ["."]

    for cat in Heck_rxn.cats:
        text = text + [str(cat)] + ["."]

    for sol in Heck_rxn.solvents:
        text = text + [str(sol)] + ["."]

    text = text[:-1] + [">"]

    text = text + [str(Heck_rxn.products[0])]

    return "".join(text)

def get_Heck_rxnfp(Heck_rxn):
    text = list()
    for reactant in Heck_rxn.reactants:
        text = text + [str(reactant)] + ["."]

    for reagent in Heck_rxn.reagents:
        text = text + [str(reagent)] + ["."]

    for cat in Heck_rxn.cats:
        text = text + [str(cat)] + ["."]

    for sol in Heck_rxn.solvents:
        text = text + [str(sol)] + ["."]

    text = text[:-1] + [">>"]

    text = text + [str(Heck_rxn.products[0])]
    text = "".join(text)

    model, tokenizer = get_default_model_and_tokenizer()
    rxnfp_generator = RXNBERTFingerprintGenerator(model, tokenizer)
    rxnfp = rxnfp_generator.convert(text)

    return rxnfp

def get_Heck_React_rxnfp(Heck_rxn):
    text = list()
    for reactant in Heck_rxn.reactants:
        text = text + [str(reactant)] + ["."]

    text = text[:-1] + [">>"]

    text = text + [str(Heck_rxn.products[0])]
    text = "".join(text)

    model, tokenizer = get_default_model_and_tokenizer()
    rxnfp_generator = RXNBERTFingerprintGenerator(model, tokenizer)
    rxnfp = rxnfp_generator.convert(text)

    return rxnfp

def get_Heck_Reagent_rxnfp(Heck_rxn):
    text = [">"]

    for reagent in Heck_rxn.reagents:
        text = text + [str(reagent)] + ["."]

    for cat in Heck_rxn.cats:
        text = text + [str(cat)] + ["."]

    for sol in Heck_rxn.solvents:
        text = text + [str(sol)] + ["."]

    text = text[:-1] + [">"]
    text = "".join(text)

    model, tokenizer = get_default_model_and_tokenizer()
    rxnfp_generator = RXNBERTFingerprintGenerator(model, tokenizer)
    rxnfp = rxnfp_generator.convert(text)

    return rxnfp

def get_Heck_drfp(Heck_rxn):
    text = list()
    for reactant in Heck_rxn.reactants:
        text = text + [str(reactant)] + ["."]

    for reagent in Heck_rxn.reagents:
        text = text + [str(reagent)] + ["."]

    for cat in Heck_rxn.cats:
        text = text + [str(cat)] + ["."]

    for sol in Heck_rxn.solvents:
        text = text + [str(sol)] + ["."]

    text = text[:-1] + [">>"]

    text = text + [str(Heck_rxn.products[0])]
    text = "".join(text)

    drfp = DrfpEncoder.encode(text)[0]

    return drfp

def get_Heck_React_drfp(Heck_rxn):
    text = list()
    for reactant in Heck_rxn.reactants:
        text = text + [str(reactant)] + ["."]

    text = text[:-1] + [">>"]

    text = text + [str(Heck_rxn.products[0])]
    text = "".join(text)

    drfp = DrfpEncoder.encode(text)[0]

    return drfp

def get_Heck_Reagent_drfp(Heck_rxn):

    text = [">"]

    for reagent in Heck_rxn.reagents:
        text = text + [str(reagent)] + ["."]

    for cat in Heck_rxn.cats:
        text = text + [str(cat)] + ["."]

    for sol in Heck_rxn.solvents:
        text = text + [str(sol)] + ["."]

    text = text[:-1] + [">"]
    text = "".join(text)

    drfp = DrfpEncoder.encode(text)[0]

    return drfp

def get_Suzuki_RxnSmi(Suzuki_HTE_df):
    reactant1 = Suzuki_HTE_df.loc["Reactant_1_Name"]
    reactant2 = Suzuki_HTE_df.loc["Reactant_2_Name"]
    cat = Suzuki_HTE_df.loc["Catalyst_1_Short_Hand"]
    ligand = Suzuki_HTE_df.loc["Ligand_Short_Hand"]
    base = Suzuki_HTE_df.loc["Reagent_1_Short_Hand"]
    sol = Suzuki_HTE_df.loc["Solvent_1_Short_Hand"]
    product = "CC1=CC=C2C(C=NN2C2CCCCO2)=C1C1C=C2C=CC=NC2=CC=1"

    if pd.isnull(ligand) and pd.isnull(base):
        text = [str(reactant1)] + ["."] + [str(reactant2)] + [">"] + [str(cat)] + [
            "."] + [str(sol)] + [">"] + [str(product)]
    if pd.isnull(ligand) and pd.isnull(base) == False:
        text = [str(reactant1)] + ["."] + [str(reactant2)] + [">"] + [str(cat)] + [
            "."] + [str(base)] + ["."] + [str(sol)] + [">"] + [str(product)]
    if pd.isnull(ligand) == False and pd.isnull(base):
        text = [str(reactant1)] + ["."] + [str(reactant2)] + [">"] + [str(cat)] + [
            "."] + [str(ligand)] + ["."] + [str(sol)] + [">"] + [str(product)]
    if pd.isnull(ligand) == False and pd.isnull(base) == False:
        text = [str(reactant1)] + ["."] + [str(reactant2)] + [">"] + [str(cat)] + [
            "."] + [str(ligand)] + ["."] + [str(base)] + ["."] + [str(sol)] + [
                   ">"] + [str(product)]

    return "".join(text)

def get_Suzuki_rxnfp(Suzuki_HTE_df):
    reactant1 = Suzuki_HTE_df.loc["Reactant_1_Name"]
    reactant2 = Suzuki_HTE_df.loc["Reactant_2_Name"]
    cat = Suzuki_HTE_df.loc["Catalyst_1_Short_Hand"]
    ligand = Suzuki_HTE_df.loc["Ligand_Short_Hand"]
    base = Suzuki_HTE_df.loc["Reagent_1_Short_Hand"]
    sol = Suzuki_HTE_df.loc["Solvent_1_Short_Hand"]
    product = "CC1=CC=C2C(C=NN2C2CCCCO2)=C1C1C=C2C=CC=NC2=CC=1"

    if pd.isnull(ligand) and pd.isnull(base):
        text = [str(reactant1)] + ["."] + [str(reactant2)] + ["."] + [str(cat)] + [
            "."] + [str(sol)] + [">>"] + [str(product)]
    if pd.isnull(ligand) and pd.isnull(base) == False:
        text = [str(reactant1)] + ["."] + [str(reactant2)] + ["."] + [str(cat)] + [
            "."] + [str(base)] + ["."] + [str(sol)] + [">>"] + [str(product)]
    if pd.isnull(ligand) == False and pd.isnull(base):
        text = [str(reactant1)] + ["."] + [str(reactant2)] + ["."] + [str(cat)] + [
            "."] + [str(ligand)] + ["."] + [str(sol)] + [">>"] + [str(product)]
    if pd.isnull(ligand) == False and pd.isnull(base) == False:
        text = [str(reactant1)] + ["."] + [str(reactant2)] + ["."] + [str(cat)] + [
            "."] + [str(ligand)] + ["."] + [str(base)] + ["."] + [str(sol)] + [
                   ">>"] + [str(product)]
    text = "".join(text)

    model, tokenizer = get_default_model_and_tokenizer()
    rxnfp_generator = RXNBERTFingerprintGenerator(model, tokenizer)
    rxnfp = rxnfp_generator.convert(text)

    return rxnfp

def get_Suzuki_drfp(Suzuki_HTE_df):
    reactant1 = Suzuki_HTE_df.loc["Reactant_1_Name"]
    reactant2 = Suzuki_HTE_df.loc["Reactant_2_Name"]
    cat = Suzuki_HTE_df.loc["Catalyst_1_Short_Hand"]
    ligand = Suzuki_HTE_df.loc["Ligand_Short_Hand"]
    base = Suzuki_HTE_df.loc["Reagent_1_Short_Hand"]
    sol = Suzuki_HTE_df.loc["Solvent_1_Short_Hand"]
    product = "CC1=CC=C2C(C=NN2C2CCCCO2)=C1C1C=C2C=CC=NC2=CC=1"

    if pd.isnull(ligand) and pd.isnull(base):
        text = [str(reactant1)] + ["."] + [str(reactant2)] + ["."] + [str(cat)] + [
            "."] + [str(sol)] + [">>"] + [str(product)]
    if pd.isnull(ligand) and pd.isnull(base) == False:
        text = [str(reactant1)] + ["."] + [str(reactant2)] + ["."] + [str(cat)] + [
            "."] + [str(base)] + ["."] + [str(sol)] + [">>"] + [str(product)]
    if pd.isnull(ligand) == False and pd.isnull(base):
        text = [str(reactant1)] + ["."] + [str(reactant2)] + ["."] + [str(cat)]+ [
            "."] + [str(ligand)] + ["."] + [str(sol)] + [">>"] + [str(product)]
    if pd.isnull(ligand) == False and pd.isnull(base) == False:
        text = [str(reactant1)] + ["."] + [str(reactant2)] + ["."] + [str(cat)] + [
            "."] + [str(ligand)] + ["."] + [str(base)] + ["."] + [str(sol)] + [
                   ">>"] + [str(product)]
    text = "".join(text)

    drfp = DrfpEncoder.encode(text)[0]

    return drfp

