"""
This Module contains some tools to handle the molecules:
    1.Convert IUPAC into SMILES
    2.SMILES Checker & Normalization
    3.utils for MMHRP model
    4.Determination of solvent polarity
"""

import re
from time import sleep
from lxml import etree
from bs4 import BeautifulSoup
from urllib.request import urlopen
from urllib.parse import quote
import torch
from torch_geometric.data import Data
import numpy as np
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import matplotlib.pyplot as plt
import torch_geometric.transforms as T

# 1. IUPAC2SMILES
def iupac_to_smi_cactus(name):
    if name == "/":
        return None
    try:
        url = 'http://cactus.nci.nih.gov/chemical/structure/' + quote(name) + '/smiles'
        smi = urlopen(url).read().decode('utf8')
        return smi
    except:
        return None

def iupac_to_smi_pubchem(name, bro):
    try:
        # open pubchem
        bro.get("https://pubchem.ncbi.nlm.nih.gov/#query=%s" % name)
        sleep(8)
        # Get current webpage information
        soup = BeautifulSoup(bro.page_source, "html.parser")
        feature_res = soup.find_all("div", id="featured-results")
        html_parser = etree.HTML(str(feature_res[0]))
        # smi for obtaining the best matching result
        label = html_parser.xpath('./descendant::*/text()')
        num = len(label)
        # If there are no search results
        if num == 0:
            return None
        for i in range(num):
            if label[i] == "Isomeric SMILES: ":
                smi = label[i + 1]
                break
            # If there is no smi
            if i == num - 1:
                return None
        return smi
    except:
        return None

class molecule():
    """
    Class of molecules
    """
    def __init__(self):
        self.name = str()
        self.smi = str()
    def record(self, name, smi):
        self.name = name
        self.smi = smi

class Mol_Manager():
    """
    Molecular name management library
    """
    def __init__(self, bro):
        self.no_smi = list() # Unable to find smi
        self.yes_smi = list() # Able to find smi
        self.bro = bro

    def get_smi(self, name):
        # in no_smi
        pd_no = False
        for mol in self.no_smi:
            if name == mol.name:
                pd_no = True
                return mol.smi

        # in yes_smi
        pd_yes = False
        for mol in self.yes_smi:
            if name == mol.name:
                pd_yes = True
                return mol.smi

        # Neither of the two
        if pd_no == False and pd_yes == False:
            # get smi
            smi = iupac_to_smi_cactus(name)
            if smi == None: # If cactus does not work, use PubChem
                smi = iupac_to_smi_pubchem(name, bro=self.bro)

            # put it into list
            mol = molecule()
            mol.record(name, smi)
            if smi == None:
                self.no_smi.append(mol)
            else:
                self.yes_smi.append(mol)
            return smi

# 2.SMILES Checker & Normalization
from rdkit import Chem
from rdkit.Chem import MolStandardize

def smi_checker(rxn):
    pd = True
    try:
        for reactant in rxn.reactants:
            if Chem.MolFromSmiles(reactant) == None:
                pd = False
                break
        for product in rxn.products:
            if Chem.MolFromSmiles(product) == None:
                pd = False
                break
        for reagent in rxn.reagents:
            if Chem.MolFromSmiles(reagent) == None:
                pd = False
                break
        for cat in rxn.cats:
            if Chem.MolFromSmiles(cat) == None:
                pd = False
                break
        for sol in rxn.solvents:
            if Chem.MolFromSmiles(sol) == None:
                pd = False
                break
    except:
        pd = False

    return pd

def Mol_Clean(smi, Uncharge=True, ChooseLargeFrag=True):
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            # normalization
            mol = MolStandardize.normalize.Normalizer().normalize(mol)
            # choose large fragment
            if ChooseLargeFrag:
                mol = MolStandardize.fragment.LargestFragmentChooser().choose(mol)
            # uncharged
            if Uncharge:
                mol = MolStandardize.charge.Uncharger().uncharge(mol)
            # metal disconnection
            mol = MolStandardize.rdMolStandardize.MetalDisconnector().Disconnect(mol)
            mol = MolStandardize.rdMolStandardize.Cleanup(mol)

            smi = Chem.MolToSmiles(mol, isomericSmiles=False, canonical=True)
            return smi
    except:
        return False

#  3.utils for MMHRP model
def AtomCharityEmbed(atom):
    C = str(atom.GetChiralTag())
    OH = [0 for i in range(3)] # [Whether it is Chiral or not, R, S]
    # No Chiral
    if C == "CHI_UNSPECIFIED":
        return OH
    # R
    if C == "CHI_TETRAHEDRAL_CW":
        OH[0] = 1
        OH[1] = 1
        return OH
    # S
    if C == "CHI_TETRAHEDRAL_CCW":
        OH[0] = 1
        OH[2] = 1
        return OH

    return OH

def AtomHybridizationEmbed(atom):
    HBZ = str(atom.GetHybridization())
    OH = [0 for i in range(6)] # [sp, sp2, sp3, sp3d, sp3d2, other]
    DICT = ["SP", "SP2", "SP3", "SP3D", "SP3D2"]
    if HBZ == "SP":
        OH[0] = 1
        return OH
    if HBZ == "SP2":
        OH[1] = 1
        return OH
    if HBZ == "SP3":
        OH[2] = 1
        return OH
    if HBZ == "SP3D":
        OH[3] = 1
        return OH
    if HBZ == "SP3D2":
        OH[4] = 1
        return OH
    if HBZ not in DICT:
        OH[5] = 1
        return OH

def smi_to_graph(smi, ChargeContribution=True, NoFeature=False):
    """
    :param smi:
    :return:
        :num: number of atoms
        :feature: feature of nodes, include AtomicNumber & FormalCharge
        :edge_index: edge information
    """
    # Get num
    mol = Chem.MolFromSmiles(smi)
    num = mol.GetNumAtoms()
    # Get label
    feature = list()
    atoms = mol.GetAtoms()
    # Get Charge contribute
    AllChem.ComputeGasteigerCharges(mol)

    for atom in atoms:
        if NoFeature == False:
            atom_feature = [
                int(atom.GetAtomicNum()), # Atomic Number,
                int(atom.GetFormalCharge()), # Formal Charge,
                int(atom.GetTotalNumHs()), # Number of connected H,
                int(atom.GetExplicitValence()), # Explicit Valence,
                int(atom.GetDegree()), # Degree of an atom
                int(atom.GetIsAromatic()), # is aromatic or not
                int(atom.IsInRing()), # is in ring or not
            ]
            # atom_feature += AtomCharityEmbed(atom) # Atom Charity
            # atom_feature += AtomHybridizationEmbed(atom) # Atom Hybridization
            if ChargeContribution:
              atom_feature.append(float(atom.GetProp("_GasteigerCharge")))  # Gasteiger Charge Contribution
            
            feature.append(atom_feature)
        else:
            atom_feature = [
                int(0),
            ]
            feature.append(atom_feature)

    # Get edge_index
    us = list()
    vs = list()
    bonds = mol.GetBonds()
    for bond in bonds:
        u = bond.GetBeginAtom().GetIdx()
        v = bond.GetEndAtom().GetIdx()
        us.append(u)
        vs.append(v)
        us.append(v)
        vs.append(u)
    edge_index = [us, vs]

    return num, feature, edge_index

def smis_to_graph(smis, ChargeContribution=True, NoFeature=False):
    temp_num = 0
    g_feature = list()
    g_us = list()
    g_vs = list()
    for smi in smis:
        num, feature, edge_index = smi_to_graph(smi, ChargeContribution=ChargeContribution, NoFeature=NoFeature)
        for atom in feature:
            g_feature.append(atom)
        for u in edge_index[0]:
            g_us.append(u + temp_num)
        for v in edge_index[1]:
            g_vs.append(v + temp_num)
        temp_num += num
    g_edge_index = [g_us, g_vs]

    x = torch.tensor(g_feature, dtype=torch.float32)
    edge_index = torch.tensor(g_edge_index, dtype=torch.long)
    data = Data(x=x, edge_index=edge_index)

    data = T.NormalizeFeatures()(data)
    data = T.AddSelfLoops()(data) # Self Loop

    return data

# SMILES Tokenizer
def smi_tokenizer(smi):
    """
    Tokenize a SMILES molecule
    """
    pattern = r"(\%\([0-9]{3}\)|\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\||\(|\)|\.|=|#|-|\+|\\|\/|:|~|@|\?|>>?|\*|\$|\%[0-9]{2}|[0-9])"
    regex = re.compile(pattern)
    tokens = [token for token in regex.findall(smi)]
    return tokens

# Format conversion: vocab_list & vocab_dict
def vocab_dict_to_txt(dict, rxn_name):
    with open("../models/%s_vocab.txt" % rxn_name, mode="w", encoding="utf-8") as txt:
        for key in dict.keys():
            txt.write("%s\t%s\t\t" % (key, dict[key]))
    txt.close()

def vocab_txt_to_dict(file):
    vocab_dict = dict()
    txt = open(file, mode="r", encoding="utf-8")
    vocab = txt.read().replace("\n", "").split("\t\t")
    vocab = vocab[:-1] # remove the last one, which is ""

    for pair in vocab:
        pair = pair.split("\t")
        vec = pair[1][1:-1].replace("  ", " ").split(" ") # delete "[", "]" and split by " "
        while "" in vec:
            vec.remove("") # remove null element
        vocab_dict[pair[0]] = np.array(vec).astype(float)

    return vocab_dict

# 4.Determination of solvent polarity
def PolarType(smi, index_addr):
    # Read the index
    sol_dict = {}
    with open(index_addr, mode="r") as f:
        for line in f.readlines():
            sol = line.split(":")[0].lower()
            type = line.split(":")[-1][:-1]
            sol_dict[sol] = type
    if smi.lower() not in sol_dict.keys():
        print("%s is not in the index." % smi)
        return None
    else:
        return sol_dict[smi.lower()]
