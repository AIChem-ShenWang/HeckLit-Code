"""
This Module contains some tools to handle the molecules:
    1.Convert IUPAC into SMILES
    2.SMILES Checker & Normalization
"""

# 1.Convert IUPAC into SMILES
import re
from time import sleep
from lxml import etree
from bs4 import BeautifulSoup
from urllib.request import urlopen
from urllib.parse import quote

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