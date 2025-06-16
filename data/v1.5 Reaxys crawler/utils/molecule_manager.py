from selenium import webdriver
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
            smi = label[i+1]
            break
        # If there is no smi
        if i == num - 1:
            return None
    return smi

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
    Molecular name and mol management library
    """

    def __init__(self, bro):
        self.no_smi = list() # Unable to find smi
        self.yes_smi = list() # Able to find smi
        self.bro = bro

    def get_smi(self, name, cactus=False):
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
            if cactus == True:
                smi = iupac_to_smi_cactus(name)
                if smi == None: # If cactus does not work, use PubChem
                    smi = iupac_to_smi_pubchem(name, bro=self.bro)
            else:
                smi = iupac_to_smi_pubchem(name, bro=self.bro)
            # put it into list
            mol = molecule()
            mol.record(name, smi)
            if smi == None:
                self.no_smi.append(mol)
            else:
                self.yes_smi.append(mol)
            return smi

# test
if __name__ == "__main__":
    bro = webdriver.Chrome("./chromedriver.exe")
    mol_manager = Mol_Manager(bro=bro)
    ans = list()
    ans.append(mol_manager.get_smi("water", cactus=True))
    ans.append(mol_manager.get_smi("water"))
    ans.append(mol_manager.get_smi("khvivi", cactus=True))
    ans.append(mol_manager.get_smi("Aspirin"))
    ans.append(mol_manager.get_smi("/", cactus=True))
    bro.quit()
    print(ans)
    print(mol_manager.no_smi)
    print(mol_manager.yes_smi)