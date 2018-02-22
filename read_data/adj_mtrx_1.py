from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
import numpy as np
import pandas as pd

def mol_adj_1(smiles_txt):

    mol = Chem.MolFromSmiles(smiles_txt)
    mol = Chem.AddHs(mol)
    
    atoms = [atom.GetSymbol()for atom in mol.GetAtoms()]
    
    adj_mt_1 = np.zeros((len(atoms),len(atoms)))
    
    
    for bond in mol.GetBonds():
        
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
    
        adj_mt_1[i][j]=1
        adj_mt_1[j][i]=1
    
    return atoms, adj_mt_1    
    
atms, mtrx1 = mol_adj_1('C=C')
    
print(atms)
print(mtrx1)
    
qm95 = pd.read_csv('qm9_5k.csv',header=None)
   
for i in range(10):
    atms, mtrx1 = mol_adj_1(str(qm95[0][i]))
    
    print(atms)
    print(mtrx1)
