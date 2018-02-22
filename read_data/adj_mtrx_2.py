from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
import numpy as np
import pandas as pd
   
def mol_adj_2(smiles_txt):
    
    mol = Chem.MolFromSmiles(smiles_txt)
    mol = Chem.AddHs(mol)
   
    atoms = [ atom.GetSymbol()for atom in mol.GetAtoms()]
   
    bond_dict ={'SINGLE':1,'DOUBLE':2,'TRIPLE':3,'AROMATIC':4}
   
    adj_mt_2 = np.zeros((len(atoms),len(atoms)))
    
    
    for bond in mol.GetBonds():
    
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
    
        val = bond_dict[str(bond.GetBondType())]
    
        adj_mt_2[i][j]=val
        adj_mt_2[j][i]=val
    
    return atoms, adj_mt_2
    
   
atms, mtrx1 = mol_adj_2('C=C')
    
print(atms)
print(mtrx1)
    
qm95 = pd.read_csv('qm9_5k.csv',header=None)
   
for i in range(10):
    atms, mtrx1 = mol_adj_1(str(qm95[0][i]))
    
    print(atms)
    print(mtrx1)
