atm_sym = set([])
max_len = 0

for i in range(len(smiles)):
    print(i)
    tst_mol = Chem.MolFromSmiles(smiles[i])
    tst_mol_h = Chem.AddHs(tst_mol)
    tst_atoms = [ atom.GetSymbol()for atom in tst_mol_h.GetAtoms()]
    #if len(tst_atoms) > max_len :
        #max_len = len(tst_atoms)
    
    #print(tst_atoms)
    
    try1 = set(tst_atoms)
    print(try1)
    atm_sym = atm_sym.union(try1)
    print(atm_sym)
    
print("unique atoms = ", atm_sym)
#print(max_len)
