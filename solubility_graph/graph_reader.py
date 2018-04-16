import numpy as np
from rdkit import Chem


class predict_solubility():

    def __init__(self, task_parameters, smile_st, solbt): 

        # parameters
        vocab_size = task_parameters['Voc']
        nb_of_opt = task_parameters['nb_clusters_target']
        atm_dict = task_parameters['atm_dict']
        
        
       
        sig_mol = Chem.MolFromSmiles(smile_st)
        sig_mol_h = Chem.AddHs(sig_mol)
        atms_sig_mol = [ atom.GetSymbol()for atom in sig_mol_h.GetAtoms()]
        
        
        #creating u
        u = np.zeros(len(atms_sig_mol))
        for i in range(len(u)):
            u[i] = atm_dict[atms_sig_mol[i]]
        
        
        u=torch.from_numpy(u)
        u=u.long()
        
        
        #creating adj matrix W
        bond_dict ={'SINGLE':1,'DOUBLE':2,'TRIPLE':3,'AROMATIC':4}
   
        W = np.zeros((len(atms_sig_mol),len(atms_sig_mol)))
    
    
        for bond in sig_mol_h.GetBonds():
    
            i = bond.GetBeginAtomIdx()
            j = bond.GetEndAtomIdx()
    
            val = bond_dict[str(bond.GetBondType())]
    
            W[i][j]=val
            W[j][i]=val
        
        


        # create the target
        target= (solbt).astype(float)
        target=torch.from_numpy(target)
        target=target.long()
        

        # mapping matrices
        W_coo=sp.coo_matrix(W)
        nb_edges=W_coo.nnz
        nb_vertices=W.shape[0]
        edge_to_starting_vertex=sp.coo_matrix( ( np.ones(nb_edges) ,(np.arange(nb_edges), W_coo.row) ),
                                               shape=(nb_edges, nb_vertices) )
        edge_to_ending_vertex=sp.coo_matrix( ( np.ones(nb_edges) ,(np.arange(nb_edges), W_coo.col) ),
                                               shape=(nb_edges, nb_vertices) )

        # attribute
        #self.adj_matrix=torch.from_numpy(W).type(default_type)   
        #self.edge_to_starting_vertex=torch.from_numpy(edge_to_starting_vertex.toarray()).type(default_type)
        #self.edge_to_ending_vertex=torch.from_numpy(edge_to_ending_vertex.toarray()).type(default_type)   
        self.adj_matrix=W  
        self.edge_to_starting_vertex=edge_to_starting_vertex
        self.edge_to_ending_vertex=edge_to_ending_vertex  
        self.signal=u
        self.target=target
