

# Xavier
# Oct 30, 2017


import block
import numpy as np
import pickle


import datetime
t = datetime.datetime.today().replace(second=0, microsecond=0)
t = str(t); t = t.replace(':', '-'); t = t.replace(' ', '_')


p = 0.5
size_subgraph = 20
Voc = 3

set_100_subgraphs = []
nb_subgraphs = 100
for s in range(nb_subgraphs):
	W0 = block.random_graph(size_subgraph,p)
	u0 = np.random.randint(Voc,size=size_subgraph)
	set_100_subgraphs.append([s,W0,u0])


name_file_result = '../data/set_100_subgraphs_p05_size20_Voc3_' + str(t) + '_.txt'
with open(name_file_result, 'wb') as fp: 
	pickle.dump(set_100_subgraphs, fp)
