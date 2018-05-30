

# Xavier
# Oct 30, 2017


import graph_generator as g
import numpy as np
import pickle


import datetime
t = datetime.datetime.today().replace(second=0, microsecond=0)
t = str(t); t = t.replace(':', '-'); t = t.replace(' ', '_')



task_parameters = {}
task_parameters['flag_task'] = 'clustering'
task_parameters['nb_communities'] = 10
task_parameters['nb_clusters_target'] = task_parameters['nb_communities']
task_parameters['Voc'] = task_parameters['nb_communities'] + 1
task_parameters['size_min'] = 5
task_parameters['size_max'] = 25
task_parameters['p'] = 0.5
task_parameters['q'] = 0.1



set_100_clustering_maps = []
nb_clustering_maps = 100
for s in range(nb_clustering_maps):
	train_x = g.graph_semi_super_clu(task_parameters)
	set_100_clustering_maps.append([s,train_x])


name_file_result = '../data/set_100_clustering_maps_p05_q01_size5_25_' + str(t) + '_.txt'
with open(name_file_result, 'wb') as fp: 
	pickle.dump(set_100_clustering_maps, fp)




