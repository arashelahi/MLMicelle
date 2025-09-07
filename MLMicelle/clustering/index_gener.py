import numpy as np
import pandas as pd
from sklearn.cluster import KMeans, DBSCAN
import os
from sklearn import metrics
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from statistics import mode


def indexing(index_file,atom_array,select_name): ## generates index file, having molecules at each cluter

    # os.system('rm  %s ' % (index_file))
    index_read=open(index_file,'a') 
    index_read.write('[cluster%d]\n' % select_name)
    step = 15 ## because index file compatible with gromacs in every line there're 15 molecule numbers
    n = 0
    while n < len(atom_array):
        line_ind = atom_array[n:n + step]
        index_read.write("".join('%8s' % str(int(x)) for x in line_ind))
        index_read.write('\n')
        n = n + step
    index_read.close()
    
def index_gen(data_final,PO_bead,index_names): ## preparing molecule, PO and all bead numbers for indexing function
    clus_P1_file=index_names+'_PO.ndx'  ## only p1 beads
    clus_all_file=index_names+'_all.ndx'  ## both E1, P1 beads
    clus_moldcule_file=index_names+'_moles.ndx' ## molecule numbers
    cluster_vec=np.unique(data_final['clusters']) ## list of names of each cluster 
    os.system('rm  %s %s %s ' % (clus_P1_file,clus_all_file,clus_moldcule_file))
    for i in cluster_vec: ## iterate over all existing clusters
        clus_cont=data_final[data_final['clusters']==i] ## extract data of i-th cluters
        P1_numb=(clus_cont[clus_cont['atom_name']==PO_bead].index+1).tolist() ## PO beads of i-th cluster
        all_numb=(clus_cont.index+1).tolist() ## all beads of i-th cluster
        molecule_numb=sorted([x for x in np.unique(clus_cont['residue_number']).tolist()]) ## molecules of i-th cluster
        indexing(clus_P1_file, P1_numb, i)
        indexing(clus_all_file, all_numb, i)
        indexing(clus_moldcule_file, molecule_numb, i)
