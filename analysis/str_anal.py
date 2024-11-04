import numpy as np
import pandas as pd
from sklearn.cluster import KMeans, DBSCAN
import os
from sklearn import metrics
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from statistics import mode


def analysis(data_final,PO_bead): ## generate  size and aggretaion number 
    mice_size_vec=[] ## hydrodynamic radius
    core_size_vec=[]
    agg_numb_vec=[] ## aggregation numbers
    cluster_vec=np.unique(data_final['clusters']) 
    for clust_iter in cluster_vec: ## iterate over each cluster
        clust_data=data_final[data_final['clusters']==clust_iter]
        PO_coord = clust_data[clust_data['atom_name']==PO_bead].loc[:,'x':'z']
        cluster_center = np.mean(PO_coord, axis=0) ## find the center of the cluster
        # PO_distance_list = np.linalg.norm(data_final[data_final['clusters']==clust_iter] - cluster_center, axis=1) ## find the distance of center to all beads in cluster i
        # P1=clus_cont[clus_cont['Bead']==PO_bead]
        # E1=pd.merge(data_final[data_final['clusters']==i],data_final[data_final['Bead']!=PO_bead],how="inner") ## it might need changes
        dist_list_EO=np.linalg.norm(clust_data[clust_data['atom_name']!=PO_bead].loc[:,'x':'z']-cluster_center,axis=1)  ## find the distance of center to all E1 in cluster i
        dist_list_PO=np.linalg.norm(clust_data[clust_data['atom_name']==PO_bead].loc[:,'x':'z']-cluster_center,axis=1)  ## find the distance of center to all E1 in cluster i        rad_mem.append(np.max(dist_list_EO)) ## the maximum distance is hydrodynamic radius
        mice_size_vec.append(np.max(dist_list_EO)) ## the maximum distance is hydrodynamic radius
        core_size_vec.append(np.max(dist_list_PO)) ## the maximum distance is core size
        agg_numb_vec.append(len(np.unique(clust_data['residue_number']))) ## number of atoms in each cluster, devided by P1, gives us aggregation number
    return mice_size_vec, core_size_vec, agg_numb_vec
    