# from traj_reading import *
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans, DBSCAN
import os
from sklearn import metrics
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from statistics import mode


def molec_extract(data_final):   ## this function extract molecule numbers for each of clusters
    data=[]
    cluster_vec=np.unique(data_final['clusters']) ## clusters assigned to all data of current frames
    for i in cluster_vec:
        clus_cont = data_final[data_final['clusters'] == i]
        resids = np.array([int(x[:-4]) for x in np.unique(clus_cont['PEPP']).tolist()])
        data.append(sorted(resids))
    return data
def fixed_atoms_cluster(frame_range,traj_file, bead_per_mol, mol_num, PO_bead, atom_nums):
    ## this function:
    # 1- generates index files of chains always stay in cluster
    # 2- output aggregation number vector for every frame, considering certain number of clusters at each frame
    data_final= refer_identifier(frame_range, traj_file, bead_per_mol, mol_num, PO_bead, atom_nums)
    reference_clus_moles = molec_extract(data_final) ## a list, that has n_cluster rows, at each row, molecules belong to that cluster exist
    cluster_vec=np.unique(data_final['clusters'])
    n_cluster=len(np.unique(cluster_vec))
    agg_frame_vec=np.zeros((n_cluster, len(frame_range))) ## an array, with row size of n_cluster as reference 
    frame_id=-1 ## sometimes frame_range is not starting from 0
    print('number of frames is %d' % len(frame_range))
    for i in frame_range:
        frame_id=frame_id+1
        print('at frame %d' % i, end= ' ' )
        coordinate, atom_type = reader(traj_file, i)
        data_final = clustering(coordinate, atom_type, bead_per_mol, mol_num, PO_bead, atom_nums)
        curr_clus_moles = molec_extract(data_final)
        for j in range(len(reference_clus_moles)):
            max_com = 0
            for k in range(len(curr_clus_moles)):
                com_item=list(set(reference_clus_moles[j]) & set(curr_clus_moles[k])) ## common items between cluster-jth of reference and cluster k-th of current frame
                if len(com_item)> max_com: ## find the most similar cluster at two frames
                    repl_line=sorted(com_item)
                    max_com=len(com_item)
                    agg_frame_vec[j, frame_id]=len(curr_clus_moles[k]) ## aggregation number of cluster k-th is set for the desired cluster agg number
            if reference_clus_moles[j]!=repl_line: ## if clusters are not the same, it means some molecules are missed
                print(', cluster %d '  %j, end = ' ' )
            reference_clus_moles[j]=repl_line ## replacing the current similar cluster to the reference vector
        print('missed molecules\n')
    ### generate index files
    fixed_ind_file = 'fixed'
    data_fixed=pd.DataFrame([])
    for i in range(len(reference_clus_moles)): ## total number of clusters
        for j in range(len(reference_clus_moles[i])): ## number of molecules at i-th cluster
            ## for all the current molecules, generate a pd.dataframe, that only includes these molecules, then it's used to make the cluster.
            data_cur=data_final[data_final['PEPP'] == '%sPEPP' % str(reference_clus_moles[i][j])].loc[:, ['PEPP', 'Bead', 'clusters']]
            data_cur['clusters']=i
            data_fixed=pd.concat([data_fixed,data_cur])
            # indexing(fixed_ind_file, reference_clus_moles[i], i)
    index_gen(data_fixed, PO_bead, fixed_ind_file)
    return agg_frame_vec

def refer_identifier(frame_range, traj_file, bead_per_mol, mol_num, PO_bead, atom_nums):
    ## this function is used to find most stable number of clusters, 
    print('determining most lasting number of clusters') 
    print('number of steps is %d ' % (int(len(frame_range)/10)))
    ave_n_clus = [] ## a list for number of cluster per frame
    # for i in range(1,int(len(frame_range)/10)):  ## in this range finds the accurate number of clusters is defined as the most happening one
    for i in range(1,int(len(frame_range)/10)):  ## in this range finds the accurate number of clusters is defined as the most happening one
    
        coordinate, atom_type = reader(traj_file, i+1) 
        data_final = clustering(coordinate, atom_type, bead_per_mol, mol_num, PO_bead, atom_nums)
        cluster_vec=np.unique(data_final['clusters'])
        n_cluster=len(np.unique(cluster_vec))
        cluster_membs=np.array([]) ## number of clusters at each frame 
        for j in cluster_vec:
            cluster_membs=np.append(cluster_membs,sum(data_final['clusters']==j)/bead_per_mol) ## finds number of molecules in j-th cluster at each frame
        print('number of clusters is %f and number of unimers is %d' % (np.sum(cluster_membs > 5) ,n_cluster-np.sum(cluster_membs > 5) )) ## print number of big micelles, that is agg_num> 5
        ave_n_clus.append(np.sum(cluster_membs > 5))   ## a list of number of  big clusters, row= frame number
    mode_clus_mem = mode(ave_n_clus)  ## most happening number of clusters, only considering big ones 
    #### finding the reference number
    for i in range (len(ave_n_clus)):  ## the frame which result in the most happening one, specifically the one both sides having same n_cluster
        if abs(ave_n_clus[i]-mode_clus_mem)+abs(ave_n_clus[i+1]-mode_clus_mem)+abs(ave_n_clus[i+2]-mode_clus_mem)==0: ## finds which frame to use
            reference_frame=i+1 ## store the frame number of one we consider as the reference frame 
            break
### making big clusters
    coordinate, atom_type = reader(traj_file, reference_frame+1) ## read the frame
    data_final = clustering(coordinate, atom_type, bead_per_mol, mol_num, PO_bead, atom_nums) ##  finds data at that frame
    data_big_micelles=pd.DataFrame(columns=data_final.columns) ## make a dataframe for big micelles
    cluster_vec_ref=np.unique(data_final['clusters'])
    for j in cluster_vec_ref:
        if sum(data_final['clusters']==j)/bead_per_mol > 5: ## only seach in big micelles
            data_big_micelles=pd.concat([data_big_micelles,data_final[data_final['clusters']==j]]) ## big micelle
    return data_big_micelles
