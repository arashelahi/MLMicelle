# # import the file
import numpy as np
# import pandas as pd
import traj_reader  as trj
import sys
import os
import clustering as clst
import analysis as mic_anal
import visualization as plts

# ##### GIVEN VALUES
PO_bead='PO40'
rad_agg_traj_file='./data/input/md_cluster_nopbc_trj.gro' ## this file is used to do the radius evaluations
num_frame_size=trj.frame_count(rad_agg_traj_file) ## determining number of steps of size evalation
n_clus_vec=[]

### find size changes with time
# for frame_iter in range (num_frame_size):
frame_iter=0
pd_data= trj.gro_reader(rad_agg_traj_file, frame_iter+1)
data_final = clst.clustering(pd_data,  PO_bead)
# plts.clusters_2D(data_final,PO_bead)
cluster_num = 1
mice_size_vec, core_size_vec, agg_numb_vec = mic_anal.analysis(data_final,PO_bead)
plts.size_2D(data_final[data_final['clusters']==cluster_num],PO_bead)
# n_clus_vec.append(data_final['clusters'].nunique())
# # print("frame number: %d      cluster number=%d" % (frame_iter+1 ,data_final['clusters'].nunique()))
# # print(np.average(n_clus_vec))
# mice_size_vec, core_size_vec, agg_numb_vec = mic_anal.analysis(data_final,PO_bead)
# print(mice_size_vec, core_size_vec, agg_numb_vec)

