from   MLMicelle.utils import traj_reader as trj
import MLMicelle.clustering as clst
import MLMicelle.analysis as mic_anal
import MLMicelle.visualization as plts
import pandas as pd
import numpy as np
import MDAnalysis as mda

# Path to the trajectory and topology files
traj_dir = './data/input/' 

## load the trajectory file using MDAnalysis
u = mda.Universe(traj_dir + '/md.tpr', traj_dir + '/md_clus.gro')
## Process the Trajectory Information to a pd.DataFrame
clustered_data = clst.clustering(u, bead_type='PO55') 
# Generating three index files, one for PO beads, one for all beads, and one for molecule numbers
clst.index_gen(clustered_data, PO_bead='PO55', index_names='./data/output/clustered') 
# Analyzing the clustered data to get micelle size, core size, aggregation number, and anisotropy
mice_size, core_size, agg_numb, anisotropy =  mic_anal.analysis(u, clustered_data) 
# Printing the Results
print('anisotropy:', anisotropy)
