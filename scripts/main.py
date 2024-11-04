from   MLMicelle.utils import traj_reader as trj
import MLMicelle.clustering as clst
import MLMicelle.analysis as mic_anal
import MLMicelle.visualization as plts
import pandas as pd
import numpy as np


traj_file='./data/input/md_cluster_nopbc_trj.gro' ## this file is used to do the radius evaluations
pd_data= trj.gro_reader(traj_file)
clustered_data = clst.clustering(pd_data)
plts.clusters_2D(clustered_data)

