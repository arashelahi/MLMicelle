import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN


def clustering(data,PO_bead):
    data_fit               = data.loc[data['atom_name']==PO_bead,'x':'z']## data for clustering is only based on PO
    clustering             = DBSCAN(eps=0.47, min_samples=3).fit(data_fit) ## DBSCAN is used to do clustering
    data_fit["clusters"]   = clustering.labels_.astype(int) ## add the clusters to data_fit
    data_final             = data.copy()
    data_final['clusters'] = data_fit['clusters'] ## add the clusters to data
    clusters_map           = data_final[data_final['atom_name'] == PO_bead].groupby('residue_number')['clusters'].first()
    data_final['clusters'] = data_final['residue_number'].map(clusters_map)

    return data_final
