import numpy as np
import pandas as pd
from typing import Literal
from sklearn.cluster import DBSCAN
import MDAnalysis as mda


def mda_to_pdb(u):
    """Return a pandas DataFrame like traj_reader.gro_reader for a given MDAnalysis Universe `u`.

    Tries to read from the underlying file if available; otherwise builds
    the DataFrame from the Universe atom attributes.
    """
    atoms = []
    for a in u.atoms:
        atoms.append({
            'residue_number': int(a.resid),
            'residue_name': a.resname,
            'atom_name': a.name,
            'atom_number': int(getattr(a, 'id', a.index + 1)),
            'x': float(a.position[0]/10),  # convert from Angstrom to nm
            'y': float(a.position[1]/10),
            'z': float(a.position[2]/10),
        })
    return pd.DataFrame(atoms)


def clustering(
    mda_data: mda.Universe,
    bead_type: str = 'PO40',
    eps: float = 0.47,
    min_samples: int = 3,
) -> pd.DataFrame:

    """Cluster bead_type beads with DBSCAN and propagate labels to residues.

    Parameters
    - data: Input DataFrame with columns at least
      ['atom_name', 'residue_number', 'x', 'y', 'z'].
    - bead_type: Bead name used for clustering (default 'PO40').
    - eps: DBSCAN neighborhood radius (in coordinate units of data).
    - min_samples: DBSCAN `min_samples` parameter.
    Returns
    - Copy of `data` with an added 'clusters' column containing integer
      labels per residue (DBSCAN labels, -1 for noise). Residues without
      a bead_type bead will have NaN.
    """
    data = mda_to_pdb(mda_data) ## convert MDAnalysis object to pandas dataframe
    # Select bead_type beads and coordinates used for clustering
    bead_mask = data['atom_name'] == bead_type
    if not bead_mask.any():
        # No bead_type beads found; return data with NaN clusters to match shape
        out = data.copy()
        out['clusters'] = np.nan
        return out
    
    bead_df = data.loc[bead_mask, ['residue_number', 'x', 'y', 'z']].copy()

    # Run DBSCAN on bead coordinates
    coords = bead_df[['x', 'y', 'z']]
    db = DBSCAN(eps=eps, min_samples=min_samples).fit(coords)
    bead_df['clusters'] = db.labels_.astype(int)

    # Aggregate labels per residue
    clusters_map = (
        bead_df.groupby('residue_number')['clusters']
        .first()
        )
    # Propagate residue labels to all atoms in that residue
    data_final = data.copy()
    data_final['clusters'] = data_final['residue_number'].map(clusters_map)

    return data_final
