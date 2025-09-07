import numpy as np
import pandas as pd
import MDAnalysis as mda
"""
This module provides functions to analyze clustered micelle data,
including calculating micelle size, core size, aggregation number,
and anisotropy based on clustering results.
"""

def analysis(
    u: mda.Universe,
    data_final: pd.DataFrame,
    PO_bead: str = 'PO40',
    ignore_noise: bool = True,
) -> tuple[list[float], list[float], list[int]]:
    """Compute micelle size (hydrodynamic radius), core size, and aggregation number per cluster.

    Parameters
    - data_final: DataFrame containing at least ['clusters', 'atom_name', 'residue_number', 'x', 'y', 'z'].
    - PO_bead: Name of the PO bead used to define the cluster center.

    Returns
    - mice_size_vec: Max distance of non-PO atoms to PO center for each cluster.
    - core_size_vec: Max distance of PO atoms to PO center for each cluster.
    - agg_numb_vec: Number of unique residues in each cluster.
    """

    required_cols = {'clusters', 'atom_name', 'residue_number', 'x', 'y', 'z'}
    missing = required_cols - set(data_final.columns)
    if missing:
        raise ValueError(f"analysis(): missing required columns: {sorted(missing)}")

    # Prepare cluster list, handle NaNs and optional noise filtering
    clusters = pd.Series(data_final['clusters']).dropna().unique()
    clusters = sorted(int(c) for c in clusters)


    mice_size_vec: list[float] = []  # hydrodynamic radius
    core_size_vec: list[float] = []
    agg_numb_vec: list[int] = []
    anisotropy_vec: list[float] = []

    for clabel in clusters:
        clust_data = data_final[data_final['clusters'] == clabel]

        po = clust_data[clust_data['atom_name'] == PO_bead].loc[:, ['x', 'y', 'z']]
        if len(po) == 0:
            # Fallback: center from all atoms in the cluster
            center = clust_data.loc[:, ['x', 'y', 'z']].mean(axis=0).to_numpy()
        else:
            center = po.mean(axis=0).to_numpy()

        # Distances
        eo_coords = clust_data[clust_data['atom_name'] != PO_bead].loc[:, ['x', 'y', 'z']]
        po_coords = po

        mice_size = float(np.max(np.linalg.norm(eo_coords - center, axis=1))) if len(eo_coords) else np.nan
        core_size = float(np.max(np.linalg.norm(po_coords - center, axis=1))) if len(po_coords) else np.nan

        mice_size_vec.append(mice_size)
        core_size_vec.append(core_size)
        agg_numb_vec.append(int(clust_data['residue_number'].nunique()))
        anisotropy_vec.append(anisotropy(u, clust_data))

    return mice_size_vec, core_size_vec, agg_numb_vec, anisotropy_vec
def anisotropy(u, clust_data):
    """Compute anisotropy from squared radii of gyration along principal axes.

    Parameters
    - evals: Iterable of three floats, squared radii of gyration along
      the major, middle, and minor principal axes.

    Returns
    - Anisotropy value (float).
    """
    # evals = [1,1,1]
    evals = rg2_principal_from_indices(u, clust_data)
    I1 = (evals[0] + evals[1] + evals[2])**2
    I2 = evals[0]*evals[1] + evals[1]*evals[2] + evals[2]*evals[0]
    if I1 == 0:
        return 0.0
    k2 = 1 - 3*I2/(I1)
    Rg2_major, Rg2_middle, Rg2_minor = evals
    Rg2_total = Rg2_major + Rg2_middle + Rg2_minor
    if Rg2_total == 0:
        return 0.0
    return 1.5 * ((Rg2_major**2 + Rg2_middle**2 + Rg2_minor**2) / (Rg2_total**2)) - 0.5
def rg2_principal_from_indices(universe, clust_data):
    """
    Compute squared radii of gyration along the principal axes
    for a given list of atom indices.
    """
    atom_indices = clust_data['atom_number'].values.tolist()
    sel_str = "index " + " ".join(map(str, atom_indices))
    ag = universe.select_atoms(sel_str)
    results = []

    for ts in universe.trajectory:
        r = ag.positions - np.average(ag.positions, axis=0, weights=ag.masses)
        M = ag.masses.sum()
        S = (ag.masses[:, None] * r).T @ r / M
        evals, evecs = np.linalg.eigh(S)
        evals = evals[::-1]  # major, middle, minor

        # results.append({
        #     "frame": ts.frame,
        #     "time_ps": ts.time,
        #     "Rg2_major": float(evals[0]),
        #     "Rg2_middle": float(evals[1]),
        #     "Rg2_minor": float(evals[2]),
        #     "Rg2_total": float(evals.sum())
        # })
    return evals