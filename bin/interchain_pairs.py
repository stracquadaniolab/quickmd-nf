#!/usr/bin/env python3
"""interchain_pairs

Find closest atom pairs at interchain interface to evidence protein secondary structure
Usage:
interchain_pairs.py [--inpdb=<pdb>]

Options:
--inpdb=<pdb>      Input PDB file of protein as obtained from previous process
"""
import logging
from docopt import docopt
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis import rms
from MDAnalysis import transformations
import pandas as pd
from biopandas.pdb import PandasPdb
from typing import Optional, Tuple, List
import numpy as np

def average_trajectory(pdb: str, pdbout: str):
    # Load the structure and trajectory
    u = mda.Universe(pdb)
    # Create a new Universe with the same topology but without coordinates
    #avg_universe = mda.Merge(u.atoms)
    ag = u.atoms
    new_dimensions = [117.0, 117.0, 117.0, 90, 90, 90]
    set_dim = transformations.boxdimensions.set_dimensions(new_dimensions)
    transform = transformations.unwrap(ag)
    center = transformations.center_in_box(ag.select_atoms('protein'), wrap=True)
    u.trajectory.add_transformations(set_dim, transform, center)
    protein = u.select_atoms("protein")
    avg_universe = mda.Merge(protein)
    avg_universe.add_TopologyAttr('tempfactors')
    #avg_coordinates = avg_universe.atoms.positions
    avg_coordinates = np.zeros_like(avg_universe.atoms.positions)
    # Loop over frames, summing up coordinates
    for ts in u.trajectory:
        avg_coordinates += protein.positions
        #avg_coordinates += u.atoms.positions
    # Compute the average
    avg_coordinates /= len(u.trajectory)
    print(len(u.trajectory))
    # Assign average coordinates back to avg_universe
    avg_universe.atoms.positions = avg_coordinates
    # Write the average structure to a PDB file
    avg_universe.atoms.write(pdbout)


def get_contact_atoms(pdbout: str, threshold: float):
    #read PDB data in pandas dataframe
    pdb_data = PandasPdb().read_pdb(pdbout)
    #pdb_df = pd.concat([pdb_data.df['ATOM'], pdb_data.df['HETATM']])
    pdb_df = pd.concat([pdb_data.df['ATOM']])
    pdb_df = pdb_df.dropna(subset=['residue_number'])
    #Strings of coordinates, chains and CA to refine dataframe  
    coord_names = ['x_coord', 'y_coord', 'z_coord']
    chain1 = "A"
    chain2 = "B"
    calpha = "CA"
    #Separate chains into separate dataframes
    df1 = pdb_df[(pdb_df['chain_id'] == chain1) & (pdb_df['atom_name'] == calpha)]
    df2 = pdb_df[(pdb_df['chain_id'] == chain2) & (pdb_df['atom_name'] == calpha)]
    #Extract coordinates to numpy
    coords1 = df1[coord_names].to_numpy()
    coords2 = df2[coord_names].to_numpy()
    #Calculate interchain distances
    dist_matrix = np.sqrt(((coords1[:, None] - coords2) ** 2).sum(axis=2))
    # Create a new dataframe containing pairs of atoms whose distance is below the threshold
    pairs = np.argwhere(dist_matrix < threshold)
    print(f"Pairs: {pairs.shape}")
    print(pairs)
    #Identify chain and redidue of atom pairs within distance threshold
    atoms1, atoms2 = df1.iloc[pairs[:, 0]], df2.iloc[pairs[:, 1]]
    distances = dist_matrix[pairs[:, 0], pairs[:, 1]]
    print(f"Length of atoms1: {len(atoms1)}")
    print(f"Length of atoms2: {len(atoms2)}")
    print(f"Length of distances: {len(distances)}")
    print(distances)
    atoms1_id = atoms1['chain_id'].map(str) + ":" + atoms1['residue_name'].map(str) + ":" + atoms1['residue_number'].map(str)
    atoms2_id = atoms2['chain_id'].map(str) + ":" + atoms2['residue_name'].map(str) + ":" + atoms2['residue_number'].map(str)
    node_pairs = np.vstack((atoms1_id.values, atoms2_id.values, distances)).T
    #node_pairs_df = pd.DataFrame({ 'Atom1_ID': atoms1['chain_id'].map(str) + ":" + atoms1['residue_name'].map(str) + ":" + atoms1 ['residue_number'].map(str), 'Atom2_ID': atoms2['chain_id'].map(str) + ":" + atoms2['residue_name'].map(str) + ":" + atoms2['residue_number'].map(str), 'Distance': distances})
    #node_pairs_df = pd.DataFrame({
    #'Atom1_ID': atoms1['chain_id'].astype(str) + ":" + atoms1['residue_name'] + ":" + atoms1['residue_number'].astype(str),
    #'Atom2_ID': atoms2['chain_id'].astype(str) + ":" + atoms2['residue_name'] + ":" + atoms2['residue_number'].astype(str),
    #'Distance': distances})

    result = pd.concat([df1.iloc[np.unique(pairs[:, 0])], df2.iloc[np.unique(pairs[:, 1])]])
    return node_pairs, result
    #return node_pairs_df, result

def main():
    arguments = docopt(__doc__, version='interchain_pairs.py')
    pdb = arguments['--inpdb']
    pdbstem = pdb.replace(".pdb","")
    pdbout = str(pdbstem + "_average.pdb") 
    csvout = str(pdbstem + "_interchain_pairs.csv")
    average_trajectory(pdb, pdbout)
    threshold = 15.0
    out = get_contact_atoms(pdbout, threshold)
    node_pairs = out[0]
    node_pairs_df = pd.DataFrame(node_pairs, columns=['Atom1', 'Atom2', 'Distance'])
    node_pairs_df.to_csv(csvout)

if __name__ == '__main__':
    arguments = docopt(__doc__, version='interchain_pairs.py')
    logging.getLogger().setLevel(logging.INFO)  
    main()
