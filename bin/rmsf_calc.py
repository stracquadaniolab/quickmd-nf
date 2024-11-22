#!/usr/bin/env python3
"""rmsf_calc

Run rmsfs via MD analysis
Usage:
rmsf_calc.py [--inpdb=<pdb>]

Options:
--inpdb=<pdb>      Load input pdb trajectory
"""

import logging
from docopt import docopt
import MDAnalysis as mda
from MDAnalysis import transformations
from MDAnalysis.analysis import rms, align
import sys
import pathlib
import pandas as pd
import numpy as np
from sys import stdout

def rmsf_analysis(pdb_traj: str, rmsfcsv: str):
    u = mda.Universe(pdb_traj)
    ag = u.atoms
    new_dimensions = [114, 114, 114, 90, 90, 90]
    set_dim = transformations.boxdimensions.set_dimensions(new_dimensions)
    transform = transformations.unwrap(ag)
    center = transformations.center_in_box(ag.select_atoms('protein'), wrap=True)
    u.trajectory.add_transformations(set_dim, transform, center)
    #u.trajectory.dt = 0.001
    average = align.AverageStructure(u, u, select='protein and name CA', ref_frame=0).run()
    ref = average.results.universe
    aligner = align.AlignTraj(u, ref, select='protein and name CA', in_memory=True).run()
    c_alphas = u.select_atoms('protein and name CA')
    R = rms.RMSF(c_alphas).run()
    residue_ids = [residue.resid for residue in c_alphas.residues]
    residue_names = [residue.resname for residue in c_alphas.residues]
    rmsf_values = R.rmsf
    df = pd.DataFrame({'Res_ID': residue_ids, 'Res_Name': residue_names, 'RMSF': rmsf_values })
    print(df)
    df.to_csv(rmsfcsv, index=False)

def rmsf_analysis_by_chain(pdb_traj: str, rmsfcsv: str):
    u = mda.Universe(pdb_traj)
    ag = u.atoms
    new_dimensions = [117, 117, 117, 90, 90, 90]
    set_dim = transformations.boxdimensions.set_dimensions(new_dimensions)
    transform = transformations.unwrap(ag)
    center = transformations.center_in_box(ag.select_atoms('protein'), wrap=True)
    u.trajectory.add_transformations(set_dim, transform, center)
    chain_ids = np.unique(u.select_atoms('protein').atoms.chainIDs)
    for chain_id in chain_ids:
        chain = u.select_atoms(f'protein and chainID {chain_id}')
        align.AlignTraj(u, chain, select=f'protein and chainID {chain_id} and name CA', in_memory=True).run()
        c_alphas = chain.select_atoms('protein and name CA')
        R = rms.RMSF(c_alphas).run()
        residue_ids = [residue.resid for residue in c_alphas.residues]
        residue_names = [residue.resname for residue in c_alphas.residues]
        rmsf_values = R.rmsf
        df = pd.DataFrame({'Res_ID': residue_ids, 'Res_Name': residue_names, 'RMSF': rmsf_values })
        print(df)
        rmsfout = str(rmsfcsv + "_rmsf1_ch" + chain_id + ".csv")
        df.to_csv(rmsfout, index=False)

def main():
    arguments = docopt(__doc__, version='rmsf_calc.py')
    inpdb = str(arguments['--inpdb'])
    stem = inpdb.replace(".pdb","")
    #rmsfout = str(stem + "_rmsf1.csv")
    #rmsf_analysis(inpdb, rmsfout)
    rmsf_analysis_by_chain(inpdb, stem)

if __name__ == '__main__':
    arguments = docopt(__doc__, version='rmsf_calc.py')
    logging.getLogger().setLevel(logging.INFO)  
    main()