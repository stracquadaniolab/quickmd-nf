#!/usr/bin/env python3
"""openmm-MD

Minimise the potential energy of the wildtype protein structure and mutated variants, according to AMBER14 potentials 
Usage:
openmmMD.py [--incsv=<incsv>] [--from-col=<col>] [--in-pdb=<pdb>] [--chain=<chain>] [--pH=<pH>] [--no_restraints] 

Options:
--incsv=<incsv>    Read data from csv file
--from-col=<col>   Select column title containing mutations
--in-pdb=<pdb>     PDB file of protein wildtype
--chain=<chain>    Select chain upon which to perform mutation
--pH=<ph>          Set pH of the protein
--no_restraints    Allow movement of all atoms
"""

import logging
from docopt import docopt
import pandas as pd
import numpy as np
import re
import sys
import pathlib
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
from sys import stdout
from pdbfixer import PDBFixer
from openmm.app import *
from openmm import *
from openmm.unit import *

def get_data(csvname: str, col: str):
    data = pd.read_csv(csvname)
    clipped_data=pd.DataFrame(data=data, columns=[col])
    divided_data = clipped_data[col].str.extract(r'([a-zA-Z]+)([^a-zA-Z]+)([a-zA-Z]+)')
    divided_cleaned = divided_data.to_string(header=None, index=False)
    amino_index = {'G': 'GLY' , 'L': 'LEU', 'R': 'ARG', 'A': 'ALA', 'Y': 'TYR', 'S': 'SER', 'M': 'MET', 'I': 'ILE', 'T': 'THR', 'C': 'CYS', 'V': 'VAL', 'P': 'PRO', 'F': 'PHE', 'W': 'TRP', 'H': 'HIS', 'K': 'LYS', 'D': 'ASP', 'E': 'GLU', 'N': 'ASN', 'Q': 'GLN', ' ': '-'}
    new_rep_codes = re.sub(r"[GLRAYSMITCVPFWHKDENQ ]", lambda x: amino_index[x.group(0)], divided_cleaned)
    new_rep_cleaned = re.sub("--","-", new_rep_codes)
    return new_rep_cleaned

def reformat_data(csvname: str, new_rep_cleaned):
    df_newcol = new_rep_cleaned.splitlines()
    pre_df = {'name': df_newcol}
    df1 = pd.DataFrame(pre_df)
    df2 = pd.read_csv(csvname)
    df3 = pd.concat([df1, df2['ΔΔG']], axis=1)
    stem = csvname.replace(".csv","")
    csv_reformat = str(stem + "_reformat.csv")
    df3.to_csv(csv_reformat)

def clean_wildtype(pdbname: str, pH: str):
    pH_fl = float(pH)
    pdb = PDBFixer(pdbname)
    #pdb.findMissingResidues()
    pdb.missingResidues = {}
    #pdb.findNonstandardResidues()
    #pdb.replaceNonstandardResidues()
    pdb.removeHeterogens(False)
    pdb.findMissingAtoms()
    pdb.addMissingAtoms()
    pdb.addMissingHydrogens(pH_fl)
    PDBFile.writeFile(pdb.topology, pdb.positions, open("wildtype_fixed.pdb", 'w'), keepIds=True)
    return pdb

def create_mutants(pdbname: str, mutant, chain: str, pH: str):
    pH_fl = float(pH)
    mutpdb = PDBFixer(pdbname)
    mutpdb.applyMutations([mutant], chain)
    #mutpdb.findMissingResidues()
    mutpdb.missingResidues = {}
    #mutpdb.findNonstandardResidues()
    #mutpdb.replaceNonstandardResidues()
    mutpdb.removeHeterogens(False)
    mutpdb.findMissingAtoms()
    mutpdb.addMissingAtoms()
    mutpdb.addMissingHydrogens(pH_fl)
    PDBFile.writeFile(mutpdb.topology, mutpdb.positions, open(mutant + "_fixed.pdb", 'w'), keepIds=True)
    return mutpdb

def setup_forcefield():
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    return forcefield

def setup_modeller(pdb):
    modeller = Modeller(pdb.topology, pdb.positions)
    return modeller

def setup_system(modeller, forcefield, solvmol: str, no_restraints: bool):
    Natoms=int(solvmol)
    temp = 300
    Rgas = 8.205736608e-5
    Na = 6.02214076e23
    Nmol = Natoms/Na
    p = 1.0
    Vol = (Nmol*Rgas*temp)/p
    x = Vol**(1./3.)
    y = Vol**(1./3.)
    z = Vol**(1./3.)
    modeller.topology.setUnitCellDimensions((x, y, z))
    modeller.addSolvent(forcefield, numAdded=solvmol)
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1.0*nanometer, constraints=HBonds)
    if not no_restraints:
        logging.info("Using restraints on backbone")
        for atom in modeller.topology.atoms():
            if atom.name in ["CA","C","N"]:
                system.setParticleMass(atom.index, 0)
    return system

def setup_simulation(modeller, system):
    integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.001*picoseconds)
    simulation = Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)
    return simulation

def energy_minimization(modeller):
    forcefield = setup_forcefield()
    solvmol = 3000
    system = setup_system(modeller, forcefield, solvmol, arguments['--no_restraints'])
    simulation = setup_simulation(modeller, system)
    init_state = simulation.context.getState(getEnergy=True, getPositions=True)
    logging.info("Starting potential energy = %.9f kcal/mol"
        % init_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole))
    simulation.minimizeEnergy
    final_state = simulation.context.getState(getEnergy=True, getPositions=True)
    logging.info("Starting potential energy = %.9f kcal/mol"
        % final_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole))
    final_pe = final_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole)
    return simulation.topology, final_state.getPositions(), final_pe, simulation, 

def md_nvt(simulation, csvname: str, totalsteps: int, reprate: int, pdbname):
    init_state = simulation.context.getState(getEnergy=True, getPositions=True)
    init_pe = init_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole)
    simulation.reporters.append(PDBReporter(pdbname, reprate))
    simulation.reporters.append(StateDataReporter(stdout, reprate, step=True, potentialEnergy=True, temperature=True, volume=True))
    prepdf = {'Step':[], 'Potential Energy (kJ/mole)':[], 'Temperature (K)':[], 'Box Volume (nm^3)':[]}
    inidf = pd.DataFrame(prepdf)
    inidf.to_csv(csvname, index=False)
    simulation.reporters.append(StateDataReporter(csvname, reprate, step=True,
        potentialEnergy=True, temperature=True, volume=True, append=True))
    simulation.step(totalsteps)
    final_state = simulation.context.getState(getEnergy=True, getPositions=True)
    final_pe = final_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole)
    df = pd.read_csv(csvname)
    av_energy = df.loc[:, 'Potential Energy (kJ/mole)'].mean()
    return init_pe, final_pe, av_energy, simulation, simulation.topology, init_state.getPositions()

#def rmsd_analysis(pdb_traj: str, pdb_ref: str):
#    trajectory = md.load_pdb(pdb_traj)
#    reference = md.load_pdb(pdb_ref)
#    trajectory.superpose(reference)
#    n_residues = trajectory.topology.n_residues
#    rmsd_per_residue = []
#    for i in range(n_residues):
#        atom_indices = trajectory.topology.residue(i).atom_indices
#        residue_rmsd = np.array([md.rmsd(trajectory, reference, frame=frame_index, atom_indices=atom_indices) for frame_index in range(trajectory.n_frames)])
#        rmsd_per_residue.append(residue_rmsd)

def rmsf_analysis(pdb_traj: str, rmsfcsv: str):
    u = mda.Universe(pdb_traj)
    average = align.AverageStructure(u, u, select='protein and name CA', ref_frame=0).run()
    ref = average.results.universe
    aligner = align.AlignTraj(u, ref, select='protein and name CA', in_memory=True).run()
    c_alphas = u.select_atoms('protein and name CA')
    R = rms.RMSF(c_alphas).run()
    residue_ids = [residue.resid for residue in c_alphas.residues]
    residue_names = [residue.resname for residue in c_alphas.residues]
    rmsf_values = R.rmsf
    df = pd.DataFrame({'Res_ID': residue_ids, 'Res_Name': residue_names, 'RMSF': rmsf_values })
    df.to_csv(rmsfcsv, index=False)

def main():
    arguments = docopt(__doc__, version='openmm-fb_MD.py')
    new_rep_cleaned = get_data(arguments['--incsv'], arguments['--from-col'])
    reformat_data(arguments['--incsv'], new_rep_cleaned)
    pdb = clean_wildtype(arguments['--in-pdb'], arguments['--pH'])
    for j in range(1,2):
        wt_pdb = energy_minimization(pdb)
        strj = str(j)
        wt_out = str("wt_minimised" + strj + ".pdb")
        PDBFile.writeFile(wt_pdb[0], wt_pdb[1], open(wt_out, "w"), keepIds=True)
        simulation = wt_pdb[3]
        csvname = str("wt_traj" + strj + ".csv")
        pdbname = str("wt_traj" + strj + ".pdb")
        sim_run = md_nvt(simulation, csvname, 2000000, 10000, pdbname)
        wt_ref = str("wt_reference" + strj + ".pdb")
        PDBFile.writeFile(sim_run[4], sim_run[5], open(wt_ref, "w"), keepIds=True)
        rmsfout = str("wt_rmsf" + strj + ".csv")
        rmsf_analysis(pdbname, rmsfout)

if __name__ == '__main__':
    arguments = docopt(__doc__, version='openmm-fb_MD.py')
    logging.getLogger().setLevel(logging.INFO)  
    main()