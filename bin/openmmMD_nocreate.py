#!/usr/bin/env python3
"""openmmMD_nocreate

Minimise the potential energy of the wildtype protein structure and mutated variants, according to AMBER14 potentials 
Usage:
openmmMD_nocreate.py [--inpdb=<pdb>] [--pH=<pH>] [--no_restraints] 

Options:
--inpdb=<pdb>      Input PDB file of protein as obtained from previous process
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
import openmm.app as app
import openmm as mm
#from openmm.app import *
#from openmm import *
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

def clean_pdb(pdbname: str, pH: str):
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
    stem = pdbname.replace(".pdb","")
    app.PDBFile.writeFile(pdb.topology, pdb.positions, open(stem + "_fixed.pdb", 'w'), keepIds=True)
    return pdb, stem

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
    app.PDBFile.writeFile(mutpdb.topology, mutpdb.positions, open(mutant + "_fixed.pdb", 'w'), keepIds=True)
    return mutpdb

def create_mutants_in_box(pdbname: str, mutant, chain: str, pH: str):
    pH_fl = float(pH)
    mutpdb = PDBFixer(pdbname)
    mutpdb.applyMutations([mutant], chain)
    mutpdb.missingResidues = {}
    mutpdb.removeHeterogens(True)
    mutpdb.findMissingAtoms()
    mutpdb.addMissingAtoms()
    mutpdb.addMissingHydrogens(pH_fl)
    app.PDBFile.writeFile(mutpdb.topology, mutpdb.positions, open(mutant + "_fixed.pdb", 'w'), keepIds=True)
    return mutpdb

def setup_forcefield():
    forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    return forcefield

def setup_modeller(pdb):
    #pdb = app.PDBFile(pdb)
    modeller = app.Modeller(pdb.topology, pdb.positions)
    return modeller

def setup_system(modeller, forcefield, solvmol: str, no_restraints: bool):
    Natoms=int(solvmol)
    temp = 300
    Rgas = 8.205736608e-5
    Na = 6.02214076e23
    Nmol = Natoms/Na
    p = 1.0
    Vol = (Nmol*Rgas*temp)/p
    x = Vol**(1./3.) + 3
    y = Vol**(1./3.) + 3
    z = Vol**(1./3.) + 3
    modeller.topology.setUnitCellDimensions((x, y, z))
    modeller.addSolvent(forcefield, padding=1.0*nanometers)
    #modeller.addSolvent(forcefield, numAdded=Natoms)
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.PME, nonbondedCutoff=1.0*nanometer, constraints=app.HBonds)
    #if not no_restraints:
    #    logging.info("Using restraints on backbone")
    #    for atom in modeller.topology.atoms():
    #        if atom.name in ["CA","C","N"]:
    #            system.setParticleMass(atom.index, 0)
    return system

def setup_system_in_box(modeller, forcefield):
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.PME, nonbondedCutoff=1.0*nanometer, constraints=app.HBonds)
    return system

def setup_simulation(modeller, system):
    integrator = mm.LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.0001*picoseconds)
    #platform = mm.Platform.getPlatformByName('CUDA')
    simulation = app.Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)
    #simulation.context.setPeriodicBoxVectors(*modeller.topology.getUnitCellVectors())
    return simulation, integrator

def energy_minimization(modeller):
    forcefield = setup_forcefield()
    solvmol = 1000
    system = setup_system(modeller, forcefield, solvmol, arguments['--no_restraints'])
    #system = system_all[0]
    #x = system_all[1]
    #y = system_all[2]
    #z = system_all[3]
    #system = setup_system(modeller, forcefield, solvmol, arguments['--no_restraints'])[0]
    #x = setup_system(modeller, forcefield, solvmol, arguments['--no_restraints'])[1]
    #y = setup_system(modeller, forcefield, solvmol, arguments['--no_restraints'])[2]
    #z = setup_system(modeller, forcefield, solvmol, arguments['--no_restraints'])[3]
    #system = setup_system(modeller, forcefield, solvmol, arguments['--no_restraints'])[0]
    #xvec = mm.Vec3(x, 0.0, 0.0)
    #yvec = mm.Vec3(0.0, y, 0.0)
    #zvec = mm.Vec3(0.0, 0.0, z)
    simulation = setup_simulation(modeller, system)[0]
    integrator = setup_simulation(modeller, system)[1]
    init_state = simulation.context.getState(getEnergy=True, getPositions=True)
    logging.info("Starting potential energy = %.9f kcal/mol"
        % init_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole))
    #simulation.context.setPeriodicBoxVectors(xvec,yvec,zvec)
    #mm.LocalEnergyMinimizer.minimize(simulation.context)
    simulation.minimizeEnergy()
    #simulation.context.setPeriodicBoxVectors(xvec,yvec,zvec)
    final_state = simulation.context.getState(getEnergy=True, getPositions=True)
    logging.info("Final potential energy = %.9f kcal/mol"
        % final_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole))
    final_pe = final_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole)
    return simulation.topology, final_state.getPositions(), final_pe, simulation, integrator, system

def energy_minimization_in_box(modeller):
    forcefield = setup_forcefield()
    system = setup_system_in_box(modeller, forcefield)
    simulation = setup_simulation(modeller, system)[0]
    integrator = setup_simulation(modeller, system)[1]
    init_state = simulation.context.getState(getEnergy=True, getPositions=True)
    logging.info("Starting potential energy = %.9f kcal/mol"
        % init_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole))
    simulation.minimizeEnergy()
    final_state = simulation.context.getState(getEnergy=True, getPositions=True)
    logging.info("Starting potential energy = %.9f kcal/mol"
        % final_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole))
    final_pe = final_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole)
    return simulation.topology, final_state.getPositions(), final_pe, simulation, integrator

def md_nvt(simulation, csvname: str, totalsteps: int, reprate: int, pdbname, integrator, system):
    init_state = simulation.context.getState(getEnergy=True, getPositions=True)
    init_pe = init_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole)
    for temp in range(0, 300, 10):
        integrator.setTemperature(temp*kelvin)
        simulation.step(10)
    simulation.context.setVelocitiesToTemperature(300*kelvin)
    #simulation.step(350000)
    new_integrator = mm.LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.001*picoseconds)
    #context = simulation.context
    #simulation.context.reinitialize(preserveState=True)
    #simulation.integrator = new_integrator
    #context.reinitialize(preserveState=True)
    state = simulation.context.getState(getPositions=True, getVelocities=True, getEnergy=True)
    simulation.integrator = new_integrator
    simulation.context = mm.Context(system, new_integrator)
    simulation.context.setPositions(state.getPositions())
    simulation.context.setVelocities(state.getVelocities())
    #simulation.context.setIntegrator(new_integrator)
    simulation.step(49690)
    simulation.reporters.append(app.PDBReporter(pdbname, reprate))
    simulation.reporters.append(app.StateDataReporter(stdout, reprate, step=True, potentialEnergy=True, temperature=True, volume=True))
    prepdf = {'Step':[], 'Potential Energy (kJ/mole)':[], 'Temperature (K)':[], 'Box Volume (nm^3)':[]}
    inidf = pd.DataFrame(prepdf)
    inidf.to_csv(csvname, index=False)
    simulation.reporters.append(app.StateDataReporter(csvname, reprate, step=True,
        potentialEnergy=True, temperature=True, volume=True, append=True))
    simulation.step(totalsteps)
    final_state = simulation.context.getState(getEnergy=True, getPositions=True)
    final_pe = final_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole)
    df = pd.read_csv(csvname)
    av_energy = df.loc[:, 'Potential Energy (kJ/mole)'].mean()
    return init_pe, final_pe, av_energy, simulation, simulation.topology, init_state.getPositions(), final_state.getPositions()

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

def main():
    arguments = docopt(__doc__, version='openmmMD_nocreate.py')
    pdb_clean = clean_pdb(arguments['--inpdb'], arguments['--pH'])
    pdb = pdb_clean[0]
    modeller = setup_modeller(pdb)
    for j in range(1,2):
        min_pdb = energy_minimization(modeller)
        strj = str(j)
        stem = pdb_clean[1]
        min_out = str(stem + strj + ".pdb")
        app.PDBFile.writeFile(min_pdb[0], min_pdb[1], open(min_out, "w"), keepIds=True)
        simulation = min_pdb[3]
        integrator = min_pdb[4]
        system = min_pdb[5]
        csvname = str(stem + "_traj" + strj + ".csv")
        pdbname = str(stem + "_traj" + strj + ".pdb")
        sim_run = md_nvt(simulation, csvname, 50000, 100, pdbname, integrator, system)
        sim_ref = str(stem + "_reference" + strj + ".pdb")
        app.PDBFile.writeFile(sim_run[4], sim_run[5], open(sim_ref, "w"), keepIds=True)
        sim_fin = str(stem + "_final" + strj + ".pdb")
        app.PDBFile.writeFile(sim_run[4], sim_run[6], open(sim_fin, "w"), keepIds=True)
        rmsfout = str(stem + "_rmsf" + strj + ".csv")
        rmsf_analysis(pdbname, rmsfout)
     

if __name__ == '__main__':
    arguments = docopt(__doc__, version='openmmMD_nocreate.py')
    logging.getLogger().setLevel(logging.INFO)  
    main()