#!/usr/bin/env python3
"""openmmMD_twophase

Minimise the potential energy of the wildtype protein structure and mutated variants, according to AMBER14 potentials 
Usage:
openmmMD_twophase.py [--inpdb=<pdb>] [--pH=<pH>] [--no_restraints] 

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

def setup_pdbfile(pdbfile):
    pdb = app.PDBFile(pdbfile)
    return pdb

def setup_system(modeller, forcefield, solvmol: str, no_restraints: bool):
    Natoms=int(solvmol)
    temp = 300
    Rgas = 8.205736608e-5
    Na = 6.02214076e23
    Nmol = Natoms/Na
    p = 1.0
    Vol = (Nmol*Rgas*temp)/p
    # Compute the center of mass of the protein
    #positions = modeller.positions
    #num_atoms = len(positions)
    #center_of_mass = mm.Vec3(0.0, 0.0, 0.0)*nanometer
    #for pos in positions:
    #    center_of_mass += pos
    #center_of_mass /= num_atoms
    #center_of_mass = center_of_mass*nanometer
    # Translate positions to center at (0,0,0)
    #new_positions = [(pos - center_of_mass) for pos in positions]
    #new_positions = [(mm.Vec3(pos.x, pos.y, pos.z) - center_of_mass) for pos in positions]
    #modeller.positions = new_positions*nanometer
    #x = Vol**(1./3.) 
    #y = Vol**(1./3.) 
    #z = Vol**(1./3.)
    #x = 130
    #y = 130
    #z = 130
    #xvec = mm.Vec3(11.44, 0.0, 0.0)
    #yvec = mm.Vec3(0.0, 11.44, 0.0)
    #zvec = mm.Vec3(0.0, 0.0, 11.44)
    xvec = mm.Vec3(11.70, 0.0, 0.0)
    yvec = mm.Vec3(0.0, 11.70, 0.0)
    zvec = mm.Vec3(0.0, 0.0, 11.70)
    modeller.topology.setPeriodicBoxVectors((xvec,yvec,zvec))
    #modeller.addSolvent(forcefield, numAdded=Natoms, padding=1.0*nanometers)
    modeller.addSolvent(forcefield, numAdded=Natoms)
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.PME, nonbondedCutoff=1.0*nanometer, constraints=app.HBonds)
    #if not no_restraints:
    #    logging.info("Using restraints on backbone")
    #    for atom in modeller.topology.atoms():
    #        if atom.name in ["CA","C","N"]:
    #            system.setParticleMass(atom.index, 0)
    return system

def setup_system_in_box(pdb, forcefield):
    x = 130
    y = 130
    z = 130
    pdb.topology.setUnitCellDimensions((x, y, z))
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.PME, nonbondedCutoff=1.0*nanometer, constraints=app.HBonds)
    return system

def setup_simulation(modeller, system):
    #integrator = mm.LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.0001*picoseconds)
    integrator = mm.LangevinMiddleIntegrator(310*kelvin, 1/picosecond, 0.001*picoseconds)
    #platform = mm.Platform.getPlatformByName('CUDA')
    simulation = app.Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)
    #simulation.context.setPeriodicBoxVectors(*modeller.topology.getUnitCellVectors())
    return simulation, integrator

def setup_simulation_in_box(pdb, system):
    integrator = mm.LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.001*picoseconds)
    #platform = mm.Platform.getPlatformByName('CUDA')
    simulation = app.Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    #simulation.context.setPeriodicBoxVectors(*modeller.topology.getUnitCellVectors())
    return simulation, integrator

def energy_minimization(modeller, info_sheet):
    forcefield = setup_forcefield()
    solvmol = 50000
    system = setup_system(modeller, forcefield, solvmol, arguments['--no_restraints'])
    xvec = mm.Vec3(11.70, 0.0, 0.0)
    yvec = mm.Vec3(0.0, 11.70, 0.0)
    zvec = mm.Vec3(0.0, 0.0, 11.70)
    #xvec = mm.Vec3(11.44, 0.0, 0.0)
    #yvec = mm.Vec3(0.0, 11.44, 0.0)
    #zvec = mm.Vec3(0.0, 0.0, 11.44)
    simulation = setup_simulation(modeller, system)[0]
    integrator = setup_simulation(modeller, system)[1]
    simulation.context.setPeriodicBoxVectors(xvec,yvec,zvec)
    init_state = simulation.context.getState(getEnergy=True, getPositions=True)
    logging.info("Starting potential energy = %.9f kcal/mol"
        % init_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole))
    init_pe = init_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole)
    #simulation.context.setPeriodicBoxVectors(xvec,yvec,zvec)
    mm.LocalEnergyMinimizer.minimize(simulation.context)
    #simulation.minimizeEnergy()
    #simulation.context.setPeriodicBoxVectors(xvec,yvec,zvec)
    final_state = simulation.context.getState(getEnergy=True, getPositions=True)
    logging.info("Final potential energy = %.9f kcal/mol"
        % final_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole))
    final_pe = final_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole)
    with open(info_sheet, 'w') as file:
        file.write(f'Initial potential energy: {init_pe}\n')
        file.write(f'Minimised potential energy: {final_pe}\n')
    return simulation.topology, final_state.getPositions(), final_pe, simulation, integrator, system

def energy_minimization_in_box(pdb):
    forcefield = setup_forcefield()
    system = setup_system_in_box(pdb, forcefield)
    simulation = setup_simulation_in_box(pdb, system)[0]
    integrator = setup_simulation_in_box(pdb, system)[1]
    init_state = simulation.context.getState(getEnergy=True, getPositions=True)
    logging.info("Starting potential energy = %.9f kcal/mol"
        % init_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole))
    #simulation.minimizeEnergy()
    final_state = simulation.context.getState(getEnergy=True, getPositions=True)
    logging.info("Starting potential energy = %.9f kcal/mol"
        % final_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole))
    final_pe = final_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole)
    return simulation.topology, final_state.getPositions(), final_pe, simulation, integrator

def md_nvt(simulation, csvname: str, totalsteps: int, reprate: int, pdbname, integrator, system):
    init_state = simulation.context.getState(getEnergy=True, getPositions=True)
    init_pe = init_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole)
    for temp in range(0, 310, 10):
        integrator.setTemperature(temp*kelvin)
        logging.info(temp)
        simulation.reporters.append(app.StateDataReporter(stdout, 100, step=True, potentialEnergy=True, temperature=True, volume=True))
        simulation.step(10)
        #simulation.step(100)
    simulation.context.setVelocitiesToTemperature(310*kelvin)
    #simulation.reporters.append(app.StateDataReporter(stdout, reprate, step=True, potentialEnergy=True, temperature=True, volume=True))
    #simulation.step(496800)
    simulation.reporters.append(app.PDBReporter(pdbname, reprate))
    simulation.reporters.append(app.StateDataReporter(stdout, reprate, step=True, potentialEnergy=True, temperature=True, volume=True))
    prepdf = {'Step':[], 'Potential Energy_kJ/mole':[], 'Temperature_K)':[], 'Box Volume_nm^3':[]}
    inidf = pd.DataFrame(prepdf)
    inidf.to_csv(csvname, index=False)
    simulation.reporters.append(app.StateDataReporter(csvname, reprate, step=True,
        potentialEnergy=True, temperature=True, volume=True, append=True))
    simulation.step(totalsteps)
    final_state = simulation.context.getState(getEnergy=True, getPositions=True)
    final_pe = final_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole)
    df = pd.read_csv(csvname)
    av_energy = df.loc[:, 'Potential Energy_kJ/mole'].mean()
    return init_pe, final_pe, av_energy, simulation, simulation.topology, init_state.getPositions(), final_state.getPositions()

def md_nvt_reset(simulation, csvname: str, totalsteps: int, reprate: int, pdbname):
    init_state = simulation.context.getState(getEnergy=True, getPositions=True)
    init_pe = init_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole)
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


def rmsf_analysis(pdb_traj: str, rmsfcsv: str):
    u = mda.Universe(pdb_traj)
    #u.trajectory.dt = 0.001
    average = align.AverageStructure(u, u, select='protein and name CA', ref_frame=0).run()
    ref = average.results.universe
    align.AlignTraj(u, ref, select='protein and name CA', in_memory=True).run()
    c_alphas = u.select_atoms('protein and name CA')
    R = rms.RMSF(c_alphas).run()
    residue_ids = [residue.resid for residue in c_alphas.residues]
    residue_names = [residue.resname for residue in c_alphas.residues]
    rmsf_values = R.rmsf
    df = pd.DataFrame({'Res_ID': residue_ids, 'Res_Name': residue_names, 'RMSF': rmsf_values })
    print(df)
    df.to_csv(rmsfcsv, index=False)

def main():
    arguments = docopt(__doc__, version='openmmMD_twophase.py')
    pdb_clean = clean_pdb(arguments['--inpdb'], arguments['--pH'])
    pdb = pdb_clean[0]
    modeller = setup_modeller(pdb)
    for j in range(1,2):
        stem = pdb_clean[1]
        info_sheet = str(stem + "_info.txt")
        min_pdb = energy_minimization(modeller, info_sheet)
        strj = str(j)
        min_out = str(stem + strj + ".pdb")
        app.PDBFile.writeFile(min_pdb[0], min_pdb[1], open(min_out, "w"), keepIds=True)
        simulation = min_pdb[3]
        integrator = min_pdb[4]
        system = min_pdb[5]
        #csvinit = str(stem + "_initial" + strj + ".csv")
        #pdbinit = str(stem + "_initial" + strj + ".pdb")
        csvname = str(stem + "_traj" + strj + ".csv")
        pdbname = str(stem + "_traj" + strj + ".pdb")
        sim_run = md_nvt(simulation, csvname, 200, 10, pdbname, integrator, system)
        #sim_run = md_nvt(simulation, csvname, 2000000, 1000, pdbname, integrator, system)
        sim_ref = str(stem + "_reference" + strj + ".pdb")
        app.PDBFile.writeFile(sim_run[4], sim_run[5], open(sim_ref, "w"), keepIds=True)
        sim_fin = str(stem + "_feedstock" + strj + ".pdb")
        app.PDBFile.writeFile(sim_run[4], sim_run[6], open(sim_fin, "w"), keepIds=True)
        #pdb_new = setup_pdbfile(sim_fin)
        #new_min = energy_minimization_in_box(pdb_new)
        #new_simulation = new_min[3]
        #sim_runII = md_nvt_reset(new_simulation, csvname, 2000, 1000, pdbname)
        rmsfout = str(stem + "_rmsf" + strj + ".csv")
        rmsf_analysis(pdbname, rmsfout)
     

if __name__ == '__main__':
    arguments = docopt(__doc__, version='openmmMD_twophase.py')
    logging.getLogger().setLevel(logging.INFO)  
    main()