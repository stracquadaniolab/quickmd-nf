#!/usr/bin/env python3
"""openmmMD_dualmin_restraint

Minimise the potential energy of the wildtype protein structure and mutated variants, according to AMBER14 potentials 
Usage:
openmmMD_dualmin_restraint.py [--inpdb=<pdb>] [--pH=<pH>] [--steps=<steps>] [--report_every=<report_every>] [--init_steps=<init_steps>] [--relax_steps=<relax_steps>] [--relax_report_every=<relax_report_every>] [--no_restraints] 

Options:
--inpdb=<pdb>                               Input PDB file of protein as obtained from previous process
--pH=<ph>                                   Set pH of the protein
--steps=<steps>                             Number of steps in final simulation proper after equilibriation and relax, used for RMSF calculation
--report_every=<report_every>               Report every # steps to final output PDB trajectory and energetic CSV file (CSV for whole reported trajectory)
--init_steps=<init_steps>                   Initial steps of solvent equilibration (fixed protein backbone) before any reported trajectory 
--relax_steps=<relax_steps>                 Initial reported trajectory steps focused on protein structure relax without restraints
--relax_report_every=<relax_report_every>   Structure sampling rate of initial reported trajectory 
--no_restraints                             Allow movement of all atoms
"""

import logging
from docopt import docopt
import pandas as pd
import numpy as np
import re
import sys
import pathlib
import MDAnalysis as mda
from MDAnalysis import transformations
from MDAnalysis.analysis import rms, align
from sys import stdout
from pdbfixer import PDBFixer
import openmm.app as app
import openmm as mm
from openmm.unit import *
from openmm.app import element

#Process to clean the original PDB, find missing atoms, remove heterogens and assign a protonation state to mimic selected pH, and output fixed PDB file
def clean_pdb(pdbname: str, pH: str):
    pH_fl = float(pH)
    pdb = PDBFixer(pdbname)
    pdb.missingResidues = {}
    pdb.removeHeterogens(False)
    pdb.findMissingAtoms()
    pdb.addMissingAtoms()
    pdb.addMissingHydrogens(pH_fl)
    stem = pdbname.replace(".pdb","")
    app.PDBFile.writeFile(pdb.topology, pdb.positions, open(stem + "_fixed.pdb", 'w'), keepIds=True)
    return pdb, stem

#Carries out minimisation on protein in vacuum to optimise not only covalent bonding structure but also important hydrogen bonds responsable for secondary structure. PDB file and final energy are outputted at the end. Similation parameters are self-contained in the process to avoid interference with later setup of solvation box.
def vacuum_minimization(pdb, pdbout: str):
    pdb = pdb
    forcefield = app.ForceField('amber14-all.xml')
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff)
    integrator = mm.LangevinMiddleIntegrator(310*kelvin, 1/picosecond, 0.001*picoseconds)
    simulation = app.Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    simulation.minimizeEnergy()
    final_state = simulation.context.getState(getEnergy=True, getPositions=True)
    logging.info("Final potential energy = %.9f kcal/mol"
        % final_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole))
    final_pe = final_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole)
    app.PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), open(pdbout, "w"), keepIds=True)
    return pdbout, final_pe

#Similar to the above, carries out minimisation on protein in vacuum to optimise not only covalent bonding structure but also important hydrogen bonds responsable for secondary structure. However coordinates, instead of PDB file, along with final energy are outputted at the end, to allow sampling of minimised states and select the most optimised. Similation parameters are self-contained in the process to avoid interference with later setup of solvation box.
def vacuum_minimization_sampling(pdb):
    pdb = pdb
    forcefield = app.ForceField('amber14-all.xml')
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff)
    integrator = mm.LangevinMiddleIntegrator(310*kelvin, 1/picosecond, 0.001*picoseconds)
    simulation = app.Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    simulation.minimizeEnergy()
    final_state = simulation.context.getState(getEnergy=True, getPositions=True)
    logging.info("Final potential energy = %.9f kcal/mol"
        % final_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole))
    final_pe = final_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole)
    positions_final = final_state.getPositions()
    return final_pe, positions_final 

#Set up forcefield with standard AMBER parameters for protein structure, and TIP3PFB
def setup_forcefield():
    forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    return forcefield

#Set up input PDB topology and positions as modeller such that the simulation box can be modified (add solvent, virtual sites etc)
def setup_modeller(pdb):
    #pdb = app.PDBFile(pdb)
    modeller = app.Modeller(pdb.topology, pdb.positions)
    return modeller

#Read PDB file
def setup_pdbfile(pdbfile):
    pdb = app.PDBFile(pdbfile)
    return pdb

#instate resraints on backbone atoms via rendering backbone massless
def apply_mass_restraints(system, topology):
    for atom in topology.atoms():
        if atom.name in ["CA", "C", "N"]:
            system.setParticleMass(atom.index, 0)

#disengage restraints on backbone atoms via reinstating default atomic mass
def remove_mass_restraints(system, topology):
    for atom in topology.atoms():
        if atom.name in ["CA", "C", "N"]:
            element_mass = element.Element.getBySymbol(atom.element.symbol).mass
            system.setParticleMass(atom.index, element_mass)


#integrate modeller and forcefied into a defined system. Before integration, solvent molecules are added to the box, with box vectors defined such that water density at a given pressure (normally standard pressure of 1atm) can be simulated. Padding is added, with excess molecules removed to achieve 1atm density
def setup_system(modeller, forcefield, solvmol: str, padding: int):
    Natoms=int(solvmol)
    xvec = mm.Vec3(11.70, 0.0, 0.0)
    yvec = mm.Vec3(0.0, 11.70, 0.0)
    zvec = mm.Vec3(0.0, 0.0, 11.70)
    modeller.topology.setPeriodicBoxVectors((xvec,yvec,zvec))
    padding = padding * nanometer
    modeller.addSolvent(forcefield, padding=padding)
    water_molecules = []
    for res in modeller.topology.residues():
        if res.name == 'HOH':
            water_molecules.append(res)
    current_water_count = len(water_molecules)
    excess_water_count = current_water_count - Natoms
    if excess_water_count > 0:
        residues_to_remove = water_molecules[:excess_water_count]
        modeller.delete(residues_to_remove)
    final_water_count = sum(1 for res in modeller.topology.residues() if res.name == 'HOH')
    logging.info("Final water count = %d ", final_water_count)
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.PME, nonbondedCutoff=1.0*nanometer, constraints=app.HBonds)
    #if not no_restraints:
    #    logging.info("Using restraints on backbone")
    #    for atom in modeller.topology.atoms():
    #        if atom.name in ["CA","C","N"]:
    #            system.setParticleMass(atom.index, 0)
    return system

#Define simulation parameter and simulation.context, combining modeller and system, and setting integrator 
def setup_simulation(modeller, system):
    integrator = mm.LangevinMiddleIntegrator(310*kelvin, 1/picosecond, 0.001*picoseconds)
    simulation = app.Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)
    return simulation, integrator

#Combine setup_forcefield and setup_system functions, redefine box vectors ahead of local minimisation. Output final minimised topology and mimimisation, and return simulation after minimisation 
def energy_minimization(modeller, info_sheet, no_restraints: bool):
    forcefield = setup_forcefield()
    solvmol = 50000
    padding = 1 #nanometer
    system = setup_system(modeller, forcefield, solvmol, padding)
    xvec = mm.Vec3(11.70, 0.0, 0.0)
    yvec = mm.Vec3(0.0, 11.70, 0.0)
    zvec = mm.Vec3(0.0, 0.0, 11.70)
    if not no_restraints:
        apply_mass_restraints(system, modeller.topology)
    simulation = setup_simulation(modeller, system)[0]
    integrator = setup_simulation(modeller, system)[1]
    simulation.context.setPeriodicBoxVectors(xvec,yvec,zvec)
    init_state = simulation.context.getState(getEnergy=True, getPositions=True)
    logging.info("Starting potential energy = %.9f kcal/mol"
        % init_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole))
    init_pe = init_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole)
    mm.LocalEnergyMinimizer.minimize(simulation.context)
    final_state = simulation.context.getState(getEnergy=True, getPositions=True)
    logging.info("Final potential energy = %.9f kcal/mol"
        % final_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole))
    final_pe = final_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole)
    with open(info_sheet, 'w') as file:
        file.write(f'Initial potential energy: {init_pe}\n')
        file.write(f'Minimised potential energy: {final_pe}\n')
    return simulation.topology, final_state.getPositions(), final_pe, simulation, integrator, system

#Takes minimised solvation box to run grand canonical ensemble (NVT), with initial temperature increase in 100-step intervals. 0,5ns of equilibration is allowed to achieve equilibrium state. After which the requested number of steps is run, with reporting corresponding to the number of steps in the --report_every option
def md_nvt(simulation, csvname: str, totalsteps: int, reprate: int, pdbname, integrator, system, initsteps: int, initpdb: str, relax_steps: int, relax_report_every: int, no_restraints: bool):
    init_state = simulation.context.getState(getEnergy=True, getPositions=True)
    init_pe = init_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole)
    for temp in range(0, 310, 10):
        integrator.setTemperature(temp*kelvin)
        simulation.reporters.append(app.StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True, volume=True))
        simulation.step(10000)
    simulation.context.setVelocitiesToTemperature(310*kelvin)
    #simulation.step(496800)
    simulation.step(initsteps)
    if not no_restraints:
        remove_mass_restraints(system, simulation.topology)
    simulation.reporters.append(app.PDBReporter(initpdb, relax_report_every))
    simulation.reporters.append(app.StateDataReporter(stdout, reprate, step=True, potentialEnergy=True, temperature=True, volume=True))
    prepdf = {'Step':[], 'Potential Energy_kJ/mole':[], 'Temperature_K)':[], 'Box Volume_nm^3':[]}
    inidf = pd.DataFrame(prepdf)
    inidf.to_csv(csvname, index=False)
    simulation.reporters.append(app.StateDataReporter(csvname, reprate, step=True,
        potentialEnergy=True, temperature=True, volume=True, append=True))
    simulation.step(relax_steps)
    simulation.reporters.append(app.PDBReporter(pdbname, reprate))
    simulation.step(totalsteps)
    final_state = simulation.context.getState(getEnergy=True, getPositions=True)
    final_pe = final_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole)
    df = pd.read_csv(csvname)
    av_energy = df.loc[:, 'Potential Energy_kJ/mole'].mean()
    return init_pe, final_pe, av_energy, simulation, simulation.topology, init_state.getPositions(), final_state.getPositions()

#output rmsf per residue per chain based on generated PDB trajectory 
def rmsf_analysis_by_chain(pdb_traj: str, rmsfcsv: str):
    u = mda.Universe(pdb_traj)
    #ag = u.atoms
    #new_dimensions = [117, 117, 117, 90, 90, 90]
    #set_dim = transformations.boxdimensions.set_dimensions(new_dimensions)
    #transform = transformations.unwrap(ag)
    #center = transformations.center_in_box(ag.select_atoms('protein'), wrap=True)
    #u.trajectory.add_transformations(set_dim, transform, center)
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
        rmsfout = str(rmsfcsv + "_rmsf_ch" + chain_id + ".csv")
        df.to_csv(rmsfout, index=False)

def main():
    arguments = docopt(__doc__, version='openmmMD_dualmin_restraint.py')
    #fix PDB input file and impose protonation state to recreate chosen pH 
    pdb_clean = clean_pdb(arguments['--inpdb'], arguments['--pH'])
    pdb = pdb_clean[0]
    stem = pdb_clean[1]
    pH_str = str(arguments['--pH'] + "pH")
    for j in range(1,2):
        strj = str(j)
        #sample 100 iterations of minimisation in a vacuum and select coordinates of lowest energy state  
        pdbout = str(stem + "_" + pH_str + "_iter_" + strj + "vac_min.pdb")
        energies = []
        positions = []
        iterations = 100
        for step in range(iterations):
            state_out = vacuum_minimization_sampling(pdb)
            energy = state_out[0]
            positions_final = state_out[1]
            energies.append(energy)
            positions.append(positions_final)
        min_energy_index = np.argmin(energies)
        lowest_energy_positions = positions[min_energy_index]
        app.PDBFile.writeFile(pdb.topology, lowest_energy_positions, open(pdbout, 'w'), keepIds=True)
        pdb = app.PDBFile(pdbout)
        max_iterations = 100
        #set up simulation of protein solvated by water, conditional on minimisation creating net potential energy state (E_tot < 0) to avoid box explosion
        info_sheet = str(stem + pH_str + "iter" + strj + "_info.txt")
        modeller = setup_modeller(pdb)
        second_iteration = 0
        while second_iteration < max_iterations:
            min_pdb = energy_minimization(modeller, info_sheet, arguments['--no_restraints'])
            outputII = min_pdb[2]
            thresholdII = 0
            if outputII < thresholdII:
                break
            second_iteration += 1
        min_out = str(stem + pH_str + "iter" + strj + ".pdb")
        app.PDBFile.writeFile(min_pdb[0], min_pdb[1], open(min_out, "w"), keepIds=True)
        #set system, simulation and integrator parameters to run main trajectory, define output CSV and trajectory PDB names, set total steps, report rate and initial steps (for equilibration), set as inputs to NVT simulation function (md_nvt)
        simulation = min_pdb[3]
        integrator = min_pdb[4]
        system = min_pdb[5]
        initpdb = str(stem + pH_str + "_structure_sampling_" + strj + ".pdb")
        csvname = str(stem + pH_str + "_traj" + strj + ".csv")
        pdbname = str(stem + pH_str + "_traj" + strj + ".pdb")
        steps =int(arguments['--steps'])
        report_every = int(arguments['--report_every'])
        init_steps = int(arguments['--init_steps']) # total unreported step during solvent equilibration/protein restrained
        temp_steps = 310000 #total number of steps for heating from 0K to 310K
        equil_steps = init_steps - temp_steps
        relax_steps = int(arguments["--relax_steps"])
        relax_report_every = int(arguments['--relax_report_every'])
        sim_run = md_nvt(simulation, csvname, steps, report_every, pdbname, integrator, system, equil_steps, initpdb, relax_steps, relax_report_every, arguments['--no_restraints'])
        av_energy = sim_run[2]
        with open(info_sheet, 'a') as file:
            file.write(f'Average trajectory potential energy: {av_energy}\n')
        #Reference (final) and feedstock (initial) frames in case of future use
        sim_ref = str(stem + pH_str + "_reference" + strj + ".pdb")
        app.PDBFile.writeFile(sim_run[4], sim_run[5], open(sim_ref, "w"), keepIds=True)
        sim_fin = str(stem + pH_str + "_feedstock" + strj + ".pdb")
        app.PDBFile.writeFile(sim_run[4], sim_run[6], open(sim_fin, "w"), keepIds=True)
        #set rmsf output name and perform rmsf analysis on output trajectory
        rmsfstem = str(stem + pH_str + "_traj" + strj)
        rmsf_analysis_by_chain(pdbname, rmsfstem)
     

if __name__ == '__main__':
    arguments = docopt(__doc__, version='openmmMD_dualmin_restraint.py')
    logging.getLogger().setLevel(logging.INFO)  
    main()