#!/usr/bin/env python3
"""openmmMD_dualmin_restraint_NPT

Minimise the potential energy of the wildtype protein structure and mutated variants, according to AMBER14 potentials 
Usage:
openmmMD_dualmin_restraint_NPT.py [--inpdb=<pdb>] [--pH=<pH>] [--steps=<steps>] [--report_every=<report_every>] [--init_steps=<init_steps>] [--relax_steps=<relax_steps>] [--relax_report_every=<relax_report_every>] [--no_restraints] [--salting=<salting>] [--test_conditions]

Options:
--inpdb=<pdb>                               Input PDB file of protein as obtained from previous process
--pH=<ph>                                   Set pH of the protein
--steps=<steps>                             Number of steps in final simulation proper after equilibriation and relax, used for RMSF calculation
--report_every=<report_every>               Report every # steps to final output PDB trajectory and energetic CSV file (CSV for whole reported trajectory)
--init_steps=<init_steps>                   Initial steps of solvent equilibration (fixed protein backbone) before any reported trajectory 
--relax_steps=<relax_steps>                 Initial reported trajectory steps focused on protein structure relax without restraints
--relax_report_every=<relax_report_every>   Structure sampling rate of initial reported trajectory 
--no_restraints                             Allow movement of all atoms
--salting=<salting>                         Add an ion concentration in Molar
--test_conditions                           Apply test conditions, reducing initial phase significantly
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

# process to clean the original PDB, find missing atoms, remove heterogens and assign a protonation state to mimic selected pH, and output fixed PDB file
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

# carries out minimisation on protein in vacuum to optimise not only covalent bonding structure but also important hydrogen bonds responsable for secondary structure. PDB file and final energy are outputted at the end. Similation parameters are self-contained in the process to avoid interference with later setup of solvation box.
def vacuum_minimization(pdb, pdbout: str):
    pdb = pdb
    forcefield = app.ForceField('amber14-all.xml')
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff)
    integrator = mm.LangevinMiddleIntegrator(310*kelvin, 1/picosecond, 0.001*picoseconds)
    simulation = app.Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    simulation.minimizeEnergy(tolerance=0.5)
    final_state = simulation.context.getState(getEnergy=True, getPositions=True)
    logging.info("Final potential energy = %.9f kcal/mol"
        % final_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole))
    final_pe = final_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole)
    app.PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), open(pdbout, "w"), keepIds=True)
    return pdbout, final_pe

# similar to the above, carries out minimisation on protein in vacuum to optimise not only covalent bonding structure but also important hydrogen bonds responsable for secondary structure. However coordinates, instead of PDB file, along with final energy are outputted at the end, to allow sampling of minimised states and select the most optimised. Similation parameters are self-contained in the process to avoid interference with later setup of solvation box.
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

# set up forcefield with standard AMBER parameters for protein structure, and TIP3PFB
def setup_forcefield():
    forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    return forcefield

# set up input PDB topology and positions as modeller such that the simulation box can be modified (add solvent, virtual sites etc)
def setup_modeller(pdb):
    #pdb = app.PDBFile(pdb)
    modeller = app.Modeller(pdb.topology, pdb.positions)
    return modeller

# read PDB file
def setup_pdbfile(pdbfile):
    pdb = app.PDBFile(pdbfile)
    return pdb

# instate mass-based resraints on backbone atoms via rendering backbone massless
def apply_mass_restraints(system, topology):
    for atom in topology.atoms():
        if atom.name in ["CA", "C", "N"]:
            system.setParticleMass(atom.index, 0)

# disengage mass-based restraints on backbone atoms via reinstating default atomic mass
def remove_mass_restraints(system, topology):
    for atom in topology.atoms():
        if atom.name in ["CA", "C", "N"]:
            element_mass = element.Element.getBySymbol(atom.element.symbol).mass
            system.setParticleMass(atom.index, element_mass)

# instate force-based restraints 
def apply_force_restraints(system, modeller):
    positions = modeller.positions 
    restraint_force = mm.CustomExternalForce('0.5 * k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)')
    restraint_force.addGlobalParameter('k', 1000.0)  # strength of restraint
    restraint_force.addPerParticleParameter('x0')
    restraint_force.addPerParticleParameter('y0')
    restraint_force.addPerParticleParameter('z0')
    for atom in modeller.topology.atoms():
        if atom.name in ["CA", "C", "N"]:
            pos = positions[atom.index]
            restraint_force.addParticle(atom.index, [pos.x, pos.y, pos.z])
    system.addForce(restraint_force)

def remove_force_restraints(simulation):
    simulation.context.setParameter('k', 0.0)
    #restraint_force.setGlobalParameterDefaultValue(0, 0.0)


# integrate modeller and forcefied into a defined system. Before integration, solvent molecules are added to the box, with box vectors defined such that water density at a given pressure (normally standard pressure of 1atm) can be simulated. Padding is added, with excess molecules removed to achieve 1atm density
def setup_system(modeller, forcefield, solvmol: str, padding: int, salting=None):
    Natoms=int(solvmol)
    xvec = mm.Vec3(11.70, 0.0, 0.0)
    yvec = mm.Vec3(0.0, 11.70, 0.0)
    zvec = mm.Vec3(0.0, 0.0, 11.70)
    modeller.topology.setPeriodicBoxVectors((xvec,yvec,zvec))
    padding = padding * nanometer
    if salting:
        ion_conc = float(salting) * molar
        modeller.addSolvent(forcefield, padding=padding, ionicStrength=ion_conc)
    else:
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
    final_water_count = sum([1 for res in modeller.topology.residues() if res.name == 'HOH'])
    logging.info("Final water count = %d ", final_water_count)
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.PME, nonbondedCutoff=1.0*nanometer, constraints=app.HBonds)
    return system

# similar to the above, integrate modeller and forcefied into a defined system. Before integration, solvent molecules are added to the box, with box vectors defined such that water density at a given pressure (normally standard pressure of 1atm) can be simulated. Padding is added. Differently to the above, excess molecules removed proportionally, be they water or ion, to achieve 1atm density, with random selection to avoid density "holes"
def setup_system_proportional(modeller, forcefield, solvmol: str, padding: int, info_sheet,  salting=None):
    Nmol=int(solvmol)
    xvec = mm.Vec3(11.70, 0.0, 0.0)
    yvec = mm.Vec3(0.0, 11.70, 0.0)
    zvec = mm.Vec3(0.0, 0.0, 11.70)
    modeller.topology.setPeriodicBoxVectors((xvec,yvec,zvec))
    padding = padding * nanometer
    if salting:
        ion_conc = float(salting) * molar
        modeller.addSolvent(forcefield, padding=padding, ionicStrength=ion_conc)
    else:
        modeller.addSolvent(forcefield, padding=padding)
    # count up how much solvent we currently have
    # count up waters
    current_water = 0
    current_na = 0
    current_cl = 0

    for res in modeller.topology.residues():
        if res.name == 'HOH':
            current_water += 1
        elif res.name == 'NA':
            current_na += 1
        elif res.name == 'CL':
            current_cl += 1

    print('current water', current_water)
    print('current na', current_na)
    print('current cl', current_cl)
    current_total_solvent = current_water + current_na + current_cl
    # calculate the total number of counter ions balancing the protonation state (there should be slightly more Na+ than Cl- above pH 7 and the opposite below so we'll take the absolute amount) and therefore how much overall neutral solvent remains after these neutral ions are factored in. This will help us also determine the amount of target neutral solvent after removing excess:
    ion_imbalance = abs(current_cl - current_na)
    print('ion imbalance', ion_imbalance)
    current_neutral_solvent = current_total_solvent - ion_imbalance
    print('neutral solvent', current_neutral_solvent)
    target_neutral_solvent = Nmol - ion_imbalance
    print(' target ', target_neutral_solvent)
    excess_solvent = current_neutral_solvent - target_neutral_solvent
    print("excess" , excess_solvent)
    # now, work out the fraction of the neutral solvent that is water and that which is ions. We will use the lower ion score of the two to exclude the counter ions from the proportions. We'll then gather up our residues and randomly select those to delete 
    if excess_solvent > 0:
        water_fraction = current_water/current_neutral_solvent
        print('water fraction', water_fraction)
        target_water = round(target_neutral_solvent * water_fraction)
        print('target water', target_water)
        remove_water = current_water - target_water
    if current_cl > current_na:
        ion_fraction = current_na/current_neutral_solvent
        target_ions = round(target_neutral_solvent * ion_fraction)
        remove_ions = current_na - target_ions
    else:
        ion_fraction = current_cl/current_neutral_solvent
        target_ions = round(target_neutral_solvent * ion_fraction)
        remove_ions = current_cl - target_ions
    print('ion fraction', ion_fraction)
    print('remove ion', remove_ions)
    calculated_solvent = ion_imbalance + target_water + (2*target_ions) 
    #calculated_solvent = ion_imbalance + current_neutral_solvent - (2*remove_ions) - remove_water
    remaining_diff = calculated_solvent - Nmol
    #if there is any difference due to rounding, we add or subtract a single water molecule, it should have little bearing on the overall concentration 
    if remaining_diff != 0:
        remove_water_updated = remove_water + remaining_diff
    else:
        remove_water_updated = remove_water
    water_residues = [res for res in modeller.topology.residues() if res.name == 'HOH']
    na_residues = [res for res in modeller.topology.residues() if res.name == 'NA']
    cl_residues = [res for res in modeller.topology.residues() if res.name == 'CL']
    na_check = len(na_residues)

    to_remove_total = []
    if remove_water > 0:
        to_remove_water = np.random.choice(water_residues, remove_water_updated, replace=False)
        to_remove_total.extend(to_remove_water)
        #modeller.delete(to_remove_water)
    if remove_ions > 0:
        to_remove_na = np.random.choice(na_residues, remove_ions, replace=False)
        to_remove_total.extend(to_remove_na)
        to_remove_cl = np.random.choice(cl_residues, remove_ions, replace=False)
        to_remove_total.extend(to_remove_cl)
    #    modeller.delete(to_remove_cl)
    #to_remove_na = np.random.choice(na_residues, remove_ions, replace=False)
    #modeller.delete(to_remove_na)
    modeller.delete(to_remove_total)
    remove_na = len(to_remove_na)
    #to_remove_cl = np.random.choice(cl_residues, remove_ions, replace=False)
    #modeller.delete(to_remove_cl)
    final_water_count = sum([1 for res in modeller.topology.residues() if res.name == 'HOH'])
    print('final water', final_water_count)
    final_na_count = sum([1 for res in modeller.topology.residues() if res.name == 'NA'])
    print('final na', final_na_count)
    final_cl_count = sum([1 for res in modeller.topology.residues() if res.name == 'CL'])
    print( 'final cl', final_cl_count)
    final_solvent_count = final_water_count + final_na_count + final_cl_count
    print(final_solvent_count)
    with open(info_sheet, 'a') as file:
        file.write(f'original water total: {current_water}\n')
        file.write(f'original Na total: {current_na}\n')
        file.write(f'original Cl total: {current_cl}\n')
        file.write(f'original total solvent: {current_total_solvent}\n')
        file.write(f'ion imbalance: {ion_imbalance}\n')
        file.write(f'original neutral solvent: {current_neutral_solvent}\n')
        file.write(f'target neutral solvent: {target_neutral_solvent}\n')
        file.write(f'excess solvent: {excess_solvent}\n')
        file.write(f'water fraction: {water_fraction}\n')
        file.write(f'water molecules to remove: {remove_water}\n')
        file.write(f'water molecules to remove after correction: {remove_water_updated}\n')
        file.write(f'ion fraction: {ion_fraction}\n')
        file.write(f'ions to remove: {remove_ions}\n')
        file.write(f'Na check: {na_check}\n')
        file.write(f'final water total: {final_water_count}\n')
        file.write(f'remove Na: {remove_na}\n')
        file.write(f'final Cl count {final_cl_count}\n')
        file.write(f'final total solvent: {final_solvent_count}\n')
    # to check we got 50000 solvent molecules back (50,001 or 49999 due to rounding errors is fine)
   
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.PME, nonbondedCutoff=1.0*nanometer, constraints=app.HBonds)
    return system

# define simulation parameter and simulation.context, combining modeller and system, and setting integrator 
def setup_simulation(modeller, system):
    system.addForce(mm.MonteCarloBarostat(1*bar, 300*kelvin))
    integrator = mm.LangevinMiddleIntegrator(310*kelvin, 1/picosecond, 0.001*picoseconds)
    simulation = app.Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)
    return simulation, integrator

# combine setup_forcefield and setup_system functions, redefine box vectors ahead of local minimisation. Output final minimised topology and mimimisation, and return simulation after minimisation 
def energy_minimization(modeller, info_sheet, no_restraints: bool, salting=None):
    forcefield = setup_forcefield()
    solvmol = 50000
    padding = 1 #nanometer
    if salting:
        system = setup_system_proportional(modeller, forcefield, solvmol, padding, info_sheet, salting)
    else:
        system = setup_system(modeller, forcefield, solvmol, padding)
    xvec = mm.Vec3(11.70, 0.0, 0.0)
    yvec = mm.Vec3(0.0, 11.70, 0.0)
    zvec = mm.Vec3(0.0, 0.0, 11.70)
    if not no_restraints:
        apply_force_restraints(system, modeller)
        #apply_mass_restraints(system, modeller.topology)
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
    with open(info_sheet, 'a') as file:
        file.write(f'Initial potential energy of solution: {init_pe}\n')
        file.write(f'Minimised potential energy of solution: {final_pe}\n')
    return simulation.topology, final_state.getPositions(), final_pe, simulation, integrator, system

# takes minimised solvation box to run grand canonical ensemble (NVT), with initial temperature increase in 100-step intervals. 0,5ns of equilibration is allowed to achieve equilibrium state. After which the requested number of steps is run, with reporting corresponding to the number of steps in the --report_every option
def md_nvt(simulation, csvname: str, totalsteps: int, reprate: int, pdbname, integrator, system, initsteps: int, initpdb: str, relax_steps: int, relax_report_every: int, no_restraints: bool, test_conditions: bool):
    init_state = simulation.context.getState(getEnergy=True, getPositions=True)
    init_pe = init_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole)
    for temp in range(0, 310, 10):
        integrator.setTemperature(temp*kelvin)
        if test_conditions:
            simulation.reporters.append(app.StateDataReporter(stdout, 100, step=True, potentialEnergy=True, temperature=True, volume=True))
            simulation.step(100)
        else:
            simulation.reporters.append(app.StateDataReporter(stdout, 100, step=True, potentialEnergy=True, temperature=True, volume=True))
            simulation.step(10000)
    simulation.context.setVelocitiesToTemperature(310*kelvin)
    #simulation.step(496800)
    simulation.step(initsteps)
    if not no_restraints:
        remove_force_restraints(simulation)
        #remove_mass_restraints(system, simulation.topology)
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

# output rmsf per residue per chain based on generated PDB trajectory 
def rmsf_analysis_by_chain(pdb_traj: str, rmsfcsv: str):
    u = mda.Universe(pdb_traj)
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
    arguments = docopt(__doc__, version='openmmMD_dualmin_restraint_NPT.py')
    #fix PDB input file and impose protonation state to recreate chosen pH 
    pdb_clean = clean_pdb(arguments['--inpdb'], arguments['--pH'])
    pdb = pdb_clean[0]
    stem = pdb_clean[1]
    if arguments["--salting"]:
        pH_str = str("_" + arguments["--salting"] + "M_" + arguments['--pH'] + "pH")
    else:
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
        lowest_energy_value = energies[min_energy_index]
        app.PDBFile.writeFile(pdb.topology, lowest_energy_positions, open(pdbout, 'w'), keepIds=True)
        pdb = app.PDBFile(pdbout)
        info_sheet = str(stem + pH_str + "iter" + strj + "_info.txt")
        with open(info_sheet, 'w') as file:
            file.write(f'Minimum potential energy of protein-in-vacuo: {lowest_energy_value}\n')
        #set up simulation of protein solvated by water, conditional on minimisation creating net potential energy state (E_tot < 0) to avoid box explosion
        max_iterations = 100
        modeller = setup_modeller(pdb)
        second_iteration = 0
        while second_iteration < max_iterations:
            min_pdb = energy_minimization(modeller, info_sheet, arguments['--no_restraints'], arguments['--salting'])
            outputII = min_pdb[2]
            thresholdII = 0
            if outputII < thresholdII:
                break
            second_iteration += 1
        min_out = str(stem + pH_str + "iter" + strj + ".pdb")
        app.PDBFile.writeFile(min_pdb[0], min_pdb[1], open(min_out, "w"), keepIds=True)
        # set system, simulation and integrator parameters to run main trajectory, define output CSV and trajectory PDB names, set total steps, report rate and initial steps (for equilibration), set as inputs to NVT simulation function (md_nvt)
        simulation = min_pdb[3]
        integrator = min_pdb[4]
        system = min_pdb[5]
        initpdb = str(stem + pH_str + "_structure_sampling_" + strj + ".pdb")
        csvname = str(stem + pH_str + "_traj" + strj + ".csv")
        pdbname = str(stem + pH_str + "_traj" + strj + ".pdb")
        steps =int(arguments['--steps'])
        report_every = int(arguments['--report_every'])
        init_steps = int(arguments['--init_steps']) # total unreported step during solvent equilibration/protein restrained
        if arguments['--test_conditions']:
            temp_steps = 3100 # total number of steps for heating from 0K to 310K
        else:
            temp_steps = 310000 # total number of steps for heating from 0K to 310K
        equil_steps = init_steps - temp_steps
        relax_steps = int(arguments["--relax_steps"])
        relax_report_every = int(arguments['--relax_report_every'])
        sim_run = md_nvt(simulation, csvname, steps, report_every, pdbname, integrator, system, equil_steps, initpdb, relax_steps, relax_report_every, arguments['--no_restraints'], arguments['--test_conditions'])
        av_energy = sim_run[2]
        with open(info_sheet, 'a') as file:
            file.write(f'Average trajectory potential energy: {av_energy}\n')
        # reference (final) and feedstock (initial) frames in case of future use
        sim_ref = str(stem + pH_str + "_reference" + strj + ".pdb")
        app.PDBFile.writeFile(sim_run[4], sim_run[5], open(sim_ref, "w"), keepIds=True)
        sim_fin = str(stem + pH_str + "_feedstock" + strj + ".pdb")
        app.PDBFile.writeFile(sim_run[4], sim_run[6], open(sim_fin, "w"), keepIds=True)
        # set rmsf output name and perform rmsf analysis on output trajectory
        rmsfstem = str(stem + pH_str + "_traj" + strj)
        rmsf_analysis_by_chain(pdbname, rmsfstem)
     

if __name__ == '__main__':
    arguments = docopt(__doc__, version='openmmMD_dualmin_restraint_NPT.py')
    logging.getLogger().setLevel(logging.INFO)  
    main()