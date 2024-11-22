#!/usr/bin/env python3
"""mutant_creator

Create variants to the wildtype PDB based on a CSV file containing a list of named variants with corresponding list of mutations for each one
Usage:
mutant_creator.py [--wtin=<wtin>] [--varlist=<varlist>] [--pH=<pH>]

Options:
--wtin=<wtin>          Wildype PDB file to be mutated
--varlist=<varlist>    Input CSV file containing variants with mutation lists
--pH=<pH>              Set pH to desired
"""
import logging
from docopt import docopt
from pdbfixer import PDBFixer
from openmm.app import PDBFile
import csv
from Bio.SeqUtils import seq3
from Bio import pairwise2
from Bio.pairwise2 import format_alignment


#function find_mutations for use when sequence data available only
def find_mutations(ref_seq, mut_seq):
    alignments = pairwise2.align.globalxx(ref_seq, mut_seq)
    aligned_ref, aligned_mut = alignments[0][0], alignments[0][1]
    mutations = []
    for i, (r, m) in enumerate(zip(aligned_ref, aligned_mut)):
        if r != m:
            mutation = f"{r}{i+1}{m}"
            mutations.append(mutation)
    return mutations

#convert mutations to a version readable by PDBFixer 
def convert_mutation_to_pdbfixer(mutation):
    original_aa = mutation[0]
    position = mutation[1:-1]
    new_aa = mutation[-1]

    # Convert to three-letter codes
    original_aa_3letter = seq3(original_aa).upper()
    new_aa_3letter = seq3(new_aa).upper()

    # Format for PDBFixer, e.g., 'GLU-48-ASP'
    pdbfixer_format = f"{original_aa_3letter}-{position}-{new_aa_3letter}"
    return pdbfixer_format

#Make aure WT is in the correct format
def clean_wildtype(pdbname: str, pH: str, pdbout: str):
    pH_fl = float(pH)
    pdb = PDBFixer(pdbname)
    #numChains = len(list(pdb.topology.chains()))
    #pdb.removeChains(range(1, numChains))
    pdb.findMissingResidues()
    #pdb.missingResidues = {}
    pdb.findNonstandardResidues()
    pdb.replaceNonstandardResidues()
    pdb.removeHeterogens(False)
    pdb.findMissingAtoms()
    pdb.addMissingAtoms()
    pdb.addMissingHydrogens(pH_fl)
    #PDBFile.writeFile(pdb.topology, pdb.positions, open("wildtype_fixed.pdb", 'w'), keepIds=True)
    PDBFile.writeFile(pdb.topology, pdb.positions, open(pdbout, 'w'))
    return pdb

#create vairants, implamenting mutations across both chains
def create_mutants(pdbname: str, mutant: list, chain: list, pH: str, pdbout: str):
    pH_fl = float(pH)
    mutpdb = PDBFixer(pdbname)
    for ch_list in chain:
        for mut_list in mutant:
            mutpdb.applyMutations([mut_list], ch_list)
    #mutpdb.applyMutations([mutant], chain)
    #mutpdb.findMissingResidues()
    mutpdb.missingResidues = {}
    #mutpdb.findNonstandardResidues()
    #mutpdb.replaceNonstandardResidues()
    mutpdb.removeHeterogens(False)
    mutpdb.findMissingAtoms()
    mutpdb.addMissingAtoms()
    mutpdb.addMissingHydrogens(pH_fl)
    PDBFile.writeFile(mutpdb.topology, mutpdb.positions, open( pdbout, 'w'), keepIds=True)
    return mutpdb

def main():
    arguments = docopt(__doc__, version='mutant_creator.py')
    wt = "wildtype_centered.pdb"
    clean_wildtype(arguments['--wtin'], arguments['--pH'], wt)
    chain = ["A", "B"]
    #Uncomment below line in case of sequence data, where mutations must be identified
    #wt_sequence = "wt sequence here"
    with open(arguments['--varlist'], 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            #comment below line in case of sequence data, uncomment line below it
            stem = row['design_id']
            #stem = row['description'] 
            pdbout = str(stem + "_centered.pdb")
            #Uncomment following two lines and comment out two lines below them if sequence data only 
            #mut_sequence = row['sequence']
            #filtered_list = find_mutations(wt_sequence, mut_sequence)
            mutation_list = row['mutations'].split(';')
            filtered_list = [mutation for mutation in mutation_list if mutation.lower() != 'c-terminal truncation']
            formatted_mutations = [convert_mutation_to_pdbfixer(mutation) for mutation in filtered_list]
            logging.info("%s %s", pdbout, ', '.join(formatted_mutations))
            create_mutants(wt, formatted_mutations, chain, arguments['--pH'], pdbout)

if __name__ == '__main__':
    arguments = docopt(__doc__, version='mutant_creator.py')
    logging.getLogger().setLevel(logging.INFO)  
    main()