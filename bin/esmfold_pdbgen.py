#!/usr/bin/env python3
"""esmfold_pdbgen

Converting FASTA file sequences to PDB files via use of ESMFold
Usage:
esmfold_pdbgen.py [--i=<fasta>]

Options:
--i=<fasta>    Input fasta file containing protein sequence
"""
import logging
from docopt import docopt
import torch
from esm import pretrained, FastaBatchedDataset

# Load the ESMFold model
model, alphabet = pretrained.esmfold_v0()
model = model.eval().cuda() if torch.cuda.is_available() else model.eval()

# Example FASTA sequence
sequence = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRIA"

# Tokenize sequence and prepare dataset
batch_converter = alphabet.get_batch_converter()
data = [("protein1", sequence)]
batch_labels, batch_strs, batch_tokens = batch_converter(data)

with torch.no_grad():
    results = model.infer_pdb(batch_tokens)

# Write the PDB structure to a file
with open("output.pdb", "w") as pdbfile:
    pdbfile.write(results)

print("PDB structure saved to output.pdb")
