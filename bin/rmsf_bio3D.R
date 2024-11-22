#!/usr/bin/env Rscript
library(docopt)
#library(ggplot2)
#library(dplyr)
library(bio3d)
doc <- "Usage:
rmsf_bio3D.R --inpdb=<pdb> 
Options:
--inpdb=<pdb>  load input pdb trajectory
"
opts <- docopt(doc)
traj <- opts$`--inpdb`
print("Loading trajectory...")
pdb <- read.pdb(traj, multi=TRUE)
print("Trajectory loaded")
#xyz <- pdb$xyz
stem <- sub(".pdb", "", traj)
csvout <- paste0(stem, "_rmsftest.csv")
#na_resno_indices <- which(is.na(pdb$atom$resno))
#na_eleno_indices <- which(is.na(pdb$atom$eleno))
# Identify N/A residues
na_residues <- pdb$atom[is.na(pdb$atom$resno), ]

# Write N/A residues to a file
#na_file <- sub(".pdb", "_NA_residues.csv", traj)
#write.csv(na_residues, file = na_file, row.names = FALSE)

# Identify N/A residues
#na_elements <- pdb$atom[is.na(pdb$atom$eleno), ]

# Write N/A residues to a file
#naelements_file <- sub(".pdb", "_NA_elements.csv", traj)
#write.csv(na_elements, file = naelements_file, row.names = FALSE)
# align trajectory
print("Starting alignment")
#aligned <- fit.xyz(xyz[1, ], xyz)
#print("Alignment finished")
# select mainchain atoms
#sele <- atom.select(pdb, elety=c("CA"))

# Align the trajectory (select CA atoms for alignment)

non_HOH_indices <- atom.select(pdb, resid = "HOH", inverse = TRUE)
pdb$atom <- pdb$atom[non_HOH_indices$atom, ]
pdb$xyz <- pdb$xyz[, non_HOH_indices$xyz]

valid_residues <- !is.na(pdb$atom$resno) & !duplicated(pdb$atom$eleno)
duplicated_indices <- duplicated(pdb$atom$eleno)

pdb$atom <- pdb$atom[!duplicated_indices, ]
pdb$xyz <- pdb$xyz[, !duplicated_indices]
#inds <- atom.select(pdb, resno = pdb$atom$resno[valid_residues])
#ca_indices <- atom.select(inds, elety = "CA")
ca_indices <- atom.select(pdb, resno = pdb$atom[valid_residues, ]$resno, elety = "CA", "protein")
#average_structure <- apply(pdb$xyz[, ca_indices$xyz], 2, mean, na.rm = TRUE)
aligned <- fit.xyz(fixed=pdb$xyz[2000, ], mobile=pdb$xyz, fixed.inds = ca_indices$xyz, mobile.inds = ca_indices$xyz)
#aligned <- fit.xyz(fixed=average_structure, mobile=pdb[, ca_indices$xyz], fixed.inds = ca_indices$xyz, mobile.inds = ca_indices$xyz)
#aligned <- fit.xyz(fixed = average_structure, mobile = pdb$xyz, fixed.inds = ca_indices$xyz, mobile.inds = ca_indices$xyz)
print("Alignment finished")
# residue numbers to group by
#resno <- pdb$atom$resno[sele$atom]
# mean rmsf value of mainchain atoms of each residue
print("RMSF calculation")
rmsf_values <- rmsf(aligned[, ca_indices$xyz], grpby = pdb$atom$resno[ca_indices$atom])
results <- data.frame(Residue = pdb$atom$resno[ca_indices$atom], RMSF = rmsf_values)
write.table(results, file = csvout, sep = ",", row.names = FALSE, col.names = c("Residue", "RMSF"))
