// enabling nextflow DSL v2

nextflow.enable.dsl=2

// process to identify mutations from a CSV file, then emits product variant PDB files from original wildtype PDB file
process FindMutations {
    publishDir "${params.resultsDir}-${params.pdb.maker.args}pH_${new Date().format('yyyy-MM-dd')}/openmm-md/", pattern: "*_centered.pdb", mode: 'copy'
    input:
        path wtin
        path varcsv
    output:
        path "*_centered.pdb", emit: mutated_pdbs
    
    shell:
    """
    mutant_creator.py --wtin $wtin --varlist $varcsv ${params.find_mut.maker.args}
    """
    errorStrategy 'ignore'
}

// original process to setup and execute NVT simulations via OpenMM featuring single mutation variant creation (no step control/report rate control/deparacated)
process OpenmmMD {

    // results/2024-11-10/<pdb>/<ph>/* 
    publishDir "${params.resultsDir}-${params.pdb.maker.args}pH_${new Date().format('yyyy-MM-dd')}/openmm-md/", pattern: "*_fixed.pdb", mode: 'copy'
    publishDir "${params.resultsDir}-${params.pdb.maker.args}pH_${new Date().format('yyyy-MM-dd')}/openmm-md/", pattern: "*_reformat.csv", mode: 'copy'
    publishDir "${params.resultsDir}-${params.pdb.maker.args}pH_${new Date().format('yyyy-MM-dd')}/openmm-md/", pattern: "*_minimised.pdb", mode: 'copy'
    publishDir "${params.resultsDir}-${params.pdb.maker.args}pH_${new Date().format('yyyy-MM-dd')}/openmm-md/", pattern: "*traj*", mode: 'copy'
    publishDir "${params.resultsDir}-${params.pdb.maker.args}pH_${new Date().format('yyyy-MM-dd')}/openmm-md/", pattern: "*rmsf*.csv", mode: 'copy'
    input:
        path incsv
        path inpdb
    output:
        path "*_fixed.pdb", emit: fixed_pdbs
        path "*_reformat.csv", emit: csv_reformat
        path '*traj*.pdb', emit: traj_pdb
        path '*traj*.csv', emit: traj_csv
        path '*rmsf*.csv', emit: rmsf_csv
        path 'wt_rmsf1.csv', emit: rmsf_wt
    
    shell:
    """
    openmmMD.py --incsv $incsv --from-col ${params.csv.col} --in-pdb $inpdb ${params.mutant.maker.args}
    """
}


// process to setup and execute NVT simulations via OpenMM, for a single pH/protonation state with set steps/reporting rate for protein solvated with 50,000 water molecules, approximating density at 1atm. PDB trajectory outputted, along with RMSF analysis performed with MDAnalysis 
process OpenmmMDNoCreate {
    publishDir "${params.resultsDir}-${params.pdb.maker.args}pH_${new Date().format('yyyy-MM-dd')}/openmm-md/", pattern: "*_fixed.pdb", mode: 'copy'
    publishDir "${params.resultsDir}-${params.pdb.maker.args}pH_${new Date().format('yyyy-MM-dd')}/openmm-md/", pattern: "*vac_min.pdb", mode: 'copy'
    publishDir "${params.resultsDir}-${params.pdb.maker.args}pH_${new Date().format('yyyy-MM-dd')}/openmm-md/", pattern: "*_minimised.pdb", mode: 'copy'
    publishDir "${params.resultsDir}-${params.pdb.maker.args}pH_${new Date().format('yyyy-MM-dd')}/openmm-md/", pattern: "*traj*", mode: 'copy'
    publishDir "${params.resultsDir}-${params.pdb.maker.args}pH_${new Date().format('yyyy-MM-dd')}/openmm-md/", pattern: "*initial*", mode: 'copy'
    publishDir "${params.resultsDir}-${params.pdb.maker.args}pH_${new Date().format('yyyy-MM-dd')}/openmm-md/", pattern: "*rmsf*.csv", mode: 'copy'
    publishDir "${params.resultsDir}-${params.pdb.maker.args}pH_${new Date().format('yyyy-MM-dd')}/openmm-md/", pattern: "*.txt", mode: 'copy'
    input:
        path inpdb
    output:
        path "*_fixed.pdb", emit: fixed_pdbs
        path '*traj*.pdb', emit: traj_pdb
        path '*traj*.csv', emit: traj_csv
        path '*rmsf*.csv', emit: rmsf_csv
    
    shell:
    """
    openmmMD_dualmin.py  --inpdb $inpdb --pH ${params.pdb.maker.args}
    """
}

// process accepting variant/pH tuple to setup and execute NVT simulations via OpenMM, for multiple pHs/protonation states with  steps/reporting rate definable in command line arguments. Sampling of protein-in-vacuo energetic state to find minimum energy, with conditional solution optimisation to ensure simulation stability. Protein solvated with 50,000 water molecules, approximating density at 1atm. PDB trajectory outputted, along with chain-by-chain RMSF analysis performed with MDAnalysis
process OpenmmMDTuple {
    publishDir "${params.resultsDir}_${params.total_length}_${new Date().format('yyyy-MM-dd')}/openmm-md/", pattern: "*_fixed.pdb", mode: 'copy'
    publishDir "${params.resultsDir}_${params.total_length}_${new Date().format('yyyy-MM-dd')}/openmm-md/", pattern: "*vac_min.pdb", mode: 'copy'
    publishDir "${params.resultsDir}_${params.total_length}_${new Date().format('yyyy-MM-dd')}/openmm-md/", pattern: "*_minimised.pdb", mode: 'copy'
    publishDir "${params.resultsDir}_${params.total_length}_${new Date().format('yyyy-MM-dd')}/openmm-md/", pattern: "*traj*", mode: 'copy'
    publishDir "${params.resultsDir}_${params.total_length}_${new Date().format('yyyy-MM-dd')}/openmm-md/", pattern: "*initial*", mode: 'copy'
    publishDir "${params.resultsDir}_${params.total_length}_${new Date().format('yyyy-MM-dd')}/openmm-md/", pattern: "*rmsf*.csv", mode: 'copy'
    publishDir "${params.resultsDir}_${params.total_length}_${new Date().format('yyyy-MM-dd')}/openmm-md/", pattern: "*.txt", mode: 'copy'
    input:
        tuple path(file), val(value)
    output:
        path "*_fixed.pdb", emit: fixed_pdbs
        path '*traj*.pdb', emit: traj_pdb
        path '*traj*.csv', emit: traj_csv
        path '*rmsf_chA*', emit: rmsf_chA
        path '*rmsf_chB*', emit: rmsf_chB
    
    shell:
    """
    openmmMD_dualmin_restraint.py  --inpdb $file --pH $value ${params.pdb.maker.args}
    """

    errorStrategy 'ignore'

    stub:
    """
    filename=\$(basename "$file")
    stem="\${filename%.*}"
    touch "${stem}_${value}_fixed.pdb"
    touch "${stem}_${value}_vac_min.pdb"
    touch "${stem}_${value}_minimised.pdb"
    touch "${stem}_${value}_initial1.pdb"
    touch "${stem}_${value}_info_sheet.txt"
    touch "${stem}_${value}_traj1.pdb"
    touch "${stem}_${value}_traj1.csv"
    touch "${stem}_${value}_rmsf_chA.csv"
    touch "${stem}_${value}_rmsf_csv.csv"
    """
    errorStrategy 'ignore'
}


// process accepting variant/pH/salting tuple to setup and execute NVT simulations via OpenMM, for multiple pHs/protonation states with  steps/reporting rate definable in command line arguments. Sampling of protein-in-vacuo energetic state to find minimum energy, with conditional solution optimisation to ensure simulation stability. Protein solvated with 50,000 water molecules, approximating density at 1atm. PDB trajectory outputted, along with chain-by-chain RMSF analysis performed with MDAnalysis
process OpenmmMDTupleSalting {
    publishDir "${params.resultsDir}_${params.total_length}_${params.salting}_${new Date().format('yyyy-MM-dd')}/openmm-md/", pattern: "*_fixed.pdb", mode: 'copy'
    publishDir "${params.resultsDir}_${params.total_length}_${params.salting}_${new Date().format('yyyy-MM-dd')}/openmm-md/", pattern: "*vac_min.pdb", mode: 'copy'
    publishDir "${params.resultsDir}_${params.total_length}_${params.salting}_${new Date().format('yyyy-MM-dd')}/openmm-md/", pattern: "*_minimised.pdb", mode: 'copy'
    publishDir "${params.resultsDir}_${params.total_length}_${params.salting}_${new Date().format('yyyy-MM-dd')}/openmm-md/", pattern: "*traj*", mode: 'copy'
    publishDir "${params.resultsDir}_${params.total_length}_${params.salting}_${new Date().format('yyyy-MM-dd')}/openmm-md/", pattern: "*rmsf*.csv", mode: 'copy'
    publishDir "${params.resultsDir}_${params.total_length}_${params.salting}_${new Date().format('yyyy-MM-dd')}/openmm-md/", pattern: "*.txt", mode: 'copy'
    input:
        tuple path(file), val(value), val(salt)
    output:
        path "*_fixed.pdb", emit: fixed_pdbs
        path '*traj*.pdb', emit: traj_pdb
        path '*traj*.csv', emit: traj_csv
        path '*rmsf_chA*', emit: rmsf_chA
        path '*rmsf_chB*', emit: rmsf_chB
        path '*_minimised.pdb'
        path '*vac_min.pdb'
        path '*.txt'
    
    errorStrategy 'ignore'

    shell:
    """
    openmmMD_dualmin_restraint.py  --inpdb $file --pH $value ${params.pdb.maker.args} --salting $salt
    """

    stub:
    """
    filename=\$(basename "$file")
    stem="\${filename%.*}"
    touch "${stem}_${value}_fixed.pdb"
    touch "${stem}_${value}_vac_min.pdb"
    touch "${stem}_${value}_minimised.pdb"
    touch "${stem}_${value}_initial1.pdb"
    touch "${stem}_${value}_info_sheet.txt"
    touch "${stem}_${value}_traj1.pdb"
    touch "${stem}_${value}_traj1.csv"
    touch "${stem}_${value}_rmsf_chA.csv"
    touch "${stem}_${value}_rmsf_csv.csv"
    """
    errorStrategy 'ignore'
}

// process accepting variant/pH/salting tuple to setup and execute NVT simulations via OpenMM, for multiple pHs/protonation states with  steps/reporting rate definable in command line arguments. Sampling of protein-in-vacuo energetic state to find minimum energy, with conditional solution optimisation to ensure simulation stability. Protein solvated with 50,000 water molecules, approximating density at 1atm. PDB trajectory outputted, along with chain-by-chain RMSF analysis performed with MDAnalysis
process OpenmmMDTuplePHSaltingGroup {
    publishDir "${params.resultsDir}_${params.total_length}_${params.salting}_${new Date().format('yyyy-MM-dd')}/openmm-md/", pattern: "*_fixed.pdb", mode: 'copy'
    publishDir "${params.resultsDir}_${params.total_length}_${params.salting}_${new Date().format('yyyy-MM-dd')}/openmm-md/", pattern: "*vac_min.pdb", mode: 'copy'
    publishDir "${params.resultsDir}_${params.total_length}_${params.salting}_${new Date().format('yyyy-MM-dd')}/openmm-md/", pattern: "*_minimised.pdb", mode: 'copy'
    publishDir "${params.resultsDir}_${params.total_length}_${params.salting}_${new Date().format('yyyy-MM-dd')}/openmm-md/", pattern: "*traj*", mode: 'copy'
    publishDir "${params.resultsDir}_${params.total_length}_${params.salting}_${new Date().format('yyyy-MM-dd')}/openmm-md/", pattern: "*rmsf*.csv", mode: 'copy'
    publishDir "${params.resultsDir}_${params.total_length}_${params.salting}_${new Date().format('yyyy-MM-dd')}/openmm-md/", pattern: "*.txt", mode: 'copy'
    input:
        tuple path(file), val(group)
    output:
        path "*_fixed.pdb", emit: fixed_pdbs
        path '*traj*.pdb', emit: traj_pdb
        path '*traj*.csv', emit: traj_csv
        path '*rmsf_chA*', emit: rmsf_chA
        path '*rmsf_chB*', emit: rmsf_chB
        path '*_minimised.pdb'
        path '*vac_min.pdb'
        path '*.txt'
    
    // errorStrategy 'ignore'

    shell:
    """
    value = $group[0]
    salt = $group[1]
    openmmMD_dualmin_restraint.py  --inpdb $file --pH $value ${params.pdb.maker.args} --salting $salt
    """

    stub:
    """
    filename=\$(basename "$file")
    stem="\${filename%.*}"
    touch "${stem}_${value}_fixed.pdb"
    touch "${stem}_${value}_vac_min.pdb"
    touch "${stem}_${value}_minimised.pdb"
    touch "${stem}_${value}_initial1.pdb"
    touch "${stem}_${value}_info_sheet.txt"
    touch "${stem}_${value}_traj1.pdb"
    touch "${stem}_${value}_traj1.csv"
    touch "${stem}_${value}_rmsf_chA.csv"
    touch "${stem}_${value}_rmsf_csv.csv"
    """
    errorStrategy 'ignore'
}


// process accepting variant/pH/salting tuple to setup and execute NPT simulations via OpenMM, for multiple pHs/protonation states with  steps/reporting rate definable in command line arguments. Sampling of protein-in-vacuo energetic state to find minimum energy, with conditional solution optimisation to ensure simulation stability. Protein solvated with 50,000 water molecules, approximating density at 1atm. PDB trajectory outputted, along with chain-by-chain RMSF analysis performed with MDAnalysis
process OpenmmMDTupleSaltingNPT {
    publishDir "${params.resultsDir}_${params.total_length}_${params.salting}_${new Date().format('yyyy-MM-dd')}/openmm-md/", pattern: "*_fixed.pdb", mode: 'copy'
    publishDir "${params.resultsDir}_${params.total_length}_${params.salting}_${new Date().format('yyyy-MM-dd')}/openmm-md/", pattern: "*vac_min.pdb", mode: 'copy'
    publishDir "${params.resultsDir}_${params.total_length}_${params.salting}_${new Date().format('yyyy-MM-dd')}/openmm-md/", pattern: "*_minimised.pdb", mode: 'copy'
    publishDir "${params.resultsDir}_${params.total_length}_${params.salting}_${new Date().format('yyyy-MM-dd')}/openmm-md/", pattern: "*traj*", mode: 'copy'
    publishDir "${params.resultsDir}_${params.total_length}_${params.salting}_${new Date().format('yyyy-MM-dd')}/openmm-md/", pattern: "*rmsf*.csv", mode: 'copy'
    publishDir "${params.resultsDir}_${params.total_length}_${params.salting}_${new Date().format('yyyy-MM-dd')}/openmm-md/", pattern: "*.txt", mode: 'copy'
    input:
        tuple path(file), val(value), val(salt)
    output:
        path "*_fixed.pdb", emit: fixed_pdbs
        path '*traj*.pdb', emit: traj_pdb
        path '*traj*.csv', emit: traj_csv
        path '*rmsf_chA*', emit: rmsf_chA
        path '*rmsf_chB*', emit: rmsf_chB
        path '*_minimised.pdb'
        path '*vac_min.pdb'
        path '*.txt'
    
    // errorStrategy 'ignore'

    shell:
    """
    openmmMD_dualmin_restraint_NPT.py  --inpdb $file --pH $value ${params.pdb.maker.args} --salting $salt
    """

    stub:
    """
    filename=\$(basename "$file")
    stem="\${filename%.*}"
    touch "${stem}_${value}_fixed.pdb"
    touch "${stem}_${value}_vac_min.pdb"
    touch "${stem}_${value}_minimised.pdb"
    touch "${stem}_${value}_initial1.pdb"
    touch "${stem}_${value}_info_sheet.txt"
    touch "${stem}_${value}_traj1.pdb"
    touch "${stem}_${value}_traj1.csv"
    touch "${stem}_${value}_rmsf_chA.csv"
    touch "${stem}_${value}_rmsf_csv.csv"
    """
    errorStrategy 'ignore'
}


// deprecated wildtype-only version of OpenmmMDNoCreate
process OpenmmMDNoCreateWT {
    publishDir "${params.resultsDir}-${params.pdb.maker.args}pH_${new Date().format('yyyy-MM-dd')}/openmm-md/", pattern: "*_fixed.pdb", mode: 'copy'
    publishDir "${params.resultsDir}-${params.pdb.maker.args}pH_${new Date().format('yyyy-MM-dd')}/openmm-md/", pattern: "*_minimised.pdb", mode: 'copy'
    publishDir "${params.resultsDir}-${params.pdb.maker.args}pH_${new Date().format('yyyy-MM-dd')}/openmm-md/", pattern: "*traj*", mode: 'copy'
    publishDir "${params.resultsDir}-${params.pdb.maker.args}pH_${new Date().format('yyyy-MM-dd')}/openmm-md/", pattern: "*initial*", mode: 'copy'
    publishDir "${params.resultsDir}-${params.pdb.maker.args}pH_${new Date().format('yyyy-MM-dd')}/openmm-md/", pattern: "*rmsf*.csv", mode: 'copy'
    publishDir "${params.resultsDir}-${params.pdb.maker.args}pH_${new Date().format('yyyy-MM-dd')}/openmm-md/", pattern: "*.txt", mode: 'copy'
    input:
        path inpdb
    output:
        path "*_fixed.pdb", emit: fixed_pdbs
        path '*traj*.pdb', emit: traj_pdb
        path '*traj*.csv', emit: traj_csv
        path 'wt_rmsf*', emit: rmsf_wt
        //optional file 'wt_rmsf1.csv', emit: rmsf_wt
    
    shell:
    """
    openmmMD_twophase.py  --inpdb $inpdb --pH ${params.pdb.maker.args}
    """
}


// process measuring the average interchain distances between residues in the respective A and B chains to predict strength of hydrogen bonding in the protein secondary structure or drift caused by polar solvent
process InterchainPairs {
    publishDir "${params.resultsDir}_${params.total_length}_${new Date().format('yyyy-MM-dd')}/interchain_pairs_all_pH/", pattern: "*_interchain_pairs.pdb", mode: 'copy'
    publishDir "${params.resultsDir}_${params.total_length}_${new Date().format('yyyy-MM-dd')}/interchain_pairs_all_pH/", pattern: "*_average.pdb", mode: 'copy'
    input:
        path inpdb 
    output:
        path "*_interchain_pairs.pdb"
    
    shell:
    """
    interchain_pairs.py --inpdb $inpdb
    """
    stub:
    """
    filename=\$(basename "$inpdb")
    stem="\${filename%.*}"
    touch "${stem}_interchain_pairs.pdb"
    touch "${stem}_average.pdb"
    """
    errorStrategy 'ignore'
}

// process to output RMSF graphics with multiple RMSFs (deprecated)
process PlotRMSF {
    publishDir "${params.resultsDir}-${params.pdb.maker.args}_${new Date().format('yyyy-MM-dd')}/plot_RMSF/", pattern: "*.png", mode: 'copy'
    input:
        path incsv
    output:
        path "*.png"
        shell:
        """
        plot_rmsf.R --incsv ${incsv.join(',')}
        """
}

// process to plot RMSF figures comparing variant to wildtype RMSF (suited to single processes only)
process PlotRMSFindividual {
    publishDir "${params.resultsDir}_${new Date().format('yyyy-MM-dd')}/plot_rmsf/", pattern: "*.pdf", mode: 'copy'
    input:
        path incsv
        path wtcsv
    output:
        path "*.pdf"
        shell:
        """
        plot_rmsf_individual.R --incsv $incsv --wtcsv $wtcsv
        """
}


// process to compare variant RMSFs to wildtype, designed to deal with full set of data for given chain/pH
process PlotRMSFindividualTuple {
    publishDir "${params.resultsDir}_${params.total_length}_${new Date().format('yyyy-MM-dd')}/plot_RMSF/", pattern: "*.pdf", mode: 'copy'
    input:
        tuple path incsv, path wtcsv
    output:
        path "*.pdf"
        shell:
        """
        plot_rmsf_individual.R --incsv $incsv --wtcsv $wtcsv
        """
        stub:
        """
        filename=\$(basename "$incsv")
        stem="\${filename%.*}"
        touch "${stem}.pdf"
        """
        errorStrategy 'ignore'
}


// process RMSFs from PDB trajectory via Bio3D
process Rmsf_Bio3d {
    publishDir "${params.resultsDir}_${new Date().format('yyyy-MM-dd')}/plot_RMSF/", pattern: "*.csv", mode: 'copy'
    input:
        path inpdb
    output:
        path "*.csv"
    shell:
    """
    rmsf_bio3D.R --inpdb $inpdb
    """

}


// calculate RMSF from PDB trajecotries via MDAnalysis separately from OpenmmMDTuple process
process Rmsf_mda {
    publishDir "${params.resultsDir}_${new Date().format('yyyy-MM-dd')}/plot_RMSF/", pattern: "*.csv", mode: 'copy'
    input:
        path inpdb
    output:
        path "*.csv"
    shell:
    """
    rmsf_calc.py --inpdb $inpdb
    """

}

// workflow creating single mutations to generate single mutations from CSV file upon wildtype PDB, then generate PDB trajectories and RMSF graphics (not currently in use)
workflow MutantMakerFlow {
    inpath_ch = channel.fromPath("${params.inputFile}")
    incsv_ch = channel.fromPath("${params.inputCsv}")
    OpenmmMD(incsv_ch, inpath_ch)
    OpenmmMD.out.rmsf_csv
        .flatten()
        .set { individual_files }
    PlotRMSFindividual(individual_files, OpenmmMD.out.rmsf_wt)
    //PlotRMSF(OpenmmMD.out.rmsf_csv)
}

// workflow taking input PDB files (wildtype + variants), no mutation creation stage. Output RMSFs passed to PlotIndividual process for graphic creation
workflow IndividualPDBFlow {
    inpath_ch = channel.fromPath("${params.inputFile}") 
    OpenmmMDNoCreate(inpath_ch)
    OpenmmMD(incsv_ch, inpath_ch)
    PlotRMSFindividual(OpenmmMDNoCreate.out.rmsf_csv)
}

// workflow to input variants + wildtype PDBs to NVT simulation, with RMSFs filtered to separate wildtype from variant RMSFs for graphic creation (not currently in use)
workflow CollectRMSFsFlow {
    //inpath_ch = channel.fromPath("${params.inputPDB}")
    inpath_ch = channel.fromPath("${params.inputCEN}")
    //inpath_ch_wt = channel.fromPath("${params.inputWT}")
    OpenmmMDNoCreate(inpath_ch)
    //OpenmmMDNoCreateWT(inpath_ch_wt)
    //OpenmmMDNoCreate.out.rmsf_csv.filter { it.name == 'wt_rmsf1.csv' } into target_ch
    //PlotRMSFindividual(OpenmmMDNoCreate.out.rmsf_csv, OpenmmMDNoCreateWT.out.rmsf_wt)
}

// workflow to generate variant PDBs from either sequence libraries or mutation lists and pass to input of OpenmmMDTuple, in which all variants (+ wildtype) are simulated at all selected protonation states to generate output PDB trajectory and RMSF CSV file. Interchain distances are then tested on the output PDB trajectories
workflow MultipleCrossFlow {
    //Uncomment below alternative inpath_ch if using CSV to generate variants
    //wtin_ch = channel.fromPath("${params.inputFile}")
    //incsv_ch = channel.fromPath("${params.inputCsv}")
    //FindMutations(wtin_ch, incsv_ch)
    //inpath_ch = FindMutations.out.mutated_pdbs
    //FindMutations.out.mutated_pdbs
    //    .flatten()
    //    .set { inpath_ch }
    inpath_ch = channel.fromPath("${params.inputCEN}")
    //inpath_ch = channel.fromPath("${params.inputCENSPT}")
    //input_pH_ch = channel.of(4.6, 7.4)
    //input_salting_ch = channel.of(0.0935, 0.1817)
    //input_pH_salting_ch = input_pH_ch.cross(input_salting_ch)
    input_pH_salting_ch = channel.of([4.6, 0.0935], [7.4, 0.1817])
    //pdb_pH_combinations = inpath_ch.combine(input_pH_ch)
    //pdb_pH_salting_combinations = pdb_pH_combinations.combine(input_salting_ch)
    pdb_pH_salting_combinations = inpath_ch.combine(input_pH_salting_ch)
    //OpenmmMDTuple(pdb_pH_combinations)
    //OpenmmMDTupleSalting(pdb_pH_salting_combinations)
    OpenmmMDTupleSaltingNPT(pdb_pH_salting_combinations)
    //OpenmmMDTuplePHSaltingGroup(pdb_pH_salting_combinations)
    //OpenmmMDTuple.out.traj_pdb
    //    .flatten()
    //    .set { individual_files }
    
    //InterchainPairs(individual_files)
    
}


workflow RmsfFigureProcessing {
    def rmsf_chA = Channel.fromPath(params.inputrmsfA)
    // Chain A RMSFs processing
    def allFiles7pHchA = rmsf_chA.filter { it.name !~ /4.0pH/ }
    def allFiles4pHchA = rmsf_chA.filter { it.name !~ /7.4pH/ }
    def otherFiles7pHchA = allFiles7pHchA.filter { it.name !~ /wildtype/ }
    def otherFiles4pHchA = allFiles4pHchA.filter { it.name !~ /wildtype/ }
    def wtFile7pHchA = allFiles7pHchA.filter { it.name =~ /wildtype/ }.first()
    def wtFile4pHchA = allFiles4pHchA.filter { it.name =~ /wildtype/ }.first()
    
    otherFiles7pHchA.map { singleFile -> tuple(singleFile, wtFile7pHchA) }
        .set { filePairs7pHchA }
    otherFiles4pHchA.map { singleFile -> tuple(singleFile, wtFile4pHchA) }
        .set { filePairs4pHchA }
    
    PlotRMSFindividualTuple(filePairs4pHchA)
    PlotRMSFindividualTuple(filePairs7pHchA)

}


// workflow to calculate RMSFs independently
workflow Rmsf_calc {
    inpath_ch = channel.fromPath("${params.inputtraj}")
    //InterchainPairs(inpath_ch)
    //Rmsf_Bio3d(inpath_ch)
    Rmsf_mda(inpath_ch)
}

workflow {
//    MutantMakerFlow()
//    Rmsf_calc()
//    CollectRMSFsFlow()
    MultipleCrossFlow()
}