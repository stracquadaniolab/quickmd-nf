// enabling nextflow DSL v2
nextflow.enable.dsl=2

process OpenmmMD {
    publishDir "${params.resultsDir}/openmm-md/", pattern: "*_fixed.pdb", mode: 'copy'
    publishDir "${params.resultsDir}/openmm-md/", pattern: "*_reformat.csv", mode: 'copy'
    publishDir "${params.resultsDir}/openmm-md/", pattern: "*_minimised.pdb", mode: 'copy'
    publishDir "${params.resultsDir}/openmm-md/", pattern: "*traj*", mode: 'copy'
    publishDir "${params.resultsDir}/openmm-md/", pattern: "*rmsf*.csv", mode: 'copy'
    input:
        path incsv
        path inpdb
    output:
        path "*_fixed.pdb", emit: fixed_pdbs
        path "*_reformat.csv", emit: csv_reformat
        path '*_minimsed.pdb', emit: minimised_pdbs
        path '*traj*.pdb', emit: traj_pdb
        path '*traj*.csv', emit: traj_csv
        path '*rmsf*.csv', emit: rmsf_csv
    
    shell:
    """
    openmmMD.py --incsv $incsv --from-col ${params.csv.col} --in-pdb $inpdb ${params.mutant.maker.args}
    """
}
workflow {
    inpath_ch = channel.fromPath("${params.inputFile}")
    incsv_ch = channel.fromPath("${params.inputCsv}")
    OpenmmMD(incsv_ch, inpath_ch)
}

