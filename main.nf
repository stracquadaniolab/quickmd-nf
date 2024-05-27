// enabling nextflow DSL v2
nextflow.enable.dsl=2

process OpenmmMD {
    publishDir "${params.resultsDir}/openmm-minimise/", pattern: "*_fixed.pdb", mode: 'copy'
    publishDir "${params.resultsDir}/openmm-minimise/", pattern: "*_reformat.csv", mode: 'copy'
    publishDir "${params.resultsDir}/openmm-minimise/", pattern: "*_minimised.pdb", mode: 'copy'
    publishDir "${params.resultsDir}/openmm-minimise/", pattern: "data*", mode: 'copy'
    input:
        path incsv
        path inpdb
    output:
        path "*_fixed.pdb", emit: fixed_pdbs
        path "*_reformat.csv", emit: csv_reformat
        path '*_minimsed.pdb', emit: minimised_pdbs
        path 'data*', emit: data
    
    shell:
    """
    openmmMD.py --incsv $incsv --from-col ${params.csv.col} --in-pdb $inpdb ${params.mutant.maker.args}
    """
}
workflow {
    channel.fromPath("${params.inputFile}") | PrintHelloOnFile
}

