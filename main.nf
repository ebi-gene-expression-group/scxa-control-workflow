#!/usr/bin/env nextflow

// Generate a config file for the main sub-workflows

exp_name = params.exp_name
sdrf_dir = params.sdrf_dir

SDRF = Channel.fromPath("${sdrf_dir}/${exp_name}.sdrf.txt", checkIfExists: true)
IDF = Channel.fromPath("${sdrf_dir}/${exp_name}.idf.txt", checkIfExists: true)

// Derive the config files

process generate_config {

    conda 'r-data.table'

    input:
        file(sdrf_file) from SDRF
        file(idf_file) from IDF

    output:
        file '*/*.conf' into CONFIG_FILES
        
    """
    export WF_BINDIR=${workflow.projectDir}/bin

    sdrfToNfConf.R \
        --sdrf=$sdrf_file \
        --idf=$idf_file \
        --name=$exp_name \
        --verbose
    """
}

// Run quantification


// Run aggregation

 
