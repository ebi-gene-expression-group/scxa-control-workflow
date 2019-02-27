#!/usr/bin/env nextflow

// Generate a config file for the main sub-workflows

exp_name = params.exp_name
sdrf_dir = params.sdrf_dir

SDRF = Channel.fromPath("${sdrf_dir}/${exp_name}.sdrf.txt", checkIfExists: true)
IDF = Channel.fromPath("${sdrf_dir}/${exp_name}.idf.txt", checkIfExists: true)

// Derive the config files

process generate_config {

    publishDir "$SCXA_CONF/study", mode: 'copy', overwrite: true

    conda 'r-optparse r-data.table r-workflowscriptscommon'

    input:
        file(sdrf_file) from SDRF
        file(idf_file) from IDF

    output:
        file '*.conf' into CONFIG_FILES
        file '*.sdrf.txt' into SDRF_FILES
        
    """
    sdrfToNfConf.R \
        --sdrf=$sdrf_file \
        --idf=$idf_file \
        --name=$exp_name \
        --verbose \
    """
}

// Mark config files with species

process mark_conf_species {
    
    input:
        file(confFile) from CONFIG_FILES

    output:
        set stdout, file (confFile) into CONF_BY_SPECIES

    """
    echo $confFile | awk -F'.' '{print \$2}'
    """
}

process mark_sdrf_species {
    
    input:
        file(sdrfFile) from SDRF_FILES

    output:
        set stdout, file (sdrfFile) into SDRF_BY_SPECIES

    """
    echo $sdrfFile | awk -F'.' '{print \$2}'
    """
}

CONF_BY_SPECIES
    .join(SDRF_BY_SPECIES)
    .set{COMBINED_CONFIG}

// Run quantification

process quantify {

    conda 'nextflow'

    storeDir "$SCXA_RESULTS/$exp_name/$species/quantification"
    
    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus == 130  ? 'retry' : 'finish' }
    maxRetries 20
    
    input:
        set val(species), file (confFile), file(sdrfFile) from COMBINED_CONFIG

    output:
        set val(species), file ("*/*.abundance.h5") into QUANT_FILES 
        file('quantifcation.log')    

    """
        grep "sc_protocol" $confFile | grep "smart-seq" > /dev/null
        if [ \$? -eq 0 ]; then
            quantification_workflow=scxa-smartseq-workflow
        else
            echo "No workflow avialable for this experiment type" 1>&2
            exit 1
        fi

        RESULTS_ROOT=\$PWD
        SUBDIR=$exp_name/$species/\$quantification_workflow     

        pushd \$SCXA_WORK/\$SUBDIR > /dev/null

        nextflow run \
            -config \$RESULTS_ROOT/$confFile \
            --sdrf \$RESULTS_ROOT/$sdrfFile \
            --resultsRoot \$RESULTS_ROOT \
            -resume \
            \$quantification_workflow \
            -work-dir $SCXA_WORK/\$SUBDIR \
            -with-report $SCXA_RESULTS/reports/\$SUBDIR/report.html \
            -N $SCXA_REPORT_EMAIL \
            -with-dag $SCXA_RESULTS/reports/\$SUBDIR/flowchart.pdf

        if [ \$? -ne 0 ]; then
            echo "Workflow failed for $exp_name - $species - \$quantification_workflow" 1>&2
            exit 1
        fi
        
        popd > /dev/null

        cp $SCXA_WORK/\$SUBDIR/.nextflow.log quantification.log
   """
}


// Run aggregation

 
