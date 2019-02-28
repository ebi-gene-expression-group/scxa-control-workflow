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
    echo $confFile | awk -F'.' '{printf "%s", \$2}'
    """
}

process mark_sdrf_species {
    
    input:
        file(sdrfFile) from SDRF_FILES

    output:
        set stdout, file (sdrfFile) into SDRF_BY_SPECIES

    """
    echo $sdrfFile | awk -F'.' '{printf "%s", \$2}'
    """
}

CONF_BY_SPECIES
    .join(SDRF_BY_SPECIES)
    .into{
        COMBINED_CONFIG_FOR_QUANTIFY
        COMBINED_CONFIG_FOR_AGGREGATION
    }

// Run quantification

process quantify {

    conda 'nextflow'

    storeDir "$SCXA_RESULTS/$exp_name/$species/quantification"
    
    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus == 130  ? 'retry' : 'finish' }
    maxRetries 20
    
    input:
        set val(species), file (confFile), file(sdrfFile) from COMBINED_CONFIG_FOR_QUANTIFY

    output:
        set val(species), file ("kallisto/*") into KALLISTO_DIRS 
        set val(species), file("reference/reference.fastq.gz") into REFERENCE_FASTA
        set val(species), file("reference/reference.gtf.gz") into REFERENCE_GTF
        set val(species), file('quantification.log')    

    """
        grep "sc_protocol" $confFile | grep "smart-seq" > /dev/null
        if [ \$? -eq 0 ]; then
            quantification_workflow=scxa-smartseq-quantification-workflow
        else
            echo "This is not a SMART-seq experiment" 1>&2
            exit 1
        fi

        RESULTS_ROOT=\$PWD
        SUBDIR="$exp_name/$species/\$quantification_workflow"     

        mkdir -p \$SCXA_WORK/\$SUBDIR
        mkdir -p $SCXA_RESULTS/\$SUBDIR/reports
        pushd \$SCXA_WORK/\$SUBDIR > /dev/null

        nextflow run \
            -config \$RESULTS_ROOT/$confFile \
            --sdrf \$RESULTS_ROOT/$sdrfFile \
            --resultsRoot \$RESULTS_ROOT \
            -resume \
            \$quantification_workflow \
            -work-dir $SCXA_WORK/\$SUBDIR \
            -with-report $SCXA_RESULTS/\$SUBDIR/reports/report.html \
            -N $SCXA_REPORT_EMAIL \
            -with-dag $SCXA_RESULTS/\$SUBDIR/reports/flowchart.pdf

        if [ \$? -ne 0 ]; then
            echo "Workflow failed for $exp_name - $species - \$quantification_workflow" 1>&2
            exit 1
        fi
        
        popd > /dev/null

        cp $SCXA_WORK/\$SUBDIR/.nextflow.log quantification.log
   """
}

// Run aggregation with https://github.com/ebi-gene-expression-group/scxa-aggregation-workflow

process aggregate {
    
    conda 'nextflow'

    storeDir "$SCXA_RESULTS/$exp_name/$species"
    
    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus == 130  ? 'retry' : 'finish' }
    maxRetries 20
    
    input:
        set val(species), file (confFile), file(sdrfFile) from COMBINED_CONFIG_FOR_AGGREGATION
        set val(species), file ('*') from KALLISTO_DIRS.collect()
        set val(species), file(referenceGtf) from REFERENCE_GTF

    output:
        set val(species), file("matrices/*_counts.zip") into KALLISTO_COUNT_MATRIX
        set val(species), file("matrices/*_tpm.zip") into KALLISTO_ABUNDANCE_MATRIX
        set val(species), file("matrices/*.stats.tsv") into KALLISTO_STATS
        set val(species), file('aggregation.log')    

    """
        RESULTS_ROOT=\$PWD
        SUBDIR="$exp_name/$species/aggregation"     

        mkdir -p \$SCXA_WORK/\$SUBDIR
        mkdir -p $SCXA_RESULTS/\$SUBDIR/reports
        pushd \$SCXA_WORK/\$SUBDIR > /dev/null

        nextflow run \
            -config \$RESULTS_ROOT/$confFile \
            --resultsRoot \$RESULTS_ROOT \
            --referenceGtf ${referenceGtf}\
            -resume \
            scxa-aggregation-workflow \
            -work-dir $SCXA_WORK/\$SUBDIR \
            -with-report $SCXA_RESULTS/\$SUBDIR/reports/report.html \
            -N $SCXA_REPORT_EMAIL \
            -with-dag $SCXA_RESULTS/\$SUBDIR/reports/flowchart.pdf

        if [ \$? -ne 0 ]; then
            echo "Workflow failed for $exp_name - $species - scxa_aggregation_workflow" 1>&2
            exit 1
        fi
        
        popd > /dev/null

        cp $SCXA_WORK/\$SUBDIR/.nextflow.log aggregation.log
   """
    
}
