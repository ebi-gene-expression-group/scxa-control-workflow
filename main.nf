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
        --verbose

    for f in *.conf; do
        echo -e "includeConfig '${baseDir}/params.config'\n"|cat - \$f > \$f.tmp && mv \$f.tmp \$f
    done
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
        COMBINED_CONFIG_FOR_SCANPY
    }

// Run quantification with https://github.com/ebi-gene-expression-group/scxa-smartseq-quantification-workflow

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
        SUBDIR="$exp_name/$species/quantification"     

        mkdir -p $SCXA_WORK/\$SUBDIR
        mkdir -p $SCXA_NEXTFLOW/\$SUBDIR
        mkdir -p $SCXA_RESULTS/\$SUBDIR/reports
        pushd $SCXA_NEXTFLOW/\$SUBDIR > /dev/null

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

        cp $SCXA_NEXTFLOW/\$SUBDIR/.nextflow.log quantification.log
   """
}

// Run aggregation with https://github.com/ebi-gene-expression-group/scxa-aggregation-workflow

process aggregate {
    
    conda 'nextflow'

    storeDir "$SCXA_RESULTS/$exp_name/$species/aggregation"
    
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
        set val(species), file('matrices/aggregation.log')    

    """
        RESULTS_ROOT=\$PWD
        SUBDIR="$exp_name/$species/aggregation"     

        mkdir -p $SCXA_WORK/\$SUBDIR
        mkdir -p $SCXA_NEXTFLOW/\$SUBDIR
        mkdir -p $SCXA_RESULTS/\$SUBDIR/reports
        pushd $SCXA_NEXTFLOW/\$SUBDIR > /dev/null

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

        cp $SCXA_NEXTFLOW/\$SUBDIR/.nextflow.log matrices/aggregation.log
   """
    
}

KALLISTO_COUNT_MATRIX.into{
    KALLISTO_COUNT_MATRIX_FOR_SCANPY
    KALLISTO_COUNT_MATRIX_FOR_BUNDLE
}

process add_reference_for_scanpy {

    input:
        set val(species), file(countMatrix) from KALLISTO_COUNT_MATRIX_FOR_SCANPY
    
    output:
        set val(species), file(countMatrix), stdout into KALLISTO_COUNT_MATRIX_FOR_SCANPY_WITH_REF

    """
    cat ${baseDir}/conf/reference/${species}.conf | grep 'gtf' | grep -o "'.*'" | sed "s/'//g" | tr -d \'\\n\'
    """

}

// Run Scanpy with https://github.com/ebi-gene-expression-group/scanpy-workflow

process scanpy {
    
    conda 'nextflow'

    storeDir "$SCXA_RESULTS/$exp_name/$species/scanpy"
    
    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus == 130  ? 'retry' : 'finish' }
    maxRetries 20
    
    input:
        set val(species), file(countMatrix), val(referenceGtf) from KALLISTO_COUNT_MATRIX_FOR_SCANPY_WITH_REF
        set val(species), file (confFile), file(sdrfFile) from COMBINED_CONFIG_FOR_SCANPY

    output:
        set val(species), file("matrices/*_filter_cells_genes.zip") into FILTERED_MATRIX
        set val(species), file("matrices/*_normalised.zip") into NORMALISED_MATRIX
        set val(species), file("pca/*") into PCA
        set val(species), file("clustering/*") into CLUSTERING
        set val(species), file("umap/*") into UMAP
        set val(species), file("tsne/*") into TSNE
        set val(species), file("markers/*") into MARKERS

    """
        RESULTS_ROOT=\$PWD
        SUBDIR="$exp_name/$species/scanpy"     

        mkdir -p $SCXA_WORK/\$SUBDIR
        mkdir -p $SCXA_NEXTFLOW/\$SUBDIR
        mkdir -p $SCXA_RESULTS/\$SUBDIR/reports
        pushd $SCXA_NEXTFLOW/\$SUBDIR > /dev/null

        nextflow run \
            -config \$RESULTS_ROOT/$confFile \
            --resultsRoot \$RESULTS_ROOT \
            --gtf \$SCXA_DATA/reference/${referenceGtf} \
            --matrix \$RESULTS_ROOT/${countMatrix} \
            -resume \
            scanpy-workflow \
            -work-dir $SCXA_WORK/\$SUBDIR \
            -with-report $SCXA_RESULTS/\$SUBDIR/reports/report.html \
            -N $SCXA_REPORT_EMAIL \
            -with-dag $SCXA_RESULTS/\$SUBDIR/reports/flowchart.pdf

        if [ \$? -ne 0 ]; then
            echo "Workflow failed for $exp_name - $species - scanpy-workflow" 1>&2
            exit 1
        fi
        
        popd > /dev/null

        cp $SCXA_NEXTFLOW/\$SUBDIR/.nextflow.log scanpy.log
   """
    
}

// Make a bundle from the Scanpy outputs

//process bundle {
    
//    input:
//        set val(species), file('*') from FILTERED_MATRIX
//        set val(species), file('*') from NORMALISED_MATRIX
//        set val(species), file("*") into CLUSTERING
//        set val(species), file("**") into TSNE
//        set val(species), file("*") into MARKERS
        
//    output:
//        file('bundle/MANIFEST')
//        file('bunlde/software.tsv')
//        
        

    
//}
