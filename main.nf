#!/usr/bin/env nextflow

sdrfDir = params.sdrfDir

// If user has supplied an experiment ID, then just run for that experiment.
// Otherwise, watch the SDRF directory for incoming SDRF files

if ( params.containsKey('expName')){
    SDRF = Channel.fromPath("${sdrfDir}/${params.expName}.sdrf.txt", checkIfExists: true)
}else{
    Channel
        .watchPath( "${sdrfDir}/*.sdrf.txt", 'create,modify' )
        .set{ SDRF }
}

// Locate matching IDF files and determine experiment ID

process find_idf {
    
    input:
        file(sdrfFile) from SDRF
   
    output:
        set stdout, file(sdrfFile), file("${sdrfFile.getSimpleName()}.idf.txt") into SDRF_IDF

    """
        expName=\$(echo $sdrfFile | awk -F'.' '{print \$1}') 
        cp $sdrfDir/\${expName}.idf.txt .
        echo \$expName | tr -d \'\\n\'
    """
}

// Derive the config files

process generate_config {

    publishDir "$SCXA_CONF/study", mode: 'copy', overwrite: true

    conda 'r-optparse r-data.table r-workflowscriptscommon'

    input:
        set val(expName), file(sdrfFile), file(idfFile) from SDRF_IDF

    output:
        set val(expName), file ('*.conf'), file('*.sdrf.txt') into CONFIG_FILES
        
    """
    sdrfToNfConf.R \
        --sdrf=$sdrfFile \
        --idf=$idfFile \
        --name=$expName \
        --verbose

    for f in *.conf; do
        echo -e "includeConfig '${baseDir}/params.config'\n"|cat - \$f > \$f.tmp && mv \$f.tmp \$f
    done
    """
}

// Mark config files with species

process mark_conf_species {
   
    conda 'pyyaml' 
 
    input:
        set val(expName), file(confFile), file(sdrfFile) from CONFIG_FILES

    output:
        set val(expName), stdout, file(confFile), file(sdrfFile) into CONF_BY_SPECIES

    """
    parseNfConfig.py --paramFile $confFile --paramKeys params,organism
    """
}

CONF_BY_SPECIES
    .into{
        COMBINED_CONFIG_FOR_REFERENCE
        COMBINED_CONFIG_FOR_QUANTIFY
        COMBINED_CONFIG_FOR_AGGREGATION
        COMBINED_CONFIG_FOR_SCANPY
    }

// Prepare a reference depending on spikes

process prepare_reference {

    conda 'pyyaml' 
    
    storeDir "$SCXA_RESULTS/$expName/$species/reference"

    errorStrategy { task.attempt<=3 ? 'retry' : 'finish' }

    cache true

    input:
        set val(expName), val(species), file(confFile), file(sdrfFile) from COMBINED_CONFIG_FOR_REFERENCE
    
    output:
        set val(expName), val(species), file("reference.fastq.gz") into REFERENCE_FASTA
        set val(expName), val(species), file("reference.gtf.gz") into REFERENCE_GTF
        set val(expName), val(species), stdout into CONTAMINATION_INDEX

    """
    species_conf=$SCXA_PRE_CONF/reference/${species}.conf
    cdna_fasta=$SCXA_DATA/reference/\$(parseNfConfig.py --paramFile \$species_conf --paramKeys params,reference,cdna)
    cdna_gtf=$SCXA_DATA/reference/\$(parseNfConfig.py --paramFile \$species_conf --paramKeys params,reference,gtf)
    spikes=\$(parseNfConfig.py --paramFile $confFile --paramKeys params,spikes)

    if [ \$spikes != 'None' ]; then
        spikes_conf="$SCXA_PRE_CONF/reference/\${spikes}.conf"
        spikes_fasta=$SCXA_DATA/reference/\$(parseNfConfig.py --paramFile \$spikes_conf --paramKeys params,reference,spikes,cdna)
        spikes_gtf=$SCXA_DATA/reference/\$(parseNfConfig.py --paramFile \$spikes_conf --paramKeys params,reference,spikes,gtf)
        
        cat \$cdna_fasta \$spikes_fasta > reference.fastq.gz
        cat \$cdna_gtf \$spikes_gtf > reference.gtf.gz
    else
        ln -s \$cdna_fasta reference.fastq.gz
        ln -s \$cdna_gtf reference.gtf.gz
    fi

    contamination_index=\$(parseNfConfig.py --paramFile \$species_conf --paramKeys params,reference,contamination_index)
    if [ \$contamination_index != 'None' ]; then
        printf $SCXA_DATA/contamination/\$contamination_index
    else
        echo None
    fi

    """
}

// Allow a forcible skipping of the quantification phase. Useful if we've
// modified the workflows in some way that would make NF recompute, but we know
// that's not necessary

if ( params.containsKey('skipQuantification') && params.skipQuantification == 'yes'){

    process spoof_quantify {

        input:
            set val(expName), val(species), file (confFile), file(sdrfFile) from COMBINED_CONFIG_FOR_QUANTIFY
            set val(expName), val(species), file(referenceFasta) from REFERENCE_FASTA
            set val(expName), val(species), val(contaminationIndex) from CONTAMINATION_INDEX

        output:
            set val(expName), val(species), file ("kallisto/*") into KALLISTO_DIRS 
            set val(expName), val(species), file ("qc/*") into QUANT_QC 

        """
            ln -s $SCXA_RESULTS/$expName/$species/quantification/kallisto .
            ln -s $SCXA_RESULTS/$expName/$species/quantification/qc .
        """
    }
}else{

    // Run quantification with https://github.com/ebi-gene-expression-group/scxa-smartseq-quantification-workflow

    process quantify {

        maxForks 1

        conda "${baseDir}/envs/nextflow.yml"

        storeDir "$SCXA_RESULTS/$expName/$species/quantification"
        
        memory { 4.GB * task.attempt }
        errorStrategy { task.exitStatus == 130  ? 'retry' : 'finish' }
        maxRetries 20
        
        input:
            set val(expName), val(species), file (confFile), file(sdrfFile) from COMBINED_CONFIG_FOR_QUANTIFY
            set val(expName), val(species), file(referenceFasta) from REFERENCE_FASTA
            set val(expName), val(species), val(contaminationIndex) from CONTAMINATION_INDEX

        output:
            set val(expName), val(species), file ("kallisto/*") into KALLISTO_DIRS 
            set val(expName), val(species), file ("qc/*") into QUANT_QC
            set val(expName), val(species), file('quantification.log')    

        """
            protocol=\$(parseNfConfig.py --paramFile $confFile --paramKeys params,sc_protocol)

            echo \$protocol | grep "smart-seq" > /dev/null
            if [ \$? -eq 0 ]; then
                quantification_workflow=scxa-smartseq-quantification-workflow
            else
                echo "This is not a SMART-seq experiment" 1>&2
                exit 1
            fi

            RESULTS_ROOT=\$PWD
            SUBDIR="$expName/$species/quantification"     

            mkdir -p $SCXA_WORK/\$SUBDIR
            mkdir -p $SCXA_NEXTFLOW/\$SUBDIR
            mkdir -p $SCXA_RESULTS/\$SUBDIR/reports
            pushd $SCXA_NEXTFLOW/\$SUBDIR > /dev/null

            nextflow run \
                -config \$RESULTS_ROOT/$confFile \
                --sdrf \$RESULTS_ROOT/$sdrfFile \
                --referenceFasta \$RESULTS_ROOT/$referenceFasta \
                --contaminationIndex $contaminationIndex \
                --resultsRoot \$RESULTS_ROOT \
                -resume \
                \$quantification_workflow \
                -work-dir $SCXA_WORK/\$SUBDIR \
                -with-report $SCXA_RESULTS/\$SUBDIR/reports/report.html \
                -N $SCXA_REPORT_EMAIL \
                -with-dag $SCXA_RESULTS/\$SUBDIR/reports/flowchart.pdf

            if [ \$? -ne 0 ]; then
                echo "Workflow failed for $expName - $species - \$quantification_workflow" 1>&2
                exit 1
            fi
            
            popd > /dev/null

            cp $SCXA_NEXTFLOW/\$SUBDIR/.nextflow.log quantification.log
       """
    }
}

// Run aggregation with https://github.com/ebi-gene-expression-group/scxa-aggregation-workflow

process aggregate {
    
    conda "${baseDir}/envs/nextflow.yml"

    storeDir "$SCXA_RESULTS/$expName/$species/aggregation"
    
    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus == 130  ? 'retry' : 'finish' }
    maxRetries 20
    
    input:
        set val(expName), val(species), file (confFile), file(sdrfFile) from COMBINED_CONFIG_FOR_AGGREGATION
        set val(expName), val(species), file ('*') from KALLISTO_DIRS.collect()
        set val(expName), val(species), file(referenceGtf) from REFERENCE_GTF

    output:
        set val(expName), val(species), file("matrices/*_counts.zip") into KALLISTO_COUNT_MATRIX
        set val(expName), val(species), file("matrices/*_tpm.zip") into KALLISTO_ABUNDANCE_MATRIX
        set val(expName), val(species), file("matrices/*.stats.tsv") into KALLISTO_STATS
        set val(expName), val(species), file('matrices/aggregation.log')    

    """
        RESULTS_ROOT=\$PWD
        SUBDIR="$expName/$species/aggregation"     

        mkdir -p $SCXA_WORK/\$SUBDIR
        mkdir -p $SCXA_NEXTFLOW/\$SUBDIR
        mkdir -p $SCXA_RESULTS/\$SUBDIR/reports
        pushd $SCXA_NEXTFLOW/\$SUBDIR > /dev/null

        nextflow run \
            -config \$RESULTS_ROOT/$confFile \
            --resultsRoot \$RESULTS_ROOT \
            --referenceGtf ${referenceGtf} \
            -resume \
            scxa-aggregation-workflow \
            -work-dir $SCXA_WORK/\$SUBDIR \
            -with-report $SCXA_RESULTS/\$SUBDIR/reports/report.html \
            -N $SCXA_REPORT_EMAIL \
            -with-dag $SCXA_RESULTS/\$SUBDIR/reports/flowchart.pdf

        if [ \$? -ne 0 ]; then
            echo "Workflow failed for $expName - $species - scxa_aggregation_workflow" 1>&2
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
        set val(expName), val(species), file(countMatrix) from KALLISTO_COUNT_MATRIX_FOR_SCANPY
    
    output:
        set val(expName), val(species), file(countMatrix), stdout into KALLISTO_COUNT_MATRIX_FOR_SCANPY_WITH_REF

    """
    cat ${baseDir}/conf/reference/${species}.conf | grep 'gtf' | grep -o "'.*'" | sed "s/'//g" | tr -d \'\\n\'
    """

}

// Run Scanpy with https://github.com/ebi-gene-expression-group/scanpy-workflow

process scanpy {
    
    conda "${baseDir}/envs/nextflow.yml"

    storeDir "$SCXA_RESULTS/$expName/$species/scanpy"
    
    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus == 130  ? 'retry' : 'finish' }
    maxRetries 20
    
    input:
        set val(expName), val(species), file(countMatrix), val(referenceGtf) from KALLISTO_COUNT_MATRIX_FOR_SCANPY_WITH_REF
        set val(expName), val(species), file (confFile), file(sdrfFile) from COMBINED_CONFIG_FOR_SCANPY

    output:
        set val(expName), val(species), file("matrices/*_filter_cells_genes.zip") into FILTERED_MATRIX
        set val(expName), val(species), file("matrices/*_normalised.zip") into NORMALISED_MATRIX
        set val(expName), val(species), file("pca") into PCA
        set val(expName), val(species), file("clustering/clusters.txt") into CLUSTERS
        set val(expName), val(species), file("umap") into UMAP
        set val(expName), val(species), file("tsne") into TSNE
        set val(expName), val(species), file("markers") into MARKERS

    """
        RESULTS_ROOT=\$PWD
        SUBDIR="$expName/$species/scanpy"     

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
            echo "Workflow failed for $expName - $species - scanpy-workflow" 1>&2
            exit 1
        fi
        
        popd > /dev/null

        cp $SCXA_NEXTFLOW/\$SUBDIR/.nextflow.log scanpy.log
   """
    
}

// Make a bundle from the Scanpy outputs

process bundle {
    
    conda "${baseDir}/envs/nextflow.yml"

    storeDir "$SCXA_RESULTS/$expName/$species"
    
    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus == 130  ? 'retry' : 'finish' }
    maxRetries 20
    
    input:
        set val(expName), val(species), file(filteredMatrix) from FILTERED_MATRIX
        set val(expName), val(species), file(normalisedMatrix) from NORMALISED_MATRIX
        set val(expName), val(species), file(tpmMatrix) from KALLISTO_ABUNDANCE_MATRIX
        set val(expName), val(species), file(clusters) from CLUSTERS
        set val(expName), val(species), file('*') from TSNE
        set val(expName), val(species), file('*') from MARKERS
        
    output:
        file('bundle/*')
        
    """    
        RESULTS_ROOT=\$PWD
        SUBDIR="$expName/$species/bundle"     

        mkdir -p $SCXA_WORK/\$SUBDIR
        mkdir -p $SCXA_NEXTFLOW/\$SUBDIR
        mkdir -p $SCXA_RESULTS/\$SUBDIR/reports
        pushd $SCXA_NEXTFLOW/\$SUBDIR > /dev/null

        nextflow run \
            --resultsRoot \$RESULTS_ROOT \
            --rawFilteredMatrix ${filteredMatrix} \
            --normalisedMatrix ${normalisedMatrix} \
            --tpmMatrix ${tpmMatrix} \
            --clusters ${clusters} \
            --tsneDir tsne \
            --markersDir markers \
            --softwareTemplate ${baseDir}/conf/smartseq.software.tsv \
            -resume \
            scxa-bundle-workflow \
            -work-dir $SCXA_WORK/\$SUBDIR \
            -with-report $SCXA_RESULTS/\$SUBDIR/reports/report.html \
            -N $SCXA_REPORT_EMAIL \
            -with-dag $SCXA_RESULTS/\$SUBDIR/reports/flowchart.pdf

        if [ \$? -ne 0 ]; then
            echo "Workflow failed for $expName - $species - scanpy-workflow" 1>&2
            exit 1
        fi
        
        popd > /dev/null

        cp $SCXA_NEXTFLOW/\$SUBDIR/.nextflow.log scanpy.log
   """
    
}
