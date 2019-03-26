#!/usr/bin/env nextflow

sdrfDir = params.sdrfDir

enaSshUser = 'null'

if ( params.containsKey('enaSshUser') ){
    enaSshUser = params.enaSshUser
}

// If user has supplied an experiment ID, then just run for that experiment.
// Otherwise, watch the SDRF directory for incoming SDRF files

if ( params.containsKey('expName')){
    SDRF = Channel.fromPath("${sdrfDir}/${params.expName}.sdrf.txt", checkIfExists: true)
}else{
    SDRF = Channel.fromPath("${sdrfDir}/*.sdrf.txt", checkIfExists: true)
}

// Determine which SDRF files have up-to-date bundles. For new/ updated
// experiments, output to the channel with the relevant IDF

process find_new_updated {

    executor 'local'

    cache 'lenient'

    input:
        file(sdrfFile) from SDRF

    output:
        set stdout, file(sdrfFile), file("${sdrfFile.getSimpleName()}.idf.txt") optional true into SDRF_IDF
        file('bundleLines.txt') optional true into OLD_BUNDLE_LINES

    """
        expName=\$(echo $sdrfFile | awk -F'.' '{print \$1}') 
        
        newExperiment=1
        bundleLogs=\$(ls \$SCXA_RESULTS/\$expName/*/bundle.log 2>/dev/null || true)
        
        if [ -n "\$bundleLogs" ]; then
            newExperiment=0
            while read -r bundleLog; do
                if [ $sdrfFile -nt "\$bundleLog" ]; then
                    newExperiment=1
                else
                    species=\$(echo \$bundleLog | awk -F'/' '{print \$(NF-1)}' | tr -d \'\\n\')
                    echo -e "\$expName\\t\$species\\t$SCXA_RESULTS/\$expName/\$species/bundle" > bundleLines.txt
                fi
            done <<< "\$(echo -e "\$bundleLogs")"
        fi

        if [ \$newExperiment -eq 1 ]; then
            cp $sdrfDir/\${expName}.idf.txt .
            echo \$expName | tr -d \'\\n\'
        fi
    """
}

// Derive the config files. We cache based on input file size and path
// ('lenient'). If we do run, we remove the downstream stored results,
// triggering the sub-workflows (not normally re-run). These will then check
// their own caches and re-run where required. 

process generate_config {

    cache 'lenient'
    
    publishDir "$SCXA_CONF/study", mode: 'copy', overwrite: true

    conda 'r-optparse r-data.table r-workflowscriptscommon'

    input:
        set val(expName), file(sdrfFile), file(idfFile) from SDRF_IDF

    output:
        set val(expName), file ('*.conf'), file('*.sdrf.txt') into CONFIG_FILES
        
    """
    for stage in quantification aggregation scanpy bundle; do
        rm -rf $SCXA_RESULTS/$expName/*/\$stage
    done

    sdrfToNfConf.R \
        --sdrf=$sdrfFile \
        --idf=$idfFile \
        --name=$expName \
        --verbose

    for f in *.conf; do
        echo "includeConfig '${baseDir}/params.config'" | cat - \$f > \$f.tmp && mv \$f.tmp \$f
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
    
    publishDir "$SCXA_RESULTS/$expName/$species/reference", mode: 'copy', overwrite: true

    errorStrategy { task.attempt<=3 ? 'retry' : 'finish' }

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

REFERENCE_GTF.into{
    REFERENCE_GTF_FOR_AGGREGATION
    REFERENCE_GTF_FOR_SCANPY
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

    // Make a configuration for the Fastq provider, and make initial assessment of
    // the available ENA download methods. 

    process configure_download {
        
        conda "${baseDir}/envs/atlas-fastq-provider.yml"
            
        publishDir "$SCXA_RESULTS/", mode: 'copy', overwrite: true
        
        output:
            file('atlas_fastq_provider_config.sh') into DOWNLOAD_CONFIG

        script:
            downloadConfig = file('atlas_fastq_provider_config.sh')
            sshUser=''
            if (params.containsKey('enaSshUser')){
                sshUser=params.enaSshUser
            }      

            """
                echo ENA_RETRIES=\\'${params.downloadRetries}\\' > download_config.sh
                echo FETCH_FREQ_MILLIS=\\'${params.fetchFreqMillis}\\' >> download_config.sh
                echo FASTQ_PROVIDER_TEMPDIR=\\'$NXF_TEMP/atlas-fastq-provider\\' >> download_config.sh
                echo ALLOWED_DOWNLOAD_METHODS=\\'${params.allowedDownloadMethods}\\' >> download_config.sh

                if [ -n "$sshUser" ]; then
                    echo ENA_SSH_USER=\\'${params.enaSshUser}\\' >> download_config.sh
                fi

                initialiseEnaProbe.sh -c download_config.sh
                cp download_config.sh atlas_fastq_provider_config.sh 
            """
    }

    // Inialise the probe file describing download method performance

    process initialise_download_probe {
        
        conda "${baseDir}/envs/atlas-fastq-provider.yml"

        cache false

        publishDir "$NXF_TEMP/atlas-fastq-provider/", mode: 'copy', overwrite: true

        input:
            file(downloadConfig) from DOWNLOAD_CONFIG

        output:
            file(downloadConfig), file("$NXF_TEMP/fastq_provider.probe") into DOWNLOAD_CONFIG_INIT

        """
        initialiseEnaProbe.sh -c ${downloadConfig} -t fastq_provider.probe
        """
    }

    // Run quantification with https://github.com/ebi-gene-expression-group/scxa-smartseq-quantification-workflow

    process quantify {

        maxForks params.maxConcurrentQuantifications

        conda "${baseDir}/envs/nextflow.yml"

        publishDir "$SCXA_RESULTS/$expName/$species/quantification", mode: 'copy', overwrite: true
        
        memory { 4.GB * task.attempt }
        errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }
        maxRetries 20
        
        input:
            set val(expName), val(species), file (confFile), file(sdrfFile) from COMBINED_CONFIG_FOR_QUANTIFY
            set val(expName), val(species), file(referenceFasta) from REFERENCE_FASTA
            set val(expName), val(species), val(contaminationIndex) from CONTAMINATION_INDEX
            set file(downloadConfig), file(downloadProbe) from DOWNLOAD_CONFIG_INIT

        output:
            set val(expName), val(species), file ("kallisto") into KALLISTO_DIRS 
            set val(expName), val(species), file ("qc") into QUANT_QC
            set val(expName), val(species), file('quantification.log')    

        """
            for stage in aggregation scanpy bundle; do
                rm -rf $SCXA_RESULTS/$expName/$species/\$stage
            done

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
                --downloadConfig $downloadConfig \
                --referenceFasta \$RESULTS_ROOT/$referenceFasta \
                --contaminationIndex $contaminationIndex \
                --resultsRoot \$RESULTS_ROOT \
                --enaSshUser $enaSshUser \
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

    publishDir "$SCXA_RESULTS/$expName/$species/aggregation", mode: 'copy', overwrite: true
    
    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }
    maxRetries 20
    
    input:
        set val(expName), val(species), file (confFile), file(sdrfFile) from COMBINED_CONFIG_FOR_AGGREGATION
        set val(expName), val(species), file ('kallisto') from KALLISTO_DIRS
        set val(expName), val(species), file(referenceGtf) from REFERENCE_GTF_FOR_AGGREGATION

    output:
        set val(expName), val(species), file("matrices/*_counts.zip") into KALLISTO_COUNT_MATRIX
        set val(expName), val(species), file("matrices/*_tpm.zip") into KALLISTO_ABUNDANCE_MATRIX
        set val(expName), val(species), file("matrices/*.stats.tsv") into KALLISTO_STATS
        set val(expName), val(species), file('matrices/aggregation.log')    

    """
        for stage in scanpy bundle; do
            rm -rf $SCXA_RESULTS/$expName/$species/\$stage
        done
        
        RESULTS_ROOT=\$PWD
        SUBDIR="$expName/$species/aggregation"     

        mkdir -p $SCXA_WORK/\$SUBDIR
        mkdir -p $SCXA_NEXTFLOW/\$SUBDIR
        mkdir -p $SCXA_RESULTS/\$SUBDIR/reports
        pushd $SCXA_NEXTFLOW/\$SUBDIR > /dev/null

        nextflow run \
            -config \$RESULTS_ROOT/$confFile \
            --resultsRoot \$RESULTS_ROOT \
            --quantDir \$RESULTS_ROOT/kallisto \
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

// Run Scanpy with https://github.com/ebi-gene-expression-group/scanpy-workflow

process scanpy {
    
    conda "${baseDir}/envs/nextflow.yml"

    publishDir "$SCXA_RESULTS/$expName/$species/scanpy", mode: 'copy', overwrite: true
    
    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }
    maxRetries 20
    
    input:
        set val(expName), val(species), file(countMatrix) from KALLISTO_COUNT_MATRIX_FOR_SCANPY
        set val(expName), val(species), file (confFile), file(sdrfFile) from COMBINED_CONFIG_FOR_SCANPY
        set val(expName), val(species), file(referenceGtf) from REFERENCE_GTF_FOR_SCANPY

    output:
        set val(expName), val(species), file("matrices/*_filter_cells_genes.zip") into FILTERED_MATRIX
        set val(expName), val(species), file("matrices/*_normalised.zip") into NORMALISED_MATRIX
        set val(expName), val(species), file("pca") into PCA
        set val(expName), val(species), file("clustering/clusters.txt") into CLUSTERS
        set val(expName), val(species), file("umap") into UMAP
        set val(expName), val(species), file("tsne") into TSNE
        set val(expName), val(species), file("markers") into MARKERS
        file('scanpy.log')

    """
        rm -rf $SCXA_RESULTS/$expName/$species/bundle
        
        RESULTS_ROOT=\$PWD
        SUBDIR="$expName/$species/scanpy"     

        mkdir -p $SCXA_WORK/\$SUBDIR
        mkdir -p $SCXA_NEXTFLOW/\$SUBDIR
        mkdir -p $SCXA_RESULTS/\$SUBDIR/reports
        pushd $SCXA_NEXTFLOW/\$SUBDIR > /dev/null

        nextflow run \
            -config \$RESULTS_ROOT/$confFile \
            --resultsRoot \$RESULTS_ROOT \
            --gtf \$RESULTS_ROOT/${referenceGtf} \
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

    publishDir "$SCXA_RESULTS/$expName/$species", mode: 'copy', overwrite: true
    
    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }
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
        file('bundle.log')
        file('bundleLines.txt') into NEW_BUNDLE_LINES
        
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

        cp $SCXA_NEXTFLOW/\$SUBDIR/.nextflow.log bundle.log

        echo -e "$expName\\t$species\\t$SCXA_RESULTS/$expName/$species/bundle" > bundleLines.txt
   """
    
}

// Record the completed bundles

OLD_BUNDLE_LINES
    .concat ( NEW_BUNDLE_LINES )
    .collectFile(name: 'all.done.txt', sort: true, storeDir: "$SCXA_RESULTS")
    .set{
        BUNDLE_LINES
    }

// For each experiment and species with a completed bundle, remove the work dir

process cleanup {
    
    input:
        file(bundleLines) from BUNDLE_LINES
    
    output:
        file('.cleaned')

    """
    cat ${bundleLines} | while read -r l; do
        expName=\$(echo "\$l" | awk '{print \$1}')
        species=\$(echo "\$l" | awk '{print \$2}')
        rm -rf $SCXA_WORK/\$expName/\$species
    done
    touch .cleaned
    """
}

