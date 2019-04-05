#!/usr/bin/env nextflow

sdrfDir = params.sdrfDir
tertiaryWorkflow = params.tertiaryWorkflow

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

    conda "${baseDir}/envs/atlas-fastq-provider.yml"
    
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

COMBINED_CONFIG_FOR_QUANTIFY
    .join( REFERENCE_FASTA, by: [0,1] )
    .join( CONTAMINATION_INDEX, by: [0,1] )
    .set { QUANTIFICATION_INPUTS }

// Select workflow for each experiment based on protocols

process select_workflow{

    conda 'pyyaml' 
    
    input:
        set val(expName), val(species), file (confFile), file(sdrfFile), file(referenceFasta), val(contaminationIndex) from QUANTIFICATION_INPUTS

    output:
        set val(expName), val(species), stdout, file (confFile), file(sdrfFile), file(referenceFasta), val(contaminationIndex) into QUANTIFICATION_INPUTS_BY_WF

    """
    protocol=\$(parseNfConfig.py --paramFile $confFile --paramKeys params,sc_protocol)

    if [ "\$protocol" == 'smart-seq' ] ||  [ "\$protocol" == 'smart-seq2' ] || [ "\$protocol" == 'smarter' ] || [ "\$protocol" == 'smart-like' ]; then
        echo smartseq
    elif [ "\$protocol" == '10xv2' ]  [ "\$protocol" == 'drop-seq' ] || [ "\$protocol" == 'smart-seq' ]; then
        echo droplet
    else
        echo "Can't currently handle \$protocol experiments" 1>&2
        exit 1
    fi
    """
}

// Separate droplet and smart-type experiments

DROPLET = Channel.create()
SMART = Channel.create()

QUANTIFICATION_INPUTS_BY_WF.choice( SMART, DROPLET ) {a -> 
    a[2] == 'droplet' ? 1 : 0
}

// Allow a forcible skipping of the quantification phase. Useful if we've
// modified the workflows in some way that would make NF recompute, but we know
// that's not necessary

if ( params.containsKey('skipQuantification') && params.skipQuantification == 'yes'){

    process spoof_smart_quantify {

        input:
            set val(expName), val(species), file(wfType), file (confFile), file(sdrfFile), file(referenceFasta), val(contaminationIndex) from SMART

        output:
            set val(expName), val(species), val(wfType), file ("kallisto/*") into SMART_KALLISTO_DIRS 
            set val(expName), val(species), val(wfType), file ("qc/*") into SMART_QUANT_QC 

        """
            ln -s $SCXA_RESULTS/$expName/$species/quantification/kallisto .
            ln -s $SCXA_RESULTS/$expName/$species/quantification/qc .
        """
    }
}else{

    // Make a configuration for the Fastq provider, and make initial assessment
    // of the available ENA download methods. This establishes a mock
    // dependency with the quantify process, as per
    // http://nextflow-io.github.io/patterns/index.html#_mock_dependency 

    process initialise_downloader {
        
        conda "${baseDir}/envs/atlas-fastq-provider.yml"

        cache false

        output:
            val true into INIT_DONE

        """
        initialiseDownload.sh ${baseDir}/params.config $NXF_TEMP/atlas-fastq-provider/download_config.sh $NXF_TEMP/atlas-fastq-provider/fastq_provider.probe 
        """
    }

    INIT_DONE
        .into{
            INIT_DONE_SMART
            INIT_DONE_DROPLET
        }

    // Run quantification with https://github.com/ebi-gene-expression-group/scxa-smartseq-quantification-workflow

    process smart_quantify {

        maxForks params.maxConcurrentQuantifications

        conda "${baseDir}/envs/nextflow.yml"

        publishDir "$SCXA_RESULTS/$expName/$species/quantification", mode: 'copy', overwrite: true
        
        memory { 40.GB * task.attempt }
        errorStrategy { task.attempt<=5 ? 'retry' : 'finish' }
        
        input:
            set val(expName), val(species), val(wfType), file(confFile), file(sdrfFile), file(referenceFasta), val(contaminationIndex) from SMART
            val flag from INIT_DONE_SMART

        output:
            set val(expName), val(species), val(wfType), file ("kallisto") into SMART_KALLISTO_DIRS 
            set val(expName), val(species), val(wfType), file ("qc") into SMART_QUANT_QC
            set val(expName), val(species), val(wfType), file('quantification.log')    

        """
            for stage in aggregation scanpy bundle; do
                rm -rf $SCXA_RESULTS/$expName/$species/\$stage
            done

            RESULTS_ROOT=\$PWD
            SUBDIR="$expName/$species/quantification"     

            mkdir -p $SCXA_WORK/\$SUBDIR
            mkdir -p $SCXA_NEXTFLOW/\$SUBDIR
            mkdir -p $SCXA_RESULTS/\$SUBDIR/reports
            pushd $SCXA_NEXTFLOW/\$SUBDIR > /dev/null

            BRANCH=''
            if [ -n "$SCXA_BRANCH" ]; then
                BRANCH="-r $SCXA_BRANCH"
            fi

            nextflow run \
                -config \$RESULTS_ROOT/$confFile \
                --sdrf \$RESULTS_ROOT/$sdrfFile \
                --referenceFasta \$RESULTS_ROOT/$referenceFasta \
                --contaminationIndex $contaminationIndex \
                --resultsRoot \$RESULTS_ROOT \
                --enaSshUser $enaSshUser \
                -resume \
                \$quantification_workflow \
                -work-dir $SCXA_WORK/\$SUBDIR \
                -with-report $SCXA_RESULTS/\$SUBDIR/reports/report.html \
                -with-trace  $SCXA_RESULTS/\$SUBDIR/reports/trace.txt \
                -N $SCXA_REPORT_EMAIL \
                -with-dag $SCXA_RESULTS/\$SUBDIR/reports/flowchart.pdf \
                \$BRANCH

            if [ \$? -ne 0 ]; then
                echo "Workflow failed for $expName - $species - \$quantification_workflow" 1>&2
                exit 1
            fi
            
            popd > /dev/null

            cp $SCXA_NEXTFLOW/\$SUBDIR/.nextflow.log quantification.log
       """
    }
    
    // Run quantification with https://github.com/ebi-gene-expression-group/scxa-droplet-quantification-workflow

    process droplet_quantify {

        maxForks params.maxConcurrentQuantifications
    
        conda "${baseDir}/envs/nextflow.yml"

        publishDir "$SCXA_RESULTS/$expName/$species/quantification", mode: 'copy', overwrite: true
        
        memory { 40.GB * task.attempt }
        errorStrategy { task.attempt<=5 ? 'retry' : 'finish' }

        input:
            set val(expName), val(species), val(wfType), file (confFile), file(sdrfFile), file(referenceFasta), val(contaminationIndex) from DROPLET
            val flag from INIT_DONE_DROPLET

        output:
            set val(expName), val(species), val(wfType), file("matrices/*_counts.zip"), file("NOTPM") into DROPLET_MATRICES

        """
        Call alevin workflow

        touch NOTPM
        """
    }

}

// Run aggregation with https://github.com/ebi-gene-expression-group/scxa-aggregation-workflow

COMBINED_CONFIG_FOR_AGGREGATION
    .join(SMART_KALLISTO_DIRS, by: [0,1])
    .join(REFERENCE_GTF_FOR_AGGREGATION, by: [0,1])
    .set{SMART_AGGREGATION_INPUTS}

process kallisto_aggregate {
    
    conda "${baseDir}/envs/nextflow.yml"

    publishDir "$SCXA_RESULTS/$expName/$species/aggregation", mode: 'copy', overwrite: true
    
    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }
    maxRetries 20
    
    input:
        set val(expName), val(species), file (confFile), file(sdrfFile), file(wfType), file ('kallisto'), file(referenceGtf) from SMART_AGGREGATION_INPUTS   

    output:
        set val(expName), val(species), val(wfType), file("matrices/*_counts.zip"), file("matrices/*_tpm.zip") into SMART_MATRICES
        set val(expName), val(species), val(wfType),  file("matrices/*.stats.tsv") into SMART_KALLISTO_STATS
        set val(expName), val(species), val(wfType),  file('matrices/aggregation.log')    

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

        BRANCH=''
        if [ -n "$SCXA_BRANCH" ]; then
            BRANCH="-r $SCXA_BRANCH"
        fi

        nextflow run \
            -config \$RESULTS_ROOT/$confFile \
            --resultsRoot \$RESULTS_ROOT \
            --quantDir \$RESULTS_ROOT/kallisto \
            --referenceGtf ${referenceGtf} \
            -resume \
            scxa-aggregation-workflow \
            -work-dir $SCXA_WORK/\$SUBDIR \
            -with-report $SCXA_RESULTS/\$SUBDIR/reports/report.html \
            -with-trace  $SCXA_RESULTS/\$SUBDIR/reports/trace.txt \
            -N $SCXA_REPORT_EMAIL \
            -with-dag $SCXA_RESULTS/\$SUBDIR/reports/flowchart.pdf \
            \$BRANCH

        if [ \$? -ne 0 ]; then
            echo "Workflow failed for $expName - $species - scxa_aggregation_workflow" 1>&2
            exit 1
        fi
        
        popd > /dev/null

        cp $SCXA_NEXTFLOW/\$SUBDIR/.nextflow.log matrices/aggregation.log
   """
    
}

SMART_MATRICES
    .concat(DROPLET_MATRICES)
    .into{
        MATRICES_FOR_SCANPY
        MATRICES_FOR_BUNDLE
    }

// Set the inputs for tertiary analysis

MATRICES_FOR_SCANPY
    .join(COMBINED_CONFIG_FOR_SCANPY, by: [0,1])
    .join(REFERENCE_GTF_FOR_SCANPY, by: [0,1])
    .set { TERTIARY_INPUTS }

// Run Scanpy with https://github.com/ebi-gene-expression-group/scanpy-workflow

if ( tertiaryWorkflow == 'scanpy-workflow'){

    process scanpy {
        
        conda "${baseDir}/envs/nextflow.yml"

        publishDir "$SCXA_RESULTS/$expName/$species/scanpy", mode: 'copy', overwrite: true
        
        memory { 4.GB * task.attempt }
        errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }
        maxRetries 20
          
        input:
            set val(expName), val(species), file(wfType), file(countMatrix), file(tpmMatrix), file (confFile), file(sdrfFile), file(referenceGtf) from TERTIARY_INPUTS

        output:
            set val(expName), val(species), file("matrices/*_filter_cells_genes.zip"), file("matrices/*_normalised.zip"), file("pca"), file("clustering/clusters.txt"), file("umap"), file("tsne"), file("markers") into TERTIARY_RESULTS
            file('scanpy.log')

        """
            rm -rf $SCXA_RESULTS/$expName/$species/bundle
            
            RESULTS_ROOT=\$PWD
            SUBDIR="$expName/$species/scanpy"     

            mkdir -p $SCXA_WORK/\$SUBDIR
            mkdir -p $SCXA_NEXTFLOW/\$SUBDIR
            mkdir -p $SCXA_RESULTS/\$SUBDIR/reports
            pushd $SCXA_NEXTFLOW/\$SUBDIR > /dev/null
                
            BRANCH=''
            if [ -n "$SCXA_BRANCH" ]; then
                BRANCH="-r $SCXA_BRANCH"
            fi

            nextflow run \
                -config \$RESULTS_ROOT/$confFile \
                --resultsRoot \$RESULTS_ROOT \
                --gtf \$RESULTS_ROOT/${referenceGtf} \
                --matrix \$RESULTS_ROOT/${countMatrix} \
                -resume \
                scanpy-workflow \
                -work-dir $SCXA_WORK/\$SUBDIR \
                -with-report $SCXA_RESULTS/\$SUBDIR/reports/report.html \
                -with-trace  $SCXA_RESULTS/\$SUBDIR/reports/trace.txt \
                -N $SCXA_REPORT_EMAIL \
                -with-dag $SCXA_RESULTS/\$SUBDIR/reports/flowchart.pdf \
                \$BRANCH

            if [ \$? -ne 0 ]; then
                echo "Workflow failed for $expName - $species - scanpy-workflow" 1>&2
                exit 1
            fi
            
            popd > /dev/null

            cp $SCXA_NEXTFLOW/\$SUBDIR/.nextflow.log scanpy.log
       """
        
    }
    
}else{
    
    process spoof_tertiary {
        
        input:
            set val(expName), val(species), val(wfType), file(countMatrix), file(tpmMatrix), file (confFile), file(sdrfFile), file(referenceGtf) from TERTIARY_INPUTS

        output:
            set val(expName), val(species), file("NOFILT"), file("NONORM"), file("NOPCA"), file("NOCLUST"), file("NOUMAP"), file("NOTSNE"), file("NOMARKERS") into TERTIARY_RESULTS

        """
            touch NOFILT NONORM NOPCA NOCLUST NOUMAP NOTSNE NOMARKERS
        """
    }    
        
}

// Make a bundle from the outputs

MATRICES_FOR_BUNDLE
    .join(TERTIARY_RESULTS, by: [0,1])
    .set { BUNDLE_INPUTS } 

process bundle {
    
    conda "${baseDir}/envs/nextflow.yml"

    publishDir "$SCXA_RESULTS/$expName/$species", mode: 'copy', overwrite: true
    
    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }
    maxRetries 20
    
    input:
        set val(expName), val(species), val(wfType), file(rawMatrix), file(tpmMatrix), file(filteredMatrix), file(normalisedMatrix), file(pca), file(clusters), file('*'), file('*'), file('*') from BUNDLE_INPUTS
        
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
    
        # Retrieve the original reference file names to report to bundle 
        species_conf=$SCXA_PRE_CONF/reference/${species}.conf
        cdna_fasta=$SCXA_DATA/reference/\$(parseNfConfig.py --paramFile \$species_conf --paramKeys params,reference,cdna)
        cdna_gtf=$SCXA_DATA/reference/\$(parseNfConfig.py --paramFile \$species_conf --paramKeys params,reference,gtf)

        BRANCH=''
        if [ -n "$SCXA_BRANCH" ]; then
            BRANCH="-r $SCXA_BRANCH"
        fi

        TPM_OPTIONS=''
        if [ "$wfType" == 'smart' ]; then
            TPM_OPTIONS="--tpmMatrix ${tpmMatrix}"
        fi

        TERTIARY_OPTIONS=''
        if [ "$tertiaryWorkflow" == 'scanpy_workflow' ]; then
            TERTIARY_OPTIONS = "--tertiaryWorkflow scxa-$tertiaryWorkflow --rawFilteredMatrix ${filteredMatrix} --normalisedMatrix ${normalisedMatrix} --clusters ${clusters} --tsneDir tsne --markersDir markers"
        fi 

        nextflow run \
            --baseWorkflow scxa-${wfType}-quantification-workflow
            --resultsRoot \$RESULTS_ROOT \
            --rawMatrix ${rawMatrix} \$TPM_OPTIONS \
            --referenceFasta \$cdna_fasta \
            --referenceGtf \$cdna_gtf \$TERTIARY_OPTIONS \
            --softwareTemplate ${baseDir}/conf/smartseq.software.tsv \
            -resume \
            scxa-bundle-workflow \
            -work-dir $SCXA_WORK/\$SUBDIR \
            -with-report $SCXA_RESULTS/\$SUBDIR/reports/report.html \
            -with-trace  $SCXA_RESULTS/\$SUBDIR/reports/trace.txt \
            -N $SCXA_REPORT_EMAIL \
            -with-dag $SCXA_RESULTS/\$SUBDIR/reports/flowchart.pdf \
            \$BRANCH

        if [ \$? -ne 0 ]; then
            echo "Workflow failed for $expName - $species - scxa-bundle-workflow" 1>&2
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
    touch $SCXA_WORK/.success
    cat ${bundleLines} | while read -r l; do
        expName=\$(echo "\$l" | awk '{print \$1}')
        species=\$(echo "\$l" | awk '{print \$2}')
        rm -rf $SCXA_WORK/\$expName/\$species
    done
    touch .cleaned
    """
}

