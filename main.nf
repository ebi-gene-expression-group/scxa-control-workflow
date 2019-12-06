#!/usr/bin/env nextflow

tertiaryWorkflow = params.tertiaryWorkflow
overwrite = params.overwrite

dropletProtocols = [ '10xv1', '10xv1a', '10xv1i', '10xv2', '10xv3', 'drop-seq' ]

enaSshUser = 'null'

if ( params.containsKey('enaSshUser') ){
    enaSshUser = params.enaSshUser
}

galaxyCredentials = ''
if ( params.containsKey('galaxyCredentials')){
    galaxyCredentials = params.galaxyCredentials
}

skipQuantification = 'no'
if ( params.containsKey('skipQuantification') && params.skipQuantification == 'yes'){
    skipQuantification = 'yes'
}

skipAggregation = 'no'
if ( params.containsKey('skipAggregation') && params.skipAggregation == 'yes'){
    skipAggregation = 'yes'
}

skipTertiary = 'no'
if ( params.containsKey('skipTertiary') && params.skipTertiary == 'yes'){
    skipTertiary = 'yes'
}

// Now we pick the SDRFs to use

GIT_SDRFS = Channel.fromPath("$SCXA_WORKFLOW_ROOT/metadata/**/*.sdrf.txt", checkIfExists: true).map{ f -> tuple("${f.simpleName}", f) }

// If user has supplied an experiment ID, then filter SDRFs to just that
// experiment. Otherwise set a blank filter to include all results

doneSuffix=''
sdrfFilter = ''
if ( params.containsKey('expName')){
    sdrfFilter="${params.expName}"
    doneSuffix=".${params.expName}"
}

// Do a join on the SDRF channels to get single tuple per experiment ID, then
// pick the first. This will uniquify by experiment ID and prioritise the Git
// SDRFs

SDRF = GIT_SDRFS
    .filter( ~/.*\/.*${sdrfFilter}\..*/ )
    .map{ r -> tuple(r[0], r[1]) }

// Determine which SDRF files have up-to-date bundles. For new/ updated
// experiments, output to the channel with the relevant IDF.
// We don't want this process cached, since we need it to return a different
// thing with the same inputs, depending on whether an existing bunle is found.
// Note that we also exclude any experiments placed manually in the excluded.txt
// file.

process find_new_updated {

    conda "${baseDir}/envs/atlas-fastq-provider.yml"
    
    executor 'local'

    cache false
        
    input:
        set val(expName), file(sdrfFile) from SDRF

    output:
        set val(expName), file(sdrfFile), file("${sdrfFile.getSimpleName()}.idf.txt") optional true into SDRF_IDF
        file('bundleLines.txt') optional true into OLD_BUNDLE_LINES

    """
        # Have we excluded this study?

        set +e
        grep -P "$expName\\t" $SCXA_RESULTS/excluded.txt > /dev/null
        excludeStatus=\$?
        set -e

        if [ \$excludeStatus -eq 0 ]; then
            exit 0
        fi

        # Start by assuming a new experiment
        
        newExperiment=1
        bundleManifests=\$(ls \$SCXA_RESULTS/$expName/*/bundle/MANIFEST 2>/dev/null || true)
        
        # If there are existing bundles for this experiment then just use
        # those, unless the related SDRFs have been updated
 
        if [ -n "\$bundleManifests" ] && [ "$overwrite" != 'yes' ]; then
            newExperiment=0
            while read -r bundleManifest; do
                if [ $sdrfFile -nt "\$bundleManifest" ]; then
                    newExperiment=1
                    rm -f bundleLines.txt
                    break
                else
                    species=\$(echo \$bundleManifest | awk -F'/' '{print \$(NF-2)}' | tr -d \'\\n\')
                    echo -e "$expName\\t\$species\\t$SCXA_RESULTS/$expName/\$species/bundle" > bundleLines.txt
                fi
            done <<< "\$(echo -e "\$bundleManifests")"
        fi

        # Check for the existence of Atlas prod lock files indicating loading
        # in progress.  Only tag new experiments when there are no loading
        # locks present
        
        bundleLocks=\$(ls \$SCXA_RESULTS/$expName/*/bundle/atlas_prod.loading 2>/dev/null || true)
      
        if [ \$newExperiment -eq 1 ] && [ -z "\$bundleLocks" ]; then
            cp \$(dirname \$(readlink $sdrfFile))/${expName}.idf.txt .
        fi
    """
}

// Derive the config files. We cache based on content (deeply). If we do run,
// we remove the downstream stored results, triggering the sub-workflows (not
// normally re-run). These will then check their own caches and re-run where
// required. 

process generate_config {

    cache 'deep'
    
    publishDir "$SCXA_CONF/study", mode: 'copy', overwrite: true
    
    errorStrategy 'ignore'

    conda 'r-optparse r-data.table r-workflowscriptscommon'

    input:
        set val(expName), file(sdrfFile), file(idfFile) from SDRF_IDF.take(params.numExpsProcessedAtOnce)

    output:
        file('*.conf') optional true into CONF_FILES
        file('*.sdrf.txt') optional true into SDRF_FILES        
        file('*.metadata.tsv') optional true into METADATA_FILES       
        file('bundleLines.txt') optional true into NORERUN_BUNDLE_LINES

    """
    
    sdrfToNfConf.R \
        --sdrf=\$(readlink $sdrfFile) \
        --idf=$idfFile \
        --name=$expName \
        --verbose \
        --out_conf try_conf

    # If there's already a MANIFEST then this is a potential re-run. Only do
    # this if the config changes impact on analysis. This will be the case if
    # there's a difference between new derived files and the old ones.

    if [ -e $SCXA_RESULTS/$expName/*/bundle/MANIFEST ] && [ "$overwrite" != 'yes' ]; then
    
        # Start by assuming an old experiment
        
        newExperiment=0

        # Check for existing copies the config and derived SDRF. If neither have
        # changed, no need to re-analyse- even if the SDRF has been edited.

        for ext in sdrf.txt conf; do
            while read -r tc; do
                if [ -e \$SCXA_CONF/study/\$(basename \$tc) ]; then
                    set +e
                    diff -I '^//' \$tc \$SCXA_CONF/study/\$(basename \$tc) > /dev/null 2>&1

                    if [ \$? -ne 0 ]; then
                        echo \$tc is different to \$SCXA_CONF/study/\$(basename \$tc)
                        newExperiment=1
                        break
                    fi

                    set -e
                else
                    newExperiment=1
                    break
                fi
            done <<< "\$(ls try_conf/*.\$ext)"
        done
    else
        newExperiment=1
    fi

    # If this is a new experiment, or one where config changes impact on
    # analysis, re-run

   if [ \$newExperiment -eq 1 ]; then

        mv try_conf/* .

        # Only remove downstream results where we're not re-using them
        reset_stages='bundle'
        if [ "$skipQuantification" == 'no' ]; then
            reset_stages="quantification aggregation scanpy \$reset_stages"
        elif [ "$skipAggregation" = 'no' ]; then 
            reset_stages="aggregation scanpy \$reset_stages"
        elif [ "$skipTertiary" = 'no' ]; then 
            reset_stages="scanpy \$reset_stages"
        fi

        for stage in \$reset_stages; do
            rm -rf $SCXA_RESULTS/$expName/*/\$stage
        done
    else
        echo "Config files have not changed, no need to re-run"

        # Update the manifest time stamps to prevent re-config next time

        touch -m $SCXA_RESULTS/$expName/*/bundle/MANIFEST
                    
        for species in \$(ls $SCXA_RESULTS/$expName); do
            echo -e "$expName\\t\$species\\t$SCXA_RESULTS/$expName/\$species/bundle" > bundleLines.txt
        done

    fi
    rm -rf try_conf 
    """
}

// Flatten configs (each experiment could return mutliple). 

CONF_FILES
    .flatten()
    .set{FLAT_CONF_FILES}

SDRF_FILES
    .flatten()
    .set{FLAT_SDRF_FILES}

METADATA_FILES
    .flatten()
    .set{FLAT_METADATA_FILES}

// Tag SDRF files with file root

process tag_sdrf{
    
    input:
        file(sdrfFile) from FLAT_SDRF_FILES

    output:
        set stdout, file(sdrfFile) into TAGGED_SDRF_FILES

    """
        echo -n $sdrfFile | sed 's/.sdrf.txt//'
    """
}

// Re-tag conf files with exp name and species

process tag_conf{
    
    input:
        file(confFile) from FLAT_CONF_FILES

    output:
        set stdout, file(confFile) into TAGGED_CONF_FILES

    """
        echo -n $confFile | sed 's/.conf//'
    """
}

// Re-tag metadata files with exp name and species

process tag_meta{
    
    input:
        file(metaFile) from FLAT_METADATA_FILES

    output:
        set stdout, file(metaFile) into TAGGED_METADATA_FILES

    """
        echo -n $metaFile | sed 's/.metadata.tsv//'
    """
}

TAGGED_CONF_FILES
    .join( TAGGED_SDRF_FILES ) 
    .join( TAGGED_METADATA_FILES ) 
    .set{ CONF_SDRF }

// Mark up files with metadata for grouping

process markup_conf_files {

    executor 'local'
    
    cache 'deep'
    
    conda 'pyyaml' 
    
    input:
        set val(tag), file(confFile), file(sdrfFile), file(metaFile) from CONF_SDRF

    output:
        set file('expName'), file('species'), file('protocol'), file("out/$confFile"), file("out/$sdrfFile"), file("out/$metaFile") into MARKUP_CONF_FILES

    """
        parseNfConfig.py --paramFile $confFile --paramKeys params,name > expName
        parseNfConfig.py --paramFile $confFile --paramKeys params,organism > species
        protocol=\$(parseNfConfig.py --paramFile $confFile --paramKeys params,protocol)
        echo -n \$protocol > protocol

        mkdir -p out
        cp -p $confFile out/$confFile    
        cp -p $sdrfFile out/$sdrfFile    
        cp -p $metaFile out/$metaFile    
        
        protocolConfig=${baseDir}/conf/protocol/\${protocol}.conf
        if [ ! -e \$protocolConfig ]; then
            echo "\$protocolConfig does not exist" 1>&2
            exit 1
        fi 

        echo "includeConfig '${baseDir}/params.config'" >> out/${confFile}.tmp
        echo "includeConfig '\$protocolConfig'" >> out/${confFile}.tmp
        cat $confFile >> out/${confFile}.tmp
        mv out/${confFile}.tmp out/${confFile}
    """
}

MARKUP_CONF_FILES
    .map{ row-> tuple( row[0].text, row[1].text, row[2].text, row[3], row[4], row[5]) }        
    .set{ CONF_BY_META }

CONF_BY_META
    .into{
        CONF_BY_META_FOR_REFERENCE
        CONF_BY_META_FOR_QUANTIFY
        CONF_BY_META_FOR_AGGREGATION
        CONF_BY_META_FOR_SCANPY
    }

// Locate reference fines

process add_reference {

    conda 'pyyaml' 
    
    cache 'deep'
    
    publishDir "$SCXA_RESULTS/$expName/$species/reference", mode: 'copy', overwrite: true

    errorStrategy { task.attempt<=3 ? 'retry' : 'finish' }

    input:
        set val(expName), val(species), val(protocol), file(confFile), file(sdrfFile), file(metaFile) from CONF_BY_META_FOR_REFERENCE
    
    output:
        set val(expName), val(species), val(protocol), file(confFile), file(sdrfFile), file(metaFile), file("*.fa.gz"), file("*.gtf.gz"), stdout into CONF_WITH_REFERENCE

    """
    # Use references from the ISL setup

    if [ -n "\$ISL_GENOMES" ] && [ "\$IRAP_CONFIG_DIR" != '' ] && [ "\$IRAP_DATA" != '' ]; then

        irap_species_conf=$IRAP_CONFIG_DIR/${species}.conf

        if [ ${params.islReferenceType} = 'newest' ]; then
            gtf_pattern=\$(basename \$(cat \$ISL_GENOMES | grep $species | awk '{print \$6}') | sed 's/RELNO/\\*/')
            cdna_pattern=\$(basename \$(cat \$ISL_GENOMES | grep $species | awk '{print \$5}') | sed 's/RELNO/\\*/' | sed 's/.fa.gz/.\\*.fa.gz/') 

            cdna_gtf=\$(ls \$IRAP_DATA/reference/${species}/\$gtf_pattern | sort -rV | head -n 1)
            cdna_fasta=\$(ls \$IRAP_DATA/reference/${species}/\$cdna_pattern | sort -rV | head -n 1)
        else
            cdna_fasta=$IRAP_DATA/reference/$species/\$(parseIslConfig.sh \$irap_species_conf cdna_file)   
            cdna_gtf=$IRAP_DATA/reference/$species/\$(parseIslConfig.sh \$irap_species_conf gtf_file)   
        fi

        if [ ! -e "\$cdna_fasta" ]; then
            echo "Fasta file \$cdna_fasta does not exist" 1>&2
            exit 1
        elif [ ! -e "\$cdna_gtf" ]; then
            echo "GTF file \$cdna_gtf does not exist" 1>&2
            exit 1
        fi

        contamination_index=\$(parseIslConfig.sh \$irap_species_conf cont_index)  
    else
        echo "All of environment variables ISL_GENOMES, IRAP_CONFIG_DIR AND IRAP_DATA must be set" 1>&2
        exit 1
    fi
   
    echo -n "\$contamination_index"
    ln -s \$cdna_fasta
    ln -s \$cdna_gtf
    """
}

CONF_WITH_REFERENCE
    .into{
        CONF_WITH_ORIG_REFERENCE_FOR_PREPARE
        CONF_WITH_ORIG_REFERENCE_FOR_TERTIARY
    }

// Prepare a reference depending on spikes

process prepare_reference {

    conda 'pyyaml' 
    
    cache 'deep'
    
    publishDir "$SCXA_RESULTS/$expName/$species/reference", mode: 'copy', overwrite: true

    errorStrategy { task.attempt<=3 ? 'retry' : 'finish' }

    input:
        set val(expName), val(species), val(protocol), file(confFile), file(sdrfFile), file(metaFile), file(referenceFasta), file(referenceGtf), val(contaminationIndex) from CONF_WITH_ORIG_REFERENCE_FOR_PREPARE
    
    output:
        set val(expName), val(species), val(protocol), file(confFile), file(sdrfFile), file(metaFile), file("out/*.fa.gz"), file("out/*.gtf.gz"), val(contaminationIndex) into CONF_WITH_PREPARED_REFERENCE

    """
    mkdir -p out
    spikes=\$(parseNfConfig.py --paramFile $confFile --paramKeys params,spikes)

    if [ \$spikes != 'None' ]; then
        spikes_conf="$SCXA_PRE_CONF/reference/\${spikes}.conf"
        spikes_fasta=$SCXA_DATA/reference/\$(parseNfConfig.py --paramFile \$spikes_conf --paramKeys params,reference,spikes,cdna)
        spikes_gtf=$SCXA_DATA/reference/\$(parseNfConfig.py --paramFile \$spikes_conf --paramKeys params,reference,spikes,gtf)
        
        cat $referenceFasta \$spikes_fasta > out/reference.fa.gz
        cat $referenceGtf \$spikes_gtf > out/reference.gtf.gz
    else
        cp -P $referenceGtf out/reference.gtf.gz
        cp -P $referenceFasta out/reference.fa.gz
    fi
    """
}

// Synchronise the GTF and the FASTA 

process synchronize_ref_files {

    conda "${baseDir}/envs/atlas-gene-annotation-manipulation.yml"
    
    cache 'deep'

    memory { 5.GB * task.attempt }

    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3  ? 'retry' : 'ignore' }
    maxRetries 3
        
    input:
        set val(expName), val(species), val(protocol), file(confFile), file(sdrfFile), file(referenceFasta), file(referenceGtf), val(contaminationIndex) from CONF_WITH_PREPARED_REFERENCE

    output:
        set val(expName), val(species), val(protocol), file(confFile), file(sdrfFile), file('cleanedCdna.fa.gz'), file(referenceGtf), val(contaminationIndex), file('transcript_to_gene.txt') into CONF_WITH_SYNCH_REFERENCE

    """
    gtf2featureAnnotation.R --gtf-file ${referenceGtf} --no-header --version-transcripts --filter-cdnas ${referenceFasta} \
        --filter-cdnas-field "transcript_id" --filter-cdnas-output cleanedCdna.fa.gz --feature-type "transcript" \
        --first-field "transcript_id" --output-file transcript_to_gene.txt --fields "transcript_id,gene_id"    
    """
}

CONF_WITH_SYNCH_REFERENCE
    .into{
        CONF_FOR_QUANT
        CONF_FOR_AGGR
    }

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


// If we just want to run tertiary using our published quantification results
// from previous runs, we can do that

if ( skipQuantification == 'yes'){

    process spoof_quantify {
        
        executor 'local'

        input:
            set val(expName), val(species), val(protocol), file(confFile), file(sdrfFile), file(metaFile), file(referenceFasta), file(referenceGtf), val(contaminationIndex) from CONF_FOR_QUANT

        output:
            set val(expName), val(species), val(protocol), file("results/*") into QUANT_RESULTS

        """
            mkdir -p results
            for PROG in alevin kallisto; do
                if [ -e $SCXA_RESULTS/$expName/$species/quantification/$protocol/\$PROG ]; then
                    ln -s $SCXA_RESULTS/$expName/$species/quantification/$protocol/\$PROG results/\$PROG
                fi
            done
        """
    }

}else{

    // Separate droplet and smart-type experiments

    DROPLET_CONF = Channel.create()
    SMART_CONF = Channel.create()

    CONF_FOR_QUANT.choice( SMART_CONF, DROPLET_CONF ) {a -> 
        dropletProtocols.contains(a[2]) ? 1 : 0
    }

    // Run quantification with https://github.com/ebi-gene-expression-group/scxa-smartseq-quantification-workflow

    process smart_quantify {

        maxForks params.maxConcurrentQuantifications

        conda "${baseDir}/envs/nextflow.yml"
        
        cache 'deep'

        publishDir "$SCXA_RESULTS/$expName/$species/quantification/$protocol", mode: 'copy', overwrite: true
        
        memory { 10.GB * task.attempt }
        errorStrategy { task.attempt<=10 ? 'retry' : 'finish' }
        maxRetries 10
        
        input:
            set val(expName), val(species), val(protocol), file(confFile), file(sdrfFile), file(metaFile), file(referenceFasta), file(referenceGtf), val(contaminationIndex) from SMART_CONF
            val flag from INIT_DONE_SMART

        output:
            set val(expName), val(species), val(protocol), file ("kallisto") into SMART_KALLISTO_DIRS 
            set val(expName), val(species), val(protocol), file ("qc") into SMART_QUANT_QC
            file('quantification.log')    

        """
            for stage in aggregation scanpy bundle; do
                rm -rf $SCXA_RESULTS/$expName/$species/\$stage/$protocol
            done

            RESULTS_ROOT=\$PWD
            SUBDIR="$expName/$species/quantification/$protocol"     

            mkdir -p $SCXA_WORK/\$SUBDIR
            mkdir -p $SCXA_NEXTFLOW/\$SUBDIR
            mkdir -p $SCXA_RESULTS/\$SUBDIR/reports
            pushd $SCXA_NEXTFLOW/\$SUBDIR > /dev/null

            CONT_INDEX=''
            if [ "$contaminationIndex" != 'None' ]; then
                CONT_INDEX="--contaminationIndex $contaminationIndex"
            fi 

            nextflow run \
                -config \$RESULTS_ROOT/$confFile \
                --sdrf \$RESULTS_ROOT/$sdrfFile \
                --referenceFasta \$RESULTS_ROOT/$referenceFasta \$CONT_INDEX \
                --resultsRoot \$RESULTS_ROOT \
                --enaSshUser $enaSshUser \
                --manualDownloadFolder $SCXA_DATA/ManuallyDownloaded/$expName \
                -resume \
                $SCXA_WORKFLOW_ROOT/workflow/scxa-workflows/w_smart-seq_quantification/main.nf \
                -work-dir $SCXA_WORK/\$SUBDIR \
                -with-report $SCXA_RESULTS/\$SUBDIR/reports/report.html \
                -with-trace  $SCXA_RESULTS/\$SUBDIR/reports/trace.txt \
                -N $SCXA_REPORT_EMAIL \
                -with-dag $SCXA_RESULTS/\$SUBDIR/reports/flowchart.pdf

            if [ \$? -ne 0 ]; then
                echo "Workflow failed for $expName - $species - scxa-smartseq-quantification-workflow" 1>&2
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

        publishDir "$SCXA_RESULTS/$expName/$species/quantification/$protocol", mode: 'copy', overwrite: true
        
        memory { 10.GB * task.attempt }
        errorStrategy { task.attempt<=5 ? 'retry' : 'finish' }

        input:
            set val(expName), val(species), val(protocol), file(confFile), file(sdrfFile), file(metaFile), file(referenceFasta), file(referenceGtf), val(contaminationIndex) from DROPLET_CONF
            val flag from INIT_DONE_DROPLET

        output:
            set val(expName), val(species), val(protocol), file("alevin") into ALEVIN_DROPLET_DIRS
            file('quantification.log')    

        """
            for stage in aggregation scanpy bundle; do
                rm -rf $SCXA_RESULTS/$expName/$species/\$stage/$protocol
            done

            RESULTS_ROOT=\$PWD
            SUBDIR="$expName/$species/quantification/$protocol"     

            mkdir -p $SCXA_WORK/\$SUBDIR
            mkdir -p $SCXA_NEXTFLOW/\$SUBDIR
            mkdir -p $SCXA_RESULTS/\$SUBDIR/reports
            pushd $SCXA_NEXTFLOW/\$SUBDIR > /dev/null

            nextflow run \
                -config \$RESULTS_ROOT/$confFile \
                --resultsRoot \$RESULTS_ROOT \
                --sdrf \$RESULTS_ROOT/$sdrfFile \
                --referenceFasta \$RESULTS_ROOT/$referenceFasta \
                --transcriptToGene \$RESULTS_ROOT/$transcriptToGene \
                --protocol $protocol \
                --manualDownloadFolder $SCXA_DATA/ManuallyDownloaded/$expName \
                -resume \
                $SCXA_WORKFLOW_ROOT/workflow/scxa-workflows/w_droplet_quantification/main.nf \
                -work-dir $SCXA_WORK/\$SUBDIR \
                -with-report $SCXA_RESULTS/\$SUBDIR/reports/report.html \
                -with-trace  $SCXA_RESULTS/\$SUBDIR/reports/trace.txt \
                -N $SCXA_REPORT_EMAIL \
                -with-dag $SCXA_RESULTS/\$SUBDIR/reports/flowchart.pdf

            if [ \$? -ne 0 ]; then
                echo "Workflow failed for $expName - $species - scxa-droplet-quantification-workflow" 1>&2
                exit 1
            fi
            
            popd > /dev/null

            cp $SCXA_NEXTFLOW/\$SUBDIR/.nextflow.log quantification.log
        """
    }

    // Collect the smart and droplet workflow results

    SMART_KALLISTO_DIRS
        .concat(ALEVIN_DROPLET_DIRS)
        .set { QUANT_RESULTS }
}


// Fetch the config to use in aggregation

CONF_FOR_AGGR
    .join ( QUANT_RESULTS, by: [0,1,2] )
    .groupTuple( by: [0,1] )
    .set{ GROUPED_QUANTIFICATION_RESULTS }

if (skipAggregation == 'yes' ){

    process spoof_aggregate {
    
        executor 'local'
        
        input:
            set val(expName), val(species), file('quant_results/??/protocol'), file('quant_results/??/*'), file('quant_results/??/*'), file('quant_results/??/*'), file('quant_results/??/*'), file('quant_results/??/*'), file('quant_results/??/*'), file('quant_results/??/*') from GROUPED_QUANTIFICATION_RESULTS   

        output:
            set val(expName), val(species), file("matrices/counts_mtx.zip") into COUNT_MATRICES
            set val(expName), val(species), file("matrices/tpm_mtx.zip") optional true into TPM_MATRICES
            set val(expName), val(species), file("matrices/stats.tsv") optional true into KALLISTO_STATS
    
        """
            ln -s $SCXA_RESULTS/$expName/$species/aggregation/matrices 
        """
    }

}else{

    // Run aggregation with https://github.com/ebi-gene-expression-group/scxa-aggregation-workflow

    process aggregate {
        
        conda "${baseDir}/envs/nextflow.yml"
        
        maxForks 6

        publishDir "$SCXA_RESULTS/$expName/$species/aggregation", mode: 'copy', overwrite: true
        
        memory { 4.GB * task.attempt }
        errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3 ? 'retry' : 'finish' }
        maxRetries 20
        
        input:
            set val(expName), val(species), file('quant_results/??/protocol'), file('quant_results/??/*'), file('quant_results/??/*'), file('quant_results/??/*'), file('quant_results/??/*'), file('quant_results/??/*'), file('quant_results/??/*'), file('quant_results/??/*') from GROUPED_QUANTIFICATION_RESULTS   

        output:
            set val(expName), val(species), file("matrices/counts_mtx.zip") into COUNT_MATRICES
            set val(expName), val(species), file("matrices/tpm_mtx.zip") optional true into TPM_MATRICES
            set val(expName), val(species), file("matrices/kallisto_stats.tsv") optional true into KALLISTO_STATS
            set val(expName), val(species), file("matrices/alevin_stats.tsv") optional true into ALEVIN_STATS

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

            # If we have a species-wise config, supply to aggregation

            species_conf=$SCXA_PRE_CONF/reference/${species}.conf
            opt_conf=
            if  [ -e \$species_conf ]; then
                opt_conf=" --config \$species_conf"
            fi            

            nextflow run \$opt_conf \
                --resultsRoot \$RESULTS_ROOT \
                --quantDir \$RESULTS_ROOT/quant_results \
                -resume \
                $SCXA_WORKFLOW_ROOT/workflow/scxa-workflows/w_aggregation/main.nf \
                -work-dir $SCXA_WORK/\$SUBDIR \
                -with-report $SCXA_RESULTS/\$SUBDIR/reports/report.html \
                -with-trace  $SCXA_RESULTS/\$SUBDIR/reports/trace.txt \
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
}

// Run Scanpy with https://github.com/ebi-gene-expression-group/scanpy-workflow
    
// Take the original reference GTF (before spikes or anything else were added),
// for use in Scanpy as a filter. We only need one of the potentially multiple
// per experiment (where multiple protocols are present), add we then need to
// connect to the count matrix outputs- hence the incantations below.
//
// Note that we record the worfklow type (smartseq or droplet) for use in the
// reporting in bundling. We assume that we don't mix SMART and droplet within
// an experiment

CONF_WITH_ORIG_REFERENCE_FOR_TERTIARY
    .groupTuple( by: [0,1] )
    .map{ row-> tuple( row[0], row[1], row[2].join(","), row[3][0], row[5][0], row[6][0]) }
    .unique()
    .join(COUNT_MATRICES, by: [0,1])
    .set { TERTIARY_INPUTS }         

if ( tertiaryWorkflow == 'scanpy-workflow'){

    process scanpy {
        
        conda "${baseDir}/envs/nextflow.yml"

        publishDir "$SCXA_RESULTS/$expName/$species/scanpy", mode: 'copy', overwrite: true
        
        memory { 4.GB * task.attempt }
        errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }
        maxRetries 20
          
        input:
            set val(expName), val(species), val(protocolList), file(confFile), file(metaFile), file(referenceGtf), file(countMatrix) from TERTIARY_INPUTS

        output:
            set val(expName), val(species), val(protocolList), file(confFile), file(referenceFasta), file(referenceGtf), file("matrices/${countMatrix}"), file("matrices/*_filter_cells_genes.zip"), file("matrices/*_normalised.zip"), file("clustering/clusters.txt"), file("umap"), file("tsne"), file("markers") into TERTIARY_RESULTS
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
                -with-trace  $SCXA_RESULTS/\$SUBDIR/reports/trace.txt \
                -N $SCXA_REPORT_EMAIL \
                -with-dag $SCXA_RESULTS/\$SUBDIR/reports/flowchart.pdf

            if [ \$? -ne 0 ]; then
                echo "Workflow failed for $expName - $species - scanpy-workflow" 1>&2
                exit 1
            fi
            
            popd > /dev/null

            cp -P ${countMatrix} matrices

            cp $SCXA_NEXTFLOW/\$SUBDIR/.nextflow.log scanpy.log
       """
        
    }
    
}else if ( tertiaryWorkflow == 'scanpy-galaxy' ) {

    if (skipTertiary == 'yes' ){

        process spoof_scanpy_galaxy {
        
            executor 'local'
            input:
                set val(expName), val(species), val(protocolList), file(confFile), file(metaFile), file(referenceGtf), file(countMatrix) from TERTIARY_INPUTS
            
            output:
                set val(expName), val(species), val(protocolList), file(confFile), file(referenceFasta), file(referenceGtf), file("matrices/${countMatrix}"), file("matrices/raw_filtered.zip"), file("matrices/filtered_normalised.zip"), file("clusters_for_bundle.txt"), file("umap"), file("tsne"), file("markers"), file('clustering_software_versions.txt') into TERTIARY_RESULTS
        
            """
                ln -s $SCXA_RESULTS/$expName/$species/scanpy/matrices 
                ln -s $SCXA_RESULTS/$expName/$species/scanpy/clusters_for_bundle.txt
                ln -s $SCXA_RESULTS/$expName/$species/scanpy/umap 
                ln -s $SCXA_RESULTS/$expName/$species/scanpy/tsne
                ln -s $SCXA_RESULTS/$expName/$species/scanpy/markers
                ln -s $SCXA_RESULTS/$expName/$species/scanpy/clustering_software_versions.txt
            """
        }

    }else{

        process scanpy_galaxy {
            
            // Exit status of 3 is just Galaxy being annoying with history
            // deletion, no cause to error

            validExitStatus 0,3
        
            maxForks params.maxConcurrentScanpyGalaxy

            conda "${baseDir}/envs/bioblend.yml"

            publishDir "$SCXA_RESULTS/$expName/$species/scanpy", mode: 'copy', overwrite: true
            
            memory { 4.GB * task.attempt }
            errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 || task.attempt<=3  ? 'retry' : 'ignore' }
            maxRetries 3
              
            input:
                set val(expName), val(species), val(protocolList), file(confFile), file(metaFile), file(referenceGtf), file(countMatrix) from TERTIARY_INPUTS

            output:
                set val(expName), val(species), val(protocolList), file(confFile), file(referenceFasta), file(referenceGtf), file("matrices/${countMatrix}"), file("matrices/raw_filtered.zip"), file("matrices/filtered_normalised.zip"), file("clusters_for_bundle.txt"), file("umap"), file("tsne"), file("markers"), file('clustering_software_versions.txt') into TERTIARY_RESULTS

            script: 

                def isDroplet='False'
                protocol_list=protocolList.split(',').toList()
                experiment_droplet_protocols=protocol_list.intersect(dropletProtocols)
                if ( experiment_droplet_protocols.size() > 0){
                    isDroplet='True'
                }

            """
                rm -rf $SCXA_RESULTS/$expName/$species/bundle

                # Galaxy workflow wants the matrix components separately

                zipdir=\$(unzip -qql ${countMatrix.getBaseName()}.zip | head -n1 | tr -s ' ' | cut -d' ' -f5- | sed 's|/||')
                unzip ${countMatrix.getBaseName()}        
                gzip \${zipdir}/matrix.mtx
                gzip \${zipdir}/genes.tsv
                gzip \${zipdir}/barcodes.tsv

                export species=$species
                export expName=$expName
                export gtf_file=$referenceGtf

                export matrix_file=\${zipdir}/matrix.mtx.gz
                export genes_file=\${zipdir}/genes.tsv.gz
                export barcodes_file=\${zipdir}/barcodes.tsv.gz
                export tpm_filtering='False'

                export create_conda_env=no
                export GALAXY_CRED_FILE=$galaxyCredentials
                export GALAXY_INSTANCE=ebi_cluster_ge_user
       
                if [ "$isDroplet" = 'True' ]; then
                    export FLAVOUR=w_droplet_clustering
                else
                    export FLAVOUR=w_smart-seq_clustering
                fi

                run_flavour_workflows.sh

                if [ \$? -eq 0 ]; then
                    mkdir -p matrices
                    cp -P ${countMatrix} matrices
                    
                    # Group associated matrix files

                    for matrix_type in raw_filtered filtered_normalised; do
                        mkdir -p matrices/\${matrix_type} 
                        
                        for file in matrix.mtx genes.tsv barcodes.tsv; do
                            if [ ! -e \${matrix_type}_\${file} ]; then
                                echo "\${matrix_type}_\${file} does not exist" 1>&2
                                exit 2
                            else
                                mv \${matrix_type}_\${file} matrices/\${matrix_type}/\${file}
                            fi
                        done

                        pushd matrices > /dev/null 
                        zip -r \${matrix_type}.zip \${matrix_type}
                        popd > /dev/null 
                    done

                    # Organise other outputs
                 
                    mkdir -p tsne && mv tsne_* tsne 
                    mkdir -p umap && mv umap_* umap

                    mkdir -p markers
                    set +e
                    marker_files=\$(ls markers_* 2>/dev/null | grep -v markers_clusters_resolution)
                    if [ \$? -ne 0 ]; then
                        echo "No marker files present"
                        touch markers/NOMARKERS
                    else
                        mv \$marker_files markers
                    fi
                fi 
           """
        }
    }

} else{
    
    process spoof_tertiary {
    
        executor 'local'
        
        input:
            set val(expName), val(species), val(protocolList), file(confFile), file(referenceFasta), file(referenceGtf), file(countMatrix) from TERTIARY_INPUTS

        output:
            set val(expName), val(species), val(protocolList), file(confFile), file(referenceFasta), file(referenceGtf), file("matrices/${countMatrix}"), file("NOFILT"), file("NONORM"), file("NOCLUST"), file("NOUMAP"), file("NOTSNE"), file("NOMARKERS"), file('NOSOFTWARE') into TERTIARY_RESULTS

        """
            mkdir -p matrices
            cp -p ${countMatrix} matrices 
            touch NOFILT NONORM NOPCA NOCLUST NOUMAP NOTSNE NOMARKERS NOSOFTWARE
        """
    }    
        
}

// Make a bundle from the outputs. We may or may not have a TPM matrix,
// depending on protocol

TERTIARY_RESULTS
    .join(TPM_MATRICES, remainder: true, by: [0,1])
    .filter { it[2] != null }
    .set { BUNDLE_INPUTS } 

process bundle {
    
    conda "${baseDir}/envs/nextflow.yml"

    publishDir "$SCXA_RESULTS/$expName/$species", mode: 'copy', overwrite: true
    
    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }
    maxRetries 20
    
    input:
        set val(expName), val(species), val(protocolList), file(confFile), file(referenceFasta), file(referenceGtf), file(rawMatrix), file(filteredMatrix), file(normalisedMatrix), file(clusters), file('*'), file('*'), file('*'), file(softwareReport), file(tpmMatrix) from BUNDLE_INPUTS
        
    output:
        file('bundle/*')
        file('bundle.log')
        file('bundleLines.txt') into NEW_BUNDLE_LINES
        
    """    
        RESULTS_ROOT=\$PWD
        SUBDIR="$expName/$species/bundle"     

        TPM_OPTIONS=''
        tpm_filesize=\$(stat --printf="%s" \$(readlink ${tpmMatrix}))
        if [ "$tpmMatrix" != 'null' ] && [ \$tpm_filesize -gt 0 ]; then
            TPM_OPTIONS="--tpmMatrix ${tpmMatrix}"
        fi

        TERTIARY_OPTIONS=''
        if [ "$tertiaryWorkflow" == 'scanpy-workflow' ] || [ "$tertiaryWorkflow" == 'scanpy-galaxy' ]; then
            TERTIARY_OPTIONS="--tertiaryWorkflow $tertiaryWorkflow --rawFilteredMatrix ${filteredMatrix} --normalisedMatrix ${normalisedMatrix} --clusters ${clusters} --tsneDir tsne --markersDir markers --tertiarySoftwareReport ${softwareReport}"
        fi 
        
        mkdir -p $SCXA_WORK/\$SUBDIR
        mkdir -p $SCXA_NEXTFLOW/\$SUBDIR
        mkdir -p $SCXA_RESULTS/\$SUBDIR/reports
        pushd $SCXA_NEXTFLOW/\$SUBDIR > /dev/null

        nextflow run \
            -config \$RESULTS_ROOT/$confFile \
            --masterWorkflow scxa-control-workflow \
            --resultsRoot \$RESULTS_ROOT \
            --protocolList ${protocolList} \
            --rawMatrix ${rawMatrix} \$TPM_OPTIONS \
            --referenceFasta $referenceFasta \
            --referenceGtf $referenceGtf \$TERTIARY_OPTIONS \
            -resume \
            $SCXA_WORKFLOW_ROOT/workflow/scxa-workflows/w_bundle/main.nf \
            -work-dir $SCXA_WORK/\$SUBDIR \
            -with-report $SCXA_RESULTS/\$SUBDIR/reports/report.html \
            -with-trace  $SCXA_RESULTS/\$SUBDIR/reports/trace.txt \
            -N $SCXA_REPORT_EMAIL \
            -with-dag $SCXA_RESULTS/\$SUBDIR/reports/flowchart.pdf

        if [ \$? -ne 0 ]; then
            echo "Workflow failed for $expName - $species - scxa-bundle-workflow" 1>&2
            exit 1
        fi
        
        popd > /dev/null

        cp $SCXA_NEXTFLOW/\$SUBDIR/.nextflow.log bundle.log

        chmod g+rwx $SCXA_RESULTS/$expName/$species/bundle
        echo -e "$expName\\t$species\\t$SCXA_RESULTS/$expName/$species/bundle" > bundleLines.txt
   """
    
}

// Record the completed bundles

OLD_BUNDLE_LINES
    .concat ( NORERUN_BUNDLE_LINES )
    .concat ( NEW_BUNDLE_LINES )
    .collectFile(name: 'these.all.done.txt', sort: true)
    .set{
        BUNDLE_LINES
    }

// For each experiment and species with a completed bundle, remove the work
// dir. This can take a little while for large experiments, so background the
// process so future runs are not delayed.

process cleanup {
    
    executor 'local'
    
    publishDir "$SCXA_RESULTS", mode: 'copy', overwrite: true
    
    input:
        file(bundleLines) from BUNDLE_LINES
    
    output:
        file("all.done${doneSuffix}.txt")

    """
    touch $SCXA_WORK/.success
    cat ${bundleLines} | while read -r l; do
        expName=\$(echo "\$l" | awk '{print \$1}')
        species=\$(echo "\$l" | awk '{print \$2}')
        nohup rm -rf $SCXA_WORK/\$expName/\$species &
    done
    find $SCXA_WORK/ -maxdepth 2 -type d -empty -delete

    cp $bundleLines all.done${doneSuffix}.txt
    """
}

