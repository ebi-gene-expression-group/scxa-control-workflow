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

GIT_IDFS = Channel.fromPath("$SCXA_WORKFLOW_ROOT/metadata/**/*.idf.txt", checkIfExists: true).map{ f -> tuple("${f.simpleName}", f) }

// If user has supplied an experiment ID, then filter SDRFs to just that
// experiment. Otherwise set a blank filter to include all results

doneSuffix=''
idfFilter = ''
if ( params.containsKey('expName')){
    idfFilter="${params.expName}"
    doneSuffix=".${params.expName}"
}

// Do a join on the SDRF channels to get single tuple per experiment ID, then
// pick the first. This will uniquify by experiment ID and prioritise the Git
// IDFs

IDF = GIT_IDFS
    .filter( ~/.*\/.*${idfFilter}\..*/ )
    .map{ r -> tuple(r[0], r[1]) }

// Working from the IDF location, derive SDRF and related files

process compile_metadata {
    
    input:
        set val(expName), file(idfFile) from IDF

    output:
        set val(expName), file(idfFile), file("*.sdrf.txt") into SDRF_IDF
        set val(expName), file("${expName}.cells.txt") optional true into CELLS 

    """
    compileExpMetadata.sh $expName $idfFile $overwrite 
    """
}

SDRF_IDF
    .join(CELLS, remainder: true)
    .into{
        COMPILED_METADATA_FOR_QUANT
        COMPILED_METADATA_FOR_TERTIARY
    }

// Check the update status of the experiment

process checkExperimentStatus {

    input:
        set val(expName), file(idfFile), file(sdrfFile), file(cellsFile) from COMPILED_METADATA_FOR_QUANT

    output:
        set val(expName), file(idfFile), file(sdrfFile), file(cellsFile), stdout into COMPILED_METADATA_WITH_STATUS
        
    """
    checkExperimentStatus.sh $expName $sdrfFile $overwrite
    """
}

// Pipe experiments to different channels depending on update status 

UPDATED_EXPERIMENTS = Channel.create()
NOT_UPDATED_EXPERIMENTS = Channel.create()

COMPILED_METADATA_WITH_STATUS.choice( UPDATED_EXPERIMENTS, NOT_UPDATED_EXPERIMENTS ) {a -> 
    a[4] == 'old' ? 1 : 0
}

process get_not_updated_bundles {

    input:
        val expName from NOT_UPDATED_EXPERIMENTS.map{r -> r[0]}

    output:
        stdout NOT_UPDATED_BUNDLES

    """
        ls \$SCXA_RESULTS/$expName/*/bundle/MANIFEST | while read -r l; do
            species=\$(echo \$l | awk -F'/' '{print \$(NF-2)}' | tr -d \'\\n\')
            echo -e "$expName\\t\$species\\t$SCXA_RESULTS/$expName/\$species/bundle"
        done

    """
}

// Generate the configuration file

process generate_config {

    cache 'deep'
    
    errorStrategy 'ignore'

    conda 'r-optparse r-data.table r-workflowscriptscommon'

    input:
        set val(expName), file(idfFile), file(sdrfFile), file(cellsFile), val(expStatus) from UPDATED_EXPERIMENTS.take(params.numExpsProcessedAtOnce)

    output:
        set val(expName), file('*.conf') into DERIVED_CONF_FILES
        set val(expName), file('*.sdrf.txt') into DERIVED_SDRF_FILES        
        set val(expName), file('*.metadata.tsv') into DERIVED_METADATA_FILES       

    """
    sdrfToNfConf.R \
        --sdrf=\$(readlink $sdrfFile) \
        --idf=$idfFile \
        --name=$expName \
        --verbose \
        --out_conf \$(pwd)
    """
}

// The above will make tuples like [ 'E-MTAB-1234', ['foo.human.conf',
// 'foo.mouse.conf'] ]. We can use transpose to convert those to experiment/
// file pairs. 

DERIVED_CONF_FILES
    .transpose()
    .map{ row-> tuple( row[0], row[1].toString().split('\\.')[1], row[1] ) }
    .set{
        TAGGED_CONF
    }        

DERIVED_SDRF_FILES
    .transpose()
    .map{ row-> tuple( row[0], row[1].toString().split('\\.')[1], row[1] ) }
    .set{
        TAGGED_SDRF
    }        

DERIVED_METADATA_FILES
    .transpose()
    .map{ row-> tuple( row[0], row[1].toString().split('\\.')[1], row[1] ) }
    .set{
        TAGGED_META
    }        

// Pull further information from the config file for analysis

process annotate_config {

    executor 'local'
    
    cache 'deep'
    
    conda 'pyyaml' 

    input:
        set val(expName), val(species), file(confFile) from TAGGED_CONF

    output:
        set val(expName), val(species), stdout, file(confFile) into ANNOTATED_CONF

    """
        parseNfConfig.py --paramFile $confFile --paramKeys params,protocol
    """
} 


// Combine back to a set of files per experiment/species pair

ANNOTATED_CONF
    .join( TAGGED_SDRF, by: [0,1] ) 
    .join( TAGGED_META, by: [0,1] ) 
    .set{ COMPILED_DERIVED_CONFIG }

// Where analysis is already present, check that the config has actually changed
// in a way that impacts on analysis. This is distinct from
// checkExperimentStatus(), which just flags where the source SDRF has a newer
// timestamp than the analysis.

process check_experiment_changed{

    input:
        set val(expName), val(species), val(protocol), file(confFile), file(sdrfFile), file(metaFile) from COMPILED_DERIVED_CONFIG

    output:
        set val(expName), val(species), val(protocol), file(confFile), file(sdrfFile), file(metaFile), stdout into COMPILED_DERIVED_CONFIG_WITH_STATUS
        

    """
    checkExperimentChanges.sh $expName $species $confFile $sdrfFile 
    """
}

NEW_OR_CHANGED_EXPERIMENTS = Channel.create()
NOT_CHANGED_EXPERIMENTS = Channel.create()

COMPILED_DERIVED_CONFIG_WITH_STATUS.choice( NEW_OR_CHANGED_EXPERIMENTS, NOT_CHANGED_EXPERIMENTS ) {a -> 
    a[5] == 'unchanged' ? 1 : 0
}

// Generate bundle lines for the things updated but not actually changed for
// analysis

process get_not_changed_bundles {

    input:
        set val(expName), val(species) from NOT_CHANGED_EXPERIMENTS.map{r -> tuple(r[0], r[1])}

    output:
        stdout NOT_CHANGED_BUNDLES

    """
        # Update the manifest time stamps to prevent re-config next time
        touch -m $SCXA_RESULTS/$expName/$species/bundle/MANIFEST

        echo -e "$expName\\t$species\\t$SCXA_RESULTS/$expName/$species/bundle"
    """
}

// For experiments we know are due for re-analysis, remove pre-existing
// results (where appropriate). This is also the time to publish the config
// files, overwriting any already in place.

process reset_experiment{
    
    publishDir "$SCXA_CONF/study", mode: 'copy', overwrite: true

    input:
        set val(expName), val(species), val(protocol), file("in/*"), file("in/*"), file("in/*"), val(expStatus) from NEW_OR_CHANGED_EXPERIMENTS

    output:
        set val(expName), val(species), val(protocol), file("*.conf"), file("*.sdrf.txt"), file("*.metadata.tsv") into NEW_OR_RESET_EXPERIMENTS
        
    """
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
            rm -rf $SCXA_RESULTS/$expName/$species/\$stage
        done


        # Now we've removed the results, link the results for publishing

        cp -P in/*.sdrf.txt in/*conf in/*.metadata.tsv .
    """
}

// Locate reference fines

process add_reference {

    conda 'pyyaml' 
    
    cache 'deep'
    
    publishDir "$SCXA_RESULTS/$expName/$species/reference", mode: 'copy', overwrite: true

    errorStrategy { task.attempt<=3 ? 'retry' : 'finish' }

    input:
        set val(expName), val(species), val(protocol), file(confFile), file(sdrfFile) from NEW_OR_RESET_EXPERIMENTS
    
    output:
        set val(expName), val(species), val(protocol), file(confFile), file(sdrfFile), file("*.fa.gz"), file("*.gtf.gz"), stdout into CONF_WITH_REFERENCE

    """
    # Use references from an IRAP config
    species_conf=$SCXA_PRE_CONF/reference/${species}.conf
    if [ -e "\$species_conf" ]; then
        cdna_fasta=$SCXA_DATA/reference/\$(parseNfConfig.py --paramFile \$species_conf --paramKeys params,reference,cdna)
        cdna_gtf=$SCXA_DATA/reference/\$(parseNfConfig.py --paramFile \$species_conf --paramKeys params,reference,gtf)
        contamination_index=\$(parseNfConfig.py --paramFile \$species_conf --paramKeys params,reference,contamination_index)
        
        if [ \$contamination_index != 'None' ]; then
            contamination_index=$SCXA_DATA/contamination/\$contamination_index
        fi
    
    elif [ "\$IRAP_CONFIG_DIR" != '' ] && [ "\$IRAP_DATA" != '' ]; then
        irap_species_conf=$IRAP_CONFIG_DIR/${species}.conf
        cdna_fasta=$IRAP_DATA/reference/$species/\$(parseIslConfig.sh \$irap_species_conf cdna_file)   
        cdna_gtf=$IRAP_DATA/reference/$species/\$(parseIslConfig.sh \$irap_species_conf gtf_file)   
        contamination_index=\$(parseIslConfig.sh \$irap_species_conf cont_index)  
        
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
        set val(expName), val(species), val(protocol), file(confFile), file(sdrfFile), file(referenceFasta), file(referenceGtf), val(contaminationIndex) from CONF_WITH_ORIG_REFERENCE_FOR_PREPARE
    
    output:
        set val(expName), val(species), val(protocol), file(confFile), file(sdrfFile), file("out/*.fa.gz"), file("out/*.gtf.gz"), val(contaminationIndex) into CONF_WITH_PREPARED_REFERENCE

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
        cp -p $referenceGtf out/reference.gtf.gz
        cp -p $referenceFasta out/reference.fa.gz
    fi
    """
}

CONF_WITH_PREPARED_REFERENCE
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
            set val(expName), val(species), val(protocol), file(confFile), file(sdrfFile), file(referenceFasta), file(referenceGtf), val(contaminationIndex) from CONF_FOR_QUANT

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
            set val(expName), val(species), val(protocol), file(confFile), file(sdrfFile), file(referenceFasta), file(referenceGtf), val(contaminationIndex) from SMART_CONF
            val flag from INIT_DONE_SMART

        output:
            set val(expName), val(species), val(protocol), file ("kallisto") into SMART_KALLISTO_DIRS 
            set val(expName), val(species), val(protocol), file ("qc") into SMART_QUANT_QC
            file('quantification.log')    

        """
        submitQuantificationWorkflow.sh 'smart-seq' "$expName" "$species" "$protocol" "$confFile" "$sdrfFile" "$referenceFasta" "$referenceGtf" "$contaminationIndex" "$enaSshUser"
        """
    }

    // Run quantification with https://github.com/ebi-gene-expression-group/scxa-droplet-quantification-workflow

    process droplet_quantify {

        maxForks params.maxConcurrentQuantifications

        conda "${baseDir}/envs/nextflow.yml"

        cache 'deep'
        
        publishDir "$SCXA_RESULTS/$expName/$species/quantification/$protocol", mode: 'copy', overwrite: true
        
        memory { 10.GB * task.attempt }
        errorStrategy { task.attempt<=5 ? 'retry' : 'finish' }

        input:
            set val(expName), val(species), val(protocol), file(confFile), file(sdrfFile), file(referenceFasta), file(referenceGtf), val(contaminationIndex) from DROPLET_CONF
            val flag from INIT_DONE_DROPLET

        output:
            set val(expName), val(species), val(protocol), file("alevin") into ALEVIN_DROPLET_DIRS
            file('quantification.log')    

        """
        submitQuantificationWorkflow.sh 'droplet' "$expName" "$species" "$protocol" "$confFile" "$sdrfFile" "$referenceFasta" "$referenceGtf" "$contaminationIndex" "$enaSshUser"
        """
    }

    // Collect the smart and droplet workflow results

    SMART_KALLISTO_DIRS
        .concat(ALEVIN_DROPLET_DIRS
        .set { QUANT_RESULTS }
}

// Fetch the config to use in aggregation

CONF_FOR_AGGR
    .join ( QUANT_RESULTS, by: [0,1,2] )
    .groupTuple( by: [0,1] )
    .set{ GROUPED_QUANTIFICATION_RESULTS }

// Use existing aggregations if specified

if (skipAggregation == 'yes' ){

    process spoof_aggregate {
    
        executor 'local'
        
        input:
            set val(expName), val(species), file('quant_results/??/protocol'), file('quant_results/??/*'), file('quant_results/??/*'), file('quant_results/??/*'), file('quant_results/??/*'), file('quant_results/??/*'), file('quant_results/??/*') from GROUPED_QUANTIFICATION_RESULTS   

        output:
            set val(expName), val(species), file("matrices/counts_mtx.zip") into COUNT_MATRICES
            set val(expName), val(species), file("matrices/tpm_mtx.zip") optional true into TPM_MATRICES
            set val(expName), val(species), file("matrices/stats.tsv") optional true into KALLISTO_STATS
    
        """
            ln -s $SCXA_RESULTS/$expName/$species/aggregation/matrices 
        """
    }

}else{

    process aggregate {
        
        conda "${baseDir}/envs/nextflow.yml"
        
        maxForks 6

        publishDir "$SCXA_RESULTS/$expName/$species/aggregation", mode: 'copy', overwrite: true
        
        memory { 4.GB * task.attempt }
        errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3 ? 'retry' : 'finish' }
        maxRetries 20
        
        input:
            set val(expName), val(species), file('quant_results/??/protocol'), file('quant_results/??/*'), file('quant_results/??/*'), file('quant_results/??/*'), file('quant_results/??/*'), file('quant_results/??/*'), file('quant_results/??/*') from GROUPED_QUANTIFICATION_RESULTS   

        output:
            set val(expName), val(species), file("matrices/counts_mtx.zip") into COUNT_MATRICES
            set val(expName), val(species), file("matrices/tpm_mtx.zip") optional true into TPM_MATRICES
            set val(expName), val(species), file("matrices/kallisto_stats.tsv") optional true into KALLISTO_STATS
            set val(expName), val(species), file("matrices/alevin_stats.tsv") optional true into ALEVIN_STATS

        """
            submitAggregationWorkflow.sh "$expName" "$species"
        """
    }
}

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
    .map{ row-> tuple( row[0], row[1], row[2].join(","), row[7][0]) }
    .unique()
    .join(COUNT_MATRICES, by: [0,1])
    .set { PRE_TERTIARY_INPUTS }   

// Mark droplet experiments

process mark_droplet {
   
    input:
        set val(expName), val(species), val(protocolList), file(confFile), file(metaFile), file(referenceGtf), file(countMatrix) from PRE_TERTIARY_INPUTS

    output:
        set val(expName), val(species), val(protocolList), file(confFile), file(metaFile), file(referenceGtf), file(countMatrix), stdout into MARKED_TERTIARY_INPUTS
    
    script: 

        def isDroplet='False'
        protocol_list=protocolList.split(',').toList()
        experiment_droplet_protocols=protocol_list.intersect(dropletProtocols)
        if ( experiment_droplet_protocols.size() > 0){
            isDroplet='True'
        }

        """
        echo -n "$isDroplet"
        """

// Separate Droplet from non-droplet

DROPLET_PRE_TERIARY_INPUTS = Channel.create()
SMART_TERTIARY_INPUTS = Channel.create()

MARKED_TERTIARY_INPUTS.map{r -> r.add(file('NO_FILE')}.choice( DROPLET_TERIARY_INPUTS, SMART_TERTIARY_INPUTS ) {a -> 
    a[7] == 'True' ? 0 : 1
}

// Make a cell-run mapping, required by the SDRF condense process to 'explode'
// the SDRF annotations

process cell_run_mapping {
   
    input:
        set val(expName), val(species), val(protocolList), file(confFile), file(metaFile), file(referenceGtf), file(countMatrix), val(isDroplet), file(emptyFile) from DROPLET_PRE_TERTIARY_INPUTS
 
    output:
        set val(expName), val(species), val(protocolList), file(confFile), file(metaFile), file(referenceGtf), file(countMatrix), val(isDroplet), file('cell_to_library.txt') into DROPLET_TERTIARY_INPUTS
 
    """
    makeCellLibraryMapping.sh $countMatrix $confFile cell_to_library.txt 
    """
}

// Now make a condensed SDRF file. For this to operate correctly with droplet data, the cell-library mappings file must be present in the  


SMART_TERTIARY_INPUTS
    .concat(DROPLET_TERTIARY_INPUTS)
    .join(COMPILED_METADATA_FOR_TERTIARY, by: [0,1])
    .set{
        CONDENSE_INPUTS
    }


process condense_sdrf {
        
    conda "${baseDir}/envs/atlas-experiment-metadata.yml"

    input:
        set val(expName), val(species), val(protocolList), file(confFile), file(metaFile), file(referenceGtf), file(countMatrix), val(isDroplet), file(cell_to_lib), file(idfFile), file(origSdrfFile), file(cellsFile) from CONDENSE_INPUTS) 

    output:
       set val(expName), val(species), file("${expName}.${species}.condensed-sdrf.tsv") 
        

    """
    single_cell_condensed_sdrf.sh -e $expName -f $idfFile -o \$(pwd)

    # Re-add this: -z \$zoomaExclusions
    """        

}

//TAGGED_CONF.subscribe { println "value: $it" }

//process test{
    
//    input:
//        set val(expName), val(species), val(confFile) from TAGGED_CONF
//
//    """
//    echo $expName $species $confFile
//    """
//}

