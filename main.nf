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

// Combine back to a set of files per experiment/species pair

TAGGED_CONF
    .join( TAGGED_SDRF, by: [0,1] ) 
    .join( TAGGED_META, by: [0,1] ) 
    .set{ COMPILED_DERIVED_CONFIG }

// Where analysis is already present, check that the config has actually changed
// in a way that impacts on analysis. This is distinct from
// checkExperimentStatus(), which just flags where the source SDRF has a newer
// timestamp than the analysis.

process check_experiment_changed{

    input:
        set val(expName), val(species), file(confFile), file(sdrfFile), file(metaFile) from COMPILED_DERIVED_CONFIG

    output:
        set val(expName), val(species), file(confFile), file(sdrfFile), file(metaFile), stdout into COMPILED_DERIVED_CONFIG_WITH_STATUS
        

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
        set val(expName), val(species), file("in/*"), file("in/*"), file("in/*"), val(expStatus) from NEW_OR_CHANGED_EXPERIMENTS

    output:
        set val(expName), val(species), file("*.conf"), file("*.sdrf.txt"), file("*.metadata.tsv") into NEW_OR_RESET_EXPERIMENTS
        
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
            rm -rf $SCXA_RESULTS/$expName/*/\$stage
        done


        # Now we've removed the results, link the results for publishing

        cp -P in/*.sdrf.txt in/*conf in/*.metadata.tsv .
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

