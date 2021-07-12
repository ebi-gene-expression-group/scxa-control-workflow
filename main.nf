#!/usr/bin/env nextflow

tertiaryWorkflow = params.tertiaryWorkflow
overwrite = params.overwrite

dropletProtocols = [ '10xv1', '10xv1a', '10xv1i', '10xv2', '10xv3', 'drop-seq', 'seq-well', '10x5prime' ]

enaSshUser = 'null'

if ( params.containsKey('enaSshUser') ){
    enaSshUser = params.enaSshUser
}

galaxyCredentials = ''
if ( params.containsKey('galaxyCredentials')){
    galaxyCredentials = params.galaxyCredentials
}
galaxyInstance = ''
if ( params.containsKey('galaxyInstance')){
    galaxyInstance = params.galaxyInstance
}

skipQuantification = 'no'
skipAggregation = 'no'
skipTertiary = 'no'

if ( params.containsKey('skipTertiary') && params.skipTertiary == 'yes'){
    
    // Skipping tertiary means that for consistency we must also skip the
    // quantification and aggregations
    
    skipTertiary = 'yes'
    skipQuantification = 'yes'
    skipAggregation = 'yes'
}
else if ( params.containsKey('skipAggregation') && params.skipAggregation == 'yes'){ 

    // Skipping aggregation implies skipping quantifiction
    
    skipQuantification = 'yes'
    skipAggregation = 'yes'

} else if ( params.containsKey('skipQuantification') && params.skipQuantification == 'yes'){
    
    skipQuantification = 'yes'
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
        set val(expName), file(idfFile), file("*.sdrf.txt") optional true into SDRF_IDF
        set val(expName), file("${expName}.cells.txt") optional true into CELLS 

    """
    compileExpMetadata.sh $expName $idfFile $overwrite "${doneSuffix}" 
    """
}

SDRF_IDF
    .join(CELLS, remainder: true)
    .into{
        COMPILED_METADATA_FOR_STATUS
        COMPILED_METADATA_FOR_CELLTYPE
    }

// Check the update status of the experiment

process checkExperimentStatus {

    cache false

    input:
        set val(expName), file(idfFile), file(sdrfFile), file(cellsFile) from COMPILED_METADATA_FOR_STATUS

    output:
        set val(expName), file(idfFile), file(sdrfFile), file(cellsFile), stdout into COMPILED_METADATA_WITH_STATUS
        
    """
    checkExperimentStatus.sh $expName $idfFile $sdrfFile $cellsFile $overwrite
    """
}

// Pipe experiments to different channels depending on update status 

UPDATED_EXPERIMENTS = Channel.create()
NOT_UPDATED_EXPERIMENTS = Channel.create()

COMPILED_METADATA_WITH_STATUS.choice( UPDATED_EXPERIMENTS, NOT_UPDATED_EXPERIMENTS ) {a -> 
    a[4] == 'old' ? 1 : 0
}

// Split the experiment by species- multiple-sepcies experiments will produce
// multiple bundles

process split_by_species {

    input:
        set val(expName), file(idfFile), file(sdrfFile), file(cellsFile), val(expStatus) from UPDATED_EXPERIMENTS.take(params.numExpsProcessedAtOnce)

    output:
        set val(expName), file(idfFile), file('split_by_species/*'), file(cellsFile) into MULTISPECIES_META
        
    """
    splitSDRFBySpecies.sh ${sdrfFile} split_by_species
    """
}

// Extract species from the split SDRF file names

MULTISPECIES_META
    .transpose()
    .map{ row -> tuple( row[0], row[2].toString().split('\\.')[1], row[1], row[2], row[3] ) }
    .into{
        META_WITH_SPECIES
        META_FOR_REF
    }    

// Need to adjust the IDF to match the SDRF. This is necessary for when we're
// doing the SDRF condensation later

process adjust_idf_for_sdrf{
    
    input:
        set val(expName), val(species), file(idfFile), file(sdrfFile), file(cellsFile) from META_WITH_SPECIES

    output:
        set val("${expName}-${species}"), file("${expName}.${species}.idf.txt"), file(sdrfFile), file(cellsFile) into META_WITH_SPECIES_IDF
        set val("${expName}-${species}"), val(expName), val(species) into ES_TAGS

    """
    outFile="${expName}.${species}.idf.txt"
    cp ${idfFile} \${outFile}.tmp
    sed -i 's/${expName}.sdrf.txt/${expName}.${species}.sdrf.txt/' \${outFile}.tmp  
    mv \${outFile}.tmp \${outFile}
    """
}

// Experiment/ species tags

ES_TAGS.unique().into{
    ES_TAGS_FOR_CONFIG
    ES_TAGS_FOR_CONDENSE
    ES_TAGS_FOR_META_MATCHING
    ES_TAGS_FOR_TERTIARY
    ES_TAGS_FOR_REUSE_TERTIARY
    ES_TAGS_FOR_REUSE_AGG
    ES_TAGS_FOR_BUNDLING
}

META_WITH_SPECIES_IDF
    .into{
        META_WITH_SPECIES_FOR_QUANT
        META_WITH_SPECIES_FOR_TERTIARY
    }        

// Generate the configuration file and parse it to derive protocol. Also derive
// an SDRF file containing only things relevant for analysis. This can be used
// to determine if any SDRF changes impact on analysis and necessitate re-run.
// The fundamental processing unit is experiment/ species/ protocol, and this
// is the first process that those combinations are determined. 
//
// NOTE: the references are included in the input just so that they get passed
// through to the protocol-wise channels, not used in the process itself.

process generate_configs {

    cache 'deep'
    
    errorStrategy 'ignore'

    conda 'r-optparse r-data.table r-workflowscriptscommon pyyaml'

    input:
        set val(esTag), val(expName), val(species), file(idfFile), file(sdrfFile), file(cellsFile) from  ES_TAGS_FOR_CONFIG.join(META_WITH_SPECIES_FOR_QUANT)

    output:
        set val(expName), val(species), file("*.${species}.${expName}.conf"), file("*.${species}.${expName}.meta_for_quant.txt"), file("*.${species}.${expName}.meta_for_tertiary.txt") into CONF_ANALYSIS_META_BY_EXP_SPECIES

    """
    # Cells file will be empty (from reminder: true above) for some experiments
    cells_options=
    cells_filesize=\$(stat --printf="%s" \$(readlink ${cellsFile}))
    if [ \$cells_filesize -gt 0 ]; then
        cells_options="--cells=$cellsFile"
    fi
    
    if [ -n "$params.cellAnalysisFields" ]; then
        cells_options="\$cells_options --cell_meta_fields=\\"$params.cellAnalysisFields\\""
        if [ -n "$params.cellTypeField" ]; then
          cells_options="\$cells_options --cell_type_fields=\\"$params.cellTypeField\\""  
        fi
    fi
    
    mkdir -p tmp
    eval "sdrfToNfConf.R \
        --sdrf=\$(readlink $sdrfFile) \
        --idf=$idfFile \$cells_options \
        --name=$expName \
        --verbose \
        --out_conf \$(pwd)"
    """
}
// The transpose step below gives us the protocol-wise configs, converting:
// [ 'E-MTAB-1234', 'rabbit', [ '10xv2.rabbit.E-MTAB-1234.conf', '10xv3.rabbit.E-MTAB-1234.conf' ], [ '10xv2.rabbit.E-MTAB-1234.meta_for_quant.txt', '10xv3.rabbit.E-MTAB-1234.meta_for_quant.txt' ], ['10xv2.rabbit.E-MTAB-1234.meta_for_tertiary.txt', '10xv3.rabbit.E-MTAB-1234.meta_for_tertiary.txt' ] ]
// to:
//
// [ 'E-MTAB-1234', 'rabbit', '10xv2', '10xv2.rabbit.E-MTAB-1234.conf', '10xv2.rabbit.E-MTAB-1234.meta_for_quant.txt', '10xv2.rabbit.E-MTAB-1234.meta_for_tertiary.txt' ]
// [ 'E-MTAB-1234', 'rabbit', '10xv3', '10xv3.rabbit.E-MTAB-1234.conf', '10xv3.rabbit.E-MTAB-1234.meta_for_quant.txt', '10xv3.rabbit.E-MTAB-1234.meta_for_tertiary.txt' ]
//
// Protocol is first part of '.' - separated file name from sdrfToNfConf, so we
// can use simpleName to extract it

CONF_ANALYSIS_META_BY_EXP_SPECIES
    .transpose()
    .map { r -> tuple(r[0], r[1], r[2].simpleName, r[2], r[3], r[4]) }
    .set{ CONF_ANALYSIS_META_BY_EXP_SPECIES_PROTOCOL }


///////////////////////////////////////////////////////////////////////////////
// In this section we use the above-generated config to learn about the
// experiments. What combinations of experiment, species and protocol do we
// have? Which ones are droplet? Have controlled access expeirments been
// properly provided for?
///////////////////////////////////////////////////////////////////////////////

// Check for controlled access status. This will exclude an experiment unless
// all requirements on user, directory etc are in place for a controlled access
// quantification. This process is intentionally not parallelised with others,
// to act as a bottle neck and cause the experiment to cease processing unless
// this check passes.
//
// We also create compact tags to be used here so we don't always have to pass
// around exp/species/protocol

process check_controlled_access{
    
    conda 'pyyaml' 
    
    cache false
    
    errorStrategy 'ignore'

    input:
        set val(expName), val(species), val(protocol), file(confFile), file(metaForQuant), file(metaForTertiary) from CONF_ANALYSIS_META_BY_EXP_SPECIES_PROTOCOL

    output:
        set val("${expName}-${species}-${protocol}"), val(expName), val(species), val(protocol) into ESP_TAGS
        set val("${expName}-${species}-${protocol}"), file(confFile) into CONF_BY_EXP_SPECIES_PROTOCOL
        set val("${expName}-${species}-${protocol}"), file(metaForQuant), file(metaForTertiary) into ANALYSIS_META_BY_EXP_SPECIES_PROTOCOL


    """
    checkForControlledAccess.sh $expName $confFile
    """
}

ANALYSIS_META_BY_EXP_SPECIES_PROTOCOL.into{
    ANALYSIS_META_FOR_CHANGEDCHECK
    ANALYSIS_META_FOR_QUANT
}

CONF_BY_EXP_SPECIES_PROTOCOL.into{
    CONF_BY_EXP_SPECIES_PROTOCOL_FOR_ES_CONFIG
    CONF_BY_EXP_SPECIES_PROTOCOL_FOR_EXTEND
}

// Experiment/ species/ protocol tags

ESP_TAGS.into{
    ESP_TAGS_FOR_AGG_ROUTING
    ESP_TAGS_FOR_ES_CONFIG
    ESP_TAGS_FOR_PROT_LIST
    ESP_TAGS_FOR_CONTAMINATION
    ESP_TAGS_FOR_CHANGEDCHECK
    ESP_TAGS_FOR_EXTEND
    ESP_TAGS_FOR_MARKING
    ESP_TAGS_FOR_REFERENCE
    ESP_TAGS_FOR_SYNCH_REFS
    ESP_TAGS_FOR_QUANT_CHECK
    ESP_TAGS_FOR_REUSE_QUANT
    ESP_TAGS_FOR_QUANT
    ESP_TAGS_FOR_AGGR
    ESP_TAGS_FOR_IS_DROPLET
    ESP_TAGS_FOR_PRIVACY_CHECK
}

// Make comma-separated list of the protocols for each exp/species combo

ESP_TAGS_FOR_PROT_LIST
    .map{r -> tuple(r[1], r[2], r[3])}
    .groupTuple( by: [0,1] )
    .map{r -> tuple((r[0] + '-' +  r[1]).toString(), r[2].join(","))}
    .set{
        PROTOCOLS_BY_EXP_SPECIES
    }

// Take the first config file in a species/ experiment combo to use in
// processes not dependent on protocol

ESP_TAGS_FOR_ES_CONFIG
    .join(CONF_BY_EXP_SPECIES_PROTOCOL_FOR_ES_CONFIG)
    .groupTuple( by: [1,2] )
    .map{ r -> tuple( (r[1] + '-' + r[2]).toString(), r[4][0] ) }
    .into{
       CONF_BY_EXP_SPECIES_FOR_CELLMETA
       CONF_BY_EXP_SPECIES_FOR_CELLTYPE
       CONF_BY_EXP_SPECIES_FOR_TERTIARY 
       CONF_BY_EXP_SPECIES_FOR_BUNDLE 
    }

// Extend config for protocol-wise info

process extend_config_for_protocol{

    input:
        set val(espTag), val(expName), val(species), val(protocol), file(confFile) from ESP_TAGS_FOR_EXTEND.join(CONF_BY_EXP_SPECIES_PROTOCOL_FOR_EXTEND)

    output:
        set val(espTag), file("${expName}.${species}.${protocol}.conf") into EXTENDED_CONF
    
    """
    # Add in protocol-specific parameters by 'include'ing external configs

    confFileOut=${expName}.${species}.${protocol}.conf

    protocolConfig=${baseDir}/conf/protocol/${protocol}.conf
    if [ ! -e \$protocolConfig ]; then
        echo "\$protocolConfig does not exist" 1>&2
        exit 1
    fi 
    
    echo "includeConfig '${baseDir}/params.config'" > \${confFileOut}.tmp
    echo "includeConfig '\$protocolConfig'" >> \${confFileOut}.tmp
    cat $confFile >> \${confFileOut}.tmp 
    mv \${confFileOut}.tmp \${confFileOut}
    """
}

EXTENDED_CONF.into{
    EXTENDED_CONF_FOR_CHANGED_CHECK
    EXTENDED_CONF_FOR_REF
    EXTENDED_CONF_FOR_QUANT
    EXTENDED_CONF_FOR_CELL_TO_LIB
}

// Mark droplet protocols

process mark_droplet {
   
    input:
        set val(espTag), val(expName), val(species), val(protocol) from ESP_TAGS_FOR_MARKING

    output:
        set val(espTag), stdout into IS_DROPLET
    
    script: 

        def isDroplet='False'
        if ( dropletProtocols.contains(protocol) ){
            isDroplet='True'
        }
        
        """
        echo -n "$isDroplet"        
        """
}

// Define droplet status at the level of exp/ species/ protocol AND exp/
// species. Not sure how to treat experiments with droplet AND non-droplet, so
// prevent it 

IS_DROPLET.into{
    IS_DROPLET_PROT
    IS_DROPLET_FOR_EXP_SPECIES_TEST
    IS_DROPLET_FOR_REUSE_QUANT
}

ESP_TAGS_FOR_IS_DROPLET
    .join(IS_DROPLET_FOR_EXP_SPECIES_TEST)
    .map{ r -> tuple ( r[1], r[2], r[4]) }
    .unique()
    .groupTuple( by: [0,1] )
    .map{ r -> tuple(r[0], r[1], r[2].unique()) }
    .set{
        EXP_SPECIES_IS_DROPLETS
    }
    
process check_droplet_status {

    input:
        set val(expName), val(species), val(isDroplets) from EXP_SPECIES_IS_DROPLETS

    output:
        set val("${expName}-${species}"), stdout into IS_DROPLET_EXP_SPECIES

    script: 

        length=isDroplets.size()
        first=isDroplets[0]
    """
    if [ $length -gt 1 ]; then
        echo "Experiment $expName has both droplet and non-droplet protocols" 1>&2
        exit 1
    else
        echo -n "$first"
    fi
    """
}

IS_DROPLET_EXP_SPECIES.into{
    IS_DROPLET_EXP_SPECIES_FOR_CELLMETA
    IS_DROPLET_EXP_SPECIES_FOR_TERTIARY
}

///////////////////////////////////////////////////////////////////////////////
// Processes and transformations in this channel will decide how data are
// re-used or re-analysed
///////////////////////////////////////////////////////////////////////////////

// Where analysis is already present, check that the config has actually
// changed in a way that impacts on analysis. This is distinct from
// checkExperimentStatus(), which just flags where the source SDRF has a newer
// timestamp than the analysis. We could have multiple protocol-wise files, but
// the status is returned at the experiment/ species level (note the use of
// groupTuple() below). Quantification is the only stage peformed at the
// protocol-wise level, and if one protocol of a set is being re-quantified we
// should probbly do the others for consistency in reference usage (we normally
// use the most up-to-date reference when quantifying). We could do reference
// checks to prevent unnecessary re-computation of protocols in multi-protocol
// experiments, but that logic seems unnecessary complication to the logic
// right now.

ESP_TAGS_FOR_CHANGEDCHECK
    .join(EXTENDED_CONF_FOR_CHANGED_CHECK)
    .join(ANALYSIS_META_FOR_CHANGEDCHECK)
    .groupTuple( by: [1,2])
    .set{
        CHANGED_CHECK_INPUTS
    }

process check_experiment_changed{

    publishDir "$SCXA_CONF/study", mode: 'copy', overwrite: true
    
    input:
        set val(espTags), val(expName), val(species), val(protocols), file('confFile/*'), file('metaForQuant/*'), file('metaForTertiary/*') from CHANGED_CHECK_INPUTS

    output:
        set val(expName), val(species), val(espTags), stdout into CHANGE_STATUS
        set file("*.conf"), file("*.meta_for_quant.txt"), file("*.meta_for_tertiary.txt") optional true
        
    """
    expStatus=\$(checkExperimentChanges.sh $expName $species confFile metaForQuant metaForTertiary $skipQuantification $skipAggregation $skipTertiary $overwrite)

    if [ "\$expStatus" != 'unchanged' ]; then
        cp -P confFile/* metaForTertiary/* metaForQuant/* .
    fi
    echo -n "\$expStatus"
    """
}

// For no changes for any of the protocol-wise secions of an experiment, we're
// re-reporting the bundle lines, so no further action required

NEW_OR_CHANGED_EXPERIMENTS = Channel.create()
NOT_CHANGED_EXPERIMENTS = Channel.create()

CHANGE_STATUS.choice( NEW_OR_CHANGED_EXPERIMENTS, NOT_CHANGED_EXPERIMENTS ){a ->
    a[3] == 'unchanged' ? 1 : 0
}

NOT_CHANGED_EXPERIMENTS
    .map{r -> tuple(r[0], r[1])}
    .set{
        NOT_CHANGED_EXPERIMENTS_FOR_BUNDLES
    }

// For experiments we know are due for re-analysis, remove pre-existing
// results (where appropriate). This is also the time to publish the config
// files, overwriting any already in place.

process reset_experiment{
    
    cache 'deep'
   
    input:
        set val(expName), val(species), val(espTags), val(expStatus) from NEW_OR_CHANGED_EXPERIMENTS

    output:
        set val(espTags), val(expStatus) into NEW_OR_RESET_EXPERIMENTS
        
    """
        # Only remove downstream results where we're not re-using them

        if [ "$expStatus" = 'changed_for_quantification' ];then
            reset_stages="quantification reference aggregation scanpy bundle"
        elif [ "$expStatus" = 'changed_for_aggregation' ];then
            reset_stages="aggregation scanpy bundle"
        elif [ "$expStatus" = 'changed_for_tertiary' ];then
            reset_stages="scanpy bundle"
        else
            reset_stages="bundle"
        fi

        rm -f $TMPDIR/${expName}.${species}.galaxystate

        for stage in \$reset_stages; do
            rm -rf $SCXA_RESULTS/$expName/$species/\$stage
        done
    """
}

NEW_OR_RESET_EXPERIMENTS
    .transpose()
    .into{
        EXPERIMENTS_FOR_REF
        EXPERIMENTS_FOR_ANALYSIS
    }

// We've done the right resets, good to go with any new analysis. Route
// experiments different ways depending on in what way a change has been
// observed

TO_RETERTIARY = Channel.create()
TO_REAGGREGATE = Channel.create()
TO_QUANTIFY = Channel.create()

NOT_CHANGED_FOR_QUANT = Channel.create()
NOT_CHANGED_FOR_AGGREGATION = Channel.create()
NOT_CHANGED_FOR_TERTIARY = Channel.create()

// Re-quantify a whole experiment if any of its protocol-wise components
// require it

EXPERIMENTS_FOR_ANALYSIS.join(ESP_TAGS_FOR_QUANT_CHECK).groupTuple( by: [2,3]).choice( TO_QUANTIFY, NOT_CHANGED_FOR_QUANT ) {a -> 
    a[1].contains('changed_for_quantification') ? 0 : 1
}

TO_QUANTIFY
    .transpose()
    .map{ r -> tuple( r[0] ) }
    .into{
        TO_QUANTIFY_FOR_QUANT
        TO_QUANTIFY_FOR_REFERENCE
        TO_QUANTIFY_FOR_PRIVACY_CHECK
    }

NOT_CHANGED_FOR_QUANT
    .transpose()
    .map{ r -> tuple( r[0], r[1]) }
    .into{
        NOT_CHANGED_FOR_QUANT_FOR_REFERENCES
        NOT_CHANGED_FOR_QUANT_FOR_QUANT
        NOT_CHANGED_FOR_QUANT_FOR_REUSE_QUANT
        NOT_CHANGED_FOR_QUANT_FOR_ROUTING
    }

// Quantification-subseqent steps have protocols merged

NOT_CHANGED_FOR_QUANT_FOR_ROUTING.choice( TO_REAGGREGATE, NOT_CHANGED_FOR_AGGREGATION ) {a -> 
    a[1] == 'changed_for_aggregation' ? 0 : 1
}

NOT_CHANGED_FOR_AGGREGATION                                 // esp_tag, expStatus
    .join( ESP_TAGS_FOR_AGG_ROUTING )                       // esp_tag, expStatus, expName, species, protocol
    .map{ r -> tuple( (r[2] + '-' + r[3]).toString(), r[1] ) }           // es_tag, expStatus
    .into{
        NOT_CHANGED_FOR_AGGREGATION_FOR_ROUTING
        NOT_CHANGED_FOR_AGGREGATION_FOR_REUSE_AGG
    }

NOT_CHANGED_FOR_AGGREGATION_FOR_ROUTING.choice( TO_RETERTIARY, NOT_CHANGED_FOR_TERTIARY ) {a -> 
    a[1] == 'changed_for_tertiary' ? 0 : 1
}

NOT_CHANGED_FOR_TERTIARY
    .map{r -> tuple(r[0])}
    .into{
        NOT_CHANGED_FOR_TERTIARY_FOR_REUSE_TERTIARY
        TO_REBUNDLE    
    }

///////////////////////////////////////////////////////////////////////////////
// Processes to derive results to re-use. These will be combined with novel
// results from other channels as appropriate
///////////////////////////////////////////////////////////////////////////////

// This process re-uses previously generated quantifcation results

process reuse_quantifications {
    
    executor 'local'

    input:
        set val(espTag), val(expStatus), val(expName), val(species), val(protocol), val(isDroplet) from NOT_CHANGED_FOR_QUANT_FOR_QUANT.join(ESP_TAGS_FOR_REUSE_QUANT).join(IS_DROPLET_FOR_REUSE_QUANT)

    output:
        set val(espTag), file("results/*") into REUSED_QUANT_RESULTS
        set val("${expName}-${species}"), file("*.fa.gz"), file("*.gtf.gz") into REUSED_REFERENCES
        set val(espTag), file('transcript_to_gene.txt') into REUSED_TRANSCRIPT_TO_GENE

    """
    retrieveStoredFiles.sh $expName $species reference "*.fa.gz *.gtf.gz transcript_to_gene.txt" 
    mkdir -p results

    if [ $isDroplet == 'True' ]; then
        retrieveStoredFiles.sh $expName $species quantification "$protocol/alevin" results/alevin 
    else
        retrieveStoredFiles.sh $expName $species quantification "$protocol/kallisto" results/kallisto
    fi
    """
}

// Use existing aggregations if specified. Note that the aggregations can only
// be reused if the quantifications are reused to maintain consistency, so the
// QUANT_RESULTS channel cannot be allowed to go to reuse_aggregation

process reuse_aggregation {

    executor 'local'
    
    input:
        set val(esTag), val(expStatus), val(expName), val(species) from NOT_CHANGED_FOR_AGGREGATION_FOR_REUSE_AGG.join(ES_TAGS_FOR_REUSE_AGG)

    output:
        set val(esTag), file("matrices/counts_mtx.zip") into REUSED_COUNT_MATRICES
        set val(esTag), file("matrices/tpm_mtx.zip") optional true into REUSED_TPM_MATRICES
        set val(esTag), file("matrices/stats.tsv") optional true into REUSED_KALLISTO_STATS
        set val(esTag), file("${expName}.${species}.condensed-sdrf.tsv") into REUSED_CONDENSED  
        set val(esTag), file("${expName}-${species}.metadata.matched.tsv") into REUSED_MATCHED_META

    """
        retrieveStoredFiles.sh $expName $species aggregation matrices
        retrieveStoredFiles.sh $expName $species metadata "${expName}-${species}.metadata.matched.tsv ${expName}.${species}.condensed-sdrf.tsv"
    """
}

REUSED_COUNT_MATRICES.into{
    REUSED_COUNT_MATRICES_FOR_TERTIARY
    REUSED_COUNT_MATRICES_FOR_BUNDLING
}

// Derive existing tertiary results for cases where we're just re-bundling

process reuse_tertiary {

    executor 'local'
    
    input:
        set val(esTag), val(expName), val(species) from NOT_CHANGED_FOR_TERTIARY_FOR_REUSE_TERTIARY.join(ES_TAGS_FOR_REUSE_TERTIARY)
    
    output:
        set val(esTag), file("matrices/raw_filtered.zip"), file("matrices/filtered_normalised.zip"), file("clusters_for_bundle.txt"), file("umap"), file("tsne"), file("markers"), file('clustering_software_versions.txt'), file('project.h5ad') into REUSED_TERTIARY_RESULTS

    """
        retrieveStoredFiles.sh $expName $species scanpy "matrices clusters_for_bundle.txt umap tsne markers clustering_software_versions.txt project.h5ad"
    """
}

////////////////////////////////////////////////////////////////
// Processes to derive new results
////////////////////////////////////////////////////////////////

// Is experiment public or private?

process check_privacy {

    executor 'local'

    input:
        set val(espTag), val(expName), val(species), val(protocol) from TO_QUANTIFY_FOR_PRIVACY_CHECK.join(ESP_TAGS_FOR_PRIVACY_CHECK)

    output:
        set val(espTag), stdout into PRIVACY_STATUS

    """
    privacyStatus=\$(wget -O - http://peach.ebi.ac.uk:8480/api/privacy.txt?acc=${expName} 2>/dev/null | tr "\\t" "\\n" | awk -F':' '/privacy/ {print \$2}')
    if [ -z "\$privacyStatus" ]; then
      privacyStatus=public
    fi    
    echo -n "\$privacyStatus"
    """
}

// Symlink to correct reference file, re-using the reference from last
// quantification where appropriate. spikes are used in quantification only.

process add_reference {

    conda "${baseDir}/envs/refgenie.yml"

    cache 'deep'
    
    errorStrategy { task.attempt<=3 ? 'retry' : 'finish' }
        
    publishDir "$SCXA_RESULTS/$expName/$species/reference", mode: 'copyNoFollow', overwrite: true, pattern: '{*.fa.gz,*.gtf.gz}'

    input:
        set val(espTag), file(confFile), val(expName), val(species), val(protocol) from EXTENDED_CONF_FOR_REF.join(TO_QUANTIFY_FOR_REFERENCE).join(ESP_TAGS_FOR_REFERENCE)

    output:
        set val(espTag), file('spiked/*.fa.gz'), file("spiked/*.gtf.gz"), file('spiked/*.idx'), file('spiked/salmon_index') into PREPARED_REFERENCES
        set val(espTag), val(expName), val(species), file("spiked/*.gtf.gz") into GTF_FOR_T2GENE
        set val("${expName}-${species}"), file("*.fa.gz"), file("*.gtf.gz") into NEW_REFERENCES_FOR_DOWNSTREAM

    """
    refgenieSeek.sh $species ${params.islReferenceType} "" \$(pwd)

    # If we have spikes, get the spikes references for use in quantification

    spikes=\$(parseNfConfig.py --paramFile $confFile --paramKeys params,spikes)
    if [ \$spikes != 'None' ]; then
        refgenieSeek.sh $species ${params.islReferenceType} "\$spikes" spiked
    else
        # Copy unspiked symlinks to spiked if no spikes are required
        mkdir spiked 
        cp -P *.gtf.gz *.fa.gz salmon_index *.idx spiked
    fi
    """
}

// This process is parameterised by expName, species, protocol, to allow us to
// control contamination index use for those species in future.

process find_contamination_index {

    conda "${baseDir}/envs/refgenie.yml"
  
    input:
        set val(espTag), val(expName), val(species), val(protocol) from ESP_TAGS_FOR_CONTAMINATION 

    output:
        set val(espTag), stdout into CONTAMINATION_INDEX
     
    """
    refgenie seek contamination/bowtie2_index | tr -d \'\\n\'
    """
}

// References used in tertiary analysis and bundling are a combination of the
// reused and new references depending on individual experiment statuses

NEW_REFERENCES_FOR_DOWNSTREAM.unique()
    .concat(REUSED_REFERENCES.unique())
    .into{
        REFERENCES_FOR_TERTIARY
        REFERENCES_FOR_BUNDLING
    }

// Synchronise the GTF and the FASTA 

process transcript_to_gene {

    publishDir "$SCXA_RESULTS/$expName/$species/reference", mode: 'copy', overwrite: true, pattern: 'transcript_to_gene.txt'
    
    conda "${baseDir}/envs/atlas-gene-annotation-manipulation.yml"
    
    cache 'deep'

    memory { 5.GB * task.attempt }

    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3  ? 'retry' : 'ignore' }
    maxRetries 3
        
    input:
        set val(espTag), val(expName), val(species), file(referenceGtf) from GTF_FOR_T2GENE

    output:
        set val(espTag), file('transcript_to_gene.txt') into TRANSCRIPT_TO_GENE

    """
    gtf2featureAnnotation.R --gtf-file ${referenceGtf} --no-header \
        --version-transcripts --feature-type "transcript" --first-field \
        "transcript_id" --output-file transcript_to_gene.txt --fields \
        "transcript_id,gene_id"    
    """
}

TRANSCRIPT_TO_GENE
    .concat(REUSED_TRANSCRIPT_TO_GENE)
    .into{
        TRANSCRIPT_TO_GENE_QUANT
        TRANSCRIPT_TO_GENE_AGGR
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

// Compile the info we need for fresh quantification

TO_QUANTIFY_FOR_QUANT
    .join(IS_DROPLET_PROT)
    .join(ESP_TAGS_FOR_QUANT)
    .join(EXTENDED_CONF_FOR_QUANT)
    .join(ANALYSIS_META_FOR_QUANT)
    .join(PREPARED_REFERENCES.map{tuple(it[0], it[1], it[2], it[3], it[4])})
    .join(CONTAMINATION_INDEX)
    .join(TRANSCRIPT_TO_GENE_QUANT)
    .join(PRIVACY_STATUS)
    .set{
        QUANT_INPUTS
    }

// Separate droplet and smart-type experiments

DROPLET_INPUTS = Channel.create()
SMART_INPUTS = Channel.create()

QUANT_INPUTS.choice( SMART_INPUTS, DROPLET_INPUTS ) {a ->
    a[1] == 'True' ? 1 : 0 
}

// Run quantification with https://github.com/ebi-gene-expression-group/scxa-smartseq-quantification-workflow

process smart_quantify {

    maxForks params.maxConcurrentQuantifications

    conda "${baseDir}/envs/nextflow.yml"
    
    cache 'deep'

    publishDir "$SCXA_RESULTS/$expName/$species/quantification/$protocol", mode: 'copy', overwrite: true
    
    memory { 10.GB * task.attempt }
    errorStrategy { task.attempt<=10 ? 'retry' : 'ignore' }
    maxRetries 10
    
    input:
        set val(espTag), val(isDroplet), val(expName), val(species), val(protocol), file(confFile), file(metaForQuant), file(metaForTertiary), file(referenceFasta), file(referenceGtf), file(kallistoIndex), file(salmonIndex), val(contaminationIndex), file(transcriptToGene), val(privacyStatus) from SMART_INPUTS
        val flag from INIT_DONE_SMART

    output:
        set val(espTag), file ("kallisto") into SMART_KALLISTO_DIRS 
        set val(espTag), file ("qc") into SMART_QUANT_QC
        file('quantification.log')    

    """
    submitQuantificationWorkflow.sh 'smart-seq' "$expName" "$species" "$protocol" "$confFile" "$metaForQuant" "$referenceFasta" "$transcriptToGene" "$contaminationIndex" "$kallistoIndex" "$salmonIndex" "$enaSshUser" "$privacyStatus"
    """
}

// Run quantification with https://github.com/ebi-gene-expression-group/scxa-droplet-quantification-workflow

process droplet_quantify {

    maxForks params.maxConcurrentQuantifications

    conda "${baseDir}/envs/nextflow.yml"

    cache 'deep'
    
    publishDir "$SCXA_RESULTS/$expName/$species/quantification/$protocol", mode: 'copy', overwrite: true
    
    memory { 10.GB * task.attempt }
    errorStrategy { task.attempt<=5 ? 'retry' : 'ignore' }

    input:
        set val(espTag), val(isDroplet), val(expName), val(species), val(protocol), file(confFile), file(metaForQuant), file(metaForTertiary), file(referenceFasta), file(referenceGtf), val(contaminationIndex), file(kallistoIndex), file(salmonIndex), file(transcriptToGene), val(privacyStatus) from DROPLET_INPUTS
        val flag from INIT_DONE_DROPLET

    output:
        set val(espTag), file("alevin") into ALEVIN_DROPLET_DIRS
        file('quantification.log')    

    """
    submitQuantificationWorkflow.sh 'droplet' "$expName" "$species" "$protocol" "$confFile" "$metaForQuant" "$referenceFasta" "$transcriptToGene" "$contaminationIndex" "$kallistoIndex" "$salmonIndex" "$enaSshUser" "$privacyStatus"
    """
}

// Collect the smart and droplet workflow results

SMART_KALLISTO_DIRS
    .concat(ALEVIN_DROPLET_DIRS)
    .set { 
        NEW_QUANT_RESULTS 
    }

// Aggregation just needs the quantification results and the transcript/ gene
// mappings. Input is the combination of new quantifications, and reused
// existing quantifications specified for reaggregation.

ESP_TAGS_FOR_AGGR
    .join(NEW_QUANT_RESULTS.concat(TO_REAGGREGATE.map{r -> tuple(r[0])}.join(REUSED_QUANT_RESULTS)))
    .join(TRANSCRIPT_TO_GENE_AGGR)
    .groupTuple( by: [1,2] )
    .map{ r -> tuple( r[1], r[2], r[3], r[4], r[5] ) }
    .set{ GROUPED_QUANTIFICATION_RESULTS }

process aggregate {
    
    conda "${baseDir}/envs/nextflow.yml"
    
    maxForks 6

    publishDir "$SCXA_RESULTS/$expName/$species/aggregation", mode: 'copy', overwrite: true
    
    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3 ? 'retry' : 'finish' }
    maxRetries 20
    
    input:
        set val(expName), val(species), file('quant_results/??/protocol'), file('quant_results/??/*'), file('quant_results/??/*') from GROUPED_QUANTIFICATION_RESULTS   

    output:
        set val("${expName}-${species}"), file("matrices/counts_mtx.zip") into NEW_COUNT_MATRICES
        set val("${expName}-${species}"), file("matrices/tpm_mtx.zip") optional true into NEW_TPM_MATRICES
        set val("${expName}-${species}"), file("matrices/kallisto_stats.tsv") optional true into NEW_KALLISTO_STATS
        set val("${expName}-${species}"), file("matrices/alevin_stats.tsv") optional true into NEW_ALEVIN_STATS

    """
        submitAggregationWorkflow.sh "$expName" "$species"
    """
}

NEW_COUNT_MATRICES
    .into{
        NEW_COUNT_MATRICES_FOR_CELLMETA
        NEW_COUNT_MATRICES_FOR_CELL_TO_LIB
        NEW_COUNT_MATRICES_FOR_META_MATCHING
        NEW_COUNT_MATRICES_FOR_TERTIARY
        NEW_COUNT_MATRICES_FOR_BUNDLING
    }

// Get channel with the first config file for each exp/ species pair. For when
// we need a config but don't need the protocol-specific bits, for example a
// process needs to know the column with technical replicates etc.

IS_DROPLET_EXP_SPECIES_FOR_CELLMETA         // esTag, isDroplet
    .join(CONF_BY_EXP_SPECIES_FOR_CELLMETA) // esTag, isDroplet, confFile
    .join(NEW_COUNT_MATRICES_FOR_CELLMETA)  // esTag, isDroplet, confFile, countMatrix
    .set{
        CELLMETA_INPUTS
    }

// Separate Droplet from non-droplet, to allow for some specialised processing

DROPLET_CELL_TO_LIB_INPUT = Channel.create()
SMART_CELL_TO_LIB_INPUT = Channel.create()

CELLMETA_INPUTS
    .map{r -> tuple(r[0], r[1], r[2], r[3], file('NO_FILE')) }
    .choice( DROPLET_CELL_TO_LIB_INPUT, SMART_CELL_TO_LIB_INPUT ) {a -> 
        a[1] == 'True' ? 0 : 1
    }

// SMART does not need the cell/lib mapping, so pass the empty value through

SMART_CELL_TO_LIB_INPUT
    .map{ r -> tuple(r[0], r[4]) }
    .set{SMART_CELL_TO_LIB}

// Make a cell-run mapping, required by the SDRF condense process to 'explode'
// the SDRF annotations

process cell_run_mapping {
   
    cache 'deep'
    
    input:
        set val(esTag), val(isDroplet), file(confFile), file(countMatrix), file(emptyFile) from DROPLET_CELL_TO_LIB_INPUT
 
    output:
        set val(esTag), file('cell_to_library.txt') into DROPLET_CELL_TO_LIB
 
    """
    makeCellLibraryMapping.sh $countMatrix $confFile cell_to_library.txt 
    """
}

// Now make a condensed SDRF file. For this to operate correctly with droplet
// data, the cell-library mappings file must be present

ES_TAGS_FOR_CONDENSE
    .join(SMART_CELL_TO_LIB.concat(DROPLET_CELL_TO_LIB))
    .join(META_WITH_SPECIES_FOR_TERTIARY)
    .set{
       CONDENSE_INPUTS
   }

process condense_sdrf {
        
    publishDir "$SCXA_RESULTS/$expName/$species/metadata", mode: 'copy', overwrite: true
    
    cache 'deep'
        
    maxForks 2
    
    conda "${baseDir}/envs/atlas-experiment-metadata.yml"
    
    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 || task.exitStatus ==  255 ? 'retry' : 'ignore' }
    maxRetries 10

    input:
        set val(esTag), val(expName), val(species), file(cell_to_lib), file(idfFile), file(origSdrfFile), file(cellsFile) from CONDENSE_INPUTS 

    output:
        set val(esTag), file("${expName}.${species}.condensed-sdrf.tsv") into CONDENSED 

    """
    cellTypeFields=
    if [ -n "$params.cellTypeField" ]; then
        cellTypeFields="-t \\"$params.cellTypeField\\""
    fi
    echo -e "exclusions: $ZOOMA_EXCLUSIONS"
    eval "single_cell_condensed_sdrf.sh -e $expName -f $idfFile -o \$(pwd) -z $ZOOMA_EXCLUSIONS \$cellTypeFields"
    mv ${expName}.condensed-sdrf.tsv "${expName}.${species}.condensed-sdrf.tsv"
    """        
}

CONDENSED.into{
    CONDENSED_FOR_META
    CONDENSED_FOR_BUNDLING
}

// 'unmelt' the condensed SDRF to get a metadata table to pass for tertiary
// analysis

process unmelt_condensed_sdrf {
        
    conda "${baseDir}/envs/atlas-experiment-metadata.yml"
    
    cache 'deep'
    
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 ? 'retry' : 'finish' }
    
    memory { 4.GB * task.attempt }
    
    maxRetries 20
    
    input:
        set val(esTag), file(condensedSdrf) from CONDENSED_FOR_META

    output:
       set val(esTag), file("${esTag}.metadata.tsv") into UNMELTED_META 
        
    """
    unmelt_condensed.R -i $condensedSdrf -o ${esTag}.metadata.tsv --retain-types --has-ontology
    """        
}

// Match the cell metadata to the expression matrix 

process match_metadata_to_cells {
    
    publishDir "$SCXA_RESULTS/$expName/$species/metadata", mode: 'copy', overwrite: true
    
    cache 'deep'
    
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 ? 'retry' : 'finish' }
    
    memory { 4.GB * task.attempt }
    
    maxRetries 20

    input:
        set val(esTag), file(cellMeta), val(expName), val(species), file(countMatrix) from UNMELTED_META.join(ES_TAGS_FOR_META_MATCHING).join(NEW_COUNT_MATRICES_FOR_META_MATCHING)

    output:
       set val(esTag), file("${esTag}.metadata.matched.tsv") into NEW_MATCHED_META 

    """
    matchMetadataToCells.sh $cellMeta $countMatrix ${esTag}.metadata.matched.tsv 
    """
}

NEW_MATCHED_META
    .concat(REUSED_MATCHED_META)
    .into{
        MATCHED_META_FOR_TERTIARY
        MATCHED_META_FOR_BUNDLING
    }

// Run tertiary analysis anew. Input is newly aggregated studies and
// pre-existing aggregations flagged for tertiary

ES_TAGS_FOR_TERTIARY                                                                                            // esTag, expName, species
    .join(CONF_BY_EXP_SPECIES_FOR_TERTIARY)                                                                     // esTag, expName, species, confFile 
    .join(NEW_COUNT_MATRICES_FOR_TERTIARY.concat(TO_RETERTIARY.map{ r -> tuple(r[0])}.join(REUSED_COUNT_MATRICES_FOR_TERTIARY)))       // esTag, expName, species, confFile, countMatrix
    .join(REFERENCES_FOR_TERTIARY)                                                                              // esTag, expName, species, confFile, countMatrix, referenceFasta, referenceGtf
    .join(MATCHED_META_FOR_TERTIARY)                                                                            // esTag, expName, species, confFile, countMatrix, referenceFasta, referenceGtf, cellMetadata
    .join(IS_DROPLET_EXP_SPECIES_FOR_TERTIARY)                                                                  // esTag, expName, species, confFile, countMatrix, referenceFasta, referenceGtf, cellMetadata, isDroplet
    .set{TERTIARY_INPUTS}

process tertiary {
    
    cache 'deep'
    
    // Exit status of 3 is just Galaxy being annoying with history
    // deletion, no cause to error

    validExitStatus 0,3

    maxForks params.maxConcurrentScanpyGalaxy

    conda "${baseDir}/envs/galaxy-workflow-executor.yml"

    publishDir "$SCXA_RESULTS/$expName/$species/scanpy", mode: 'copy', overwrite: true
    
    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 || task.attempt<=3  ? 'retry' : 'finish' }
    maxRetries 3
      
    input:
        set val(esTag), val(expName), val(species), file(confFile), file(countMatrix), file(referenceFasta), file(referenceGtf), file(cellMetadata), val(isDroplet) from TERTIARY_INPUTS

    output:
        set val(esTag), file("matrices/raw_filtered.zip"), file("matrices/filtered_normalised.zip"), file("clusters_for_bundle.txt"), file("umap"), file("tsne"), file("markers"), file('clustering_software_versions.txt'), file('project.h5ad') into NEW_TERTIARY_RESULTS

    script:

        """
            submitTertiaryWorkflow.sh "$expName" "$species" "$confFile" "$countMatrix" "$referenceGtf" "$cellMetadata" "$isDroplet" "$galaxyCredentials" "$galaxyInstance"
        """
}

// Bundling requires pretty much everything, so we need to gather a few
// channels. Input experiments are those with new tertiary analysis, plus those
// with re-used tertiary analysis flagged for re-bundling

NEW_TERTIARY_RESULTS                                                
    .concat(TO_REBUNDLE.join(REUSED_TERTIARY_RESULTS))                                  // esTag, filteredMatrix, normalisedMatrix, clusters, umap, tsne, markers, softwareReport, projectFile
    .join(ES_TAGS_FOR_BUNDLING)                                                         // esTag, filteredMatrix, normalisedMatrix, clusters, umap, tsne, markers, softwareReport, projectFile, expName, species
    .join(PROTOCOLS_BY_EXP_SPECIES)                                                     // esTag, filteredMatrix, normalisedMatrix, clusters, umap, tsne, markers, softwareReport, projectFile, expName, species, protocolList
    .join(CONF_BY_EXP_SPECIES_FOR_BUNDLE)                                               // esTag, filteredMatrix, normalisedMatrix, clusters, umap, tsne, markers, softwareReport, projectFile, expName, species, protocolList, confFile
    .join(NEW_COUNT_MATRICES_FOR_BUNDLING.concat(REUSED_COUNT_MATRICES_FOR_BUNDLING))   // esTag, filteredMatrix, normalisedMatrix, clusters, umap, tsne, markers, softwareReport, projectFile, expName, species, protocolList, confFile, rawMatrix
    .join(NEW_TPM_MATRICES.concat(REUSED_TPM_MATRICES), remainder: true)                                 // esTag, filteredMatrix, normalisedMatrix, clusters, umap, tsne, markers, softwareReport, projectFile, expName, species, protocolList, confFile, rawMatrix, tpmMatrix
    .join(REFERENCES_FOR_BUNDLING)                                                      // esTag, filteredMatrix, normalisedMatrix, clusters, umap, tsne, markers, softwareReport, projectFile, expName, species, protocolList, confFile, rawMatrix, tpmMatrix, referenceFasta, referenceGtf
    .join(CONDENSED_FOR_BUNDLING.concat(REUSED_CONDENSED))                              // esTag, filteredMatrix, normalisedMatrix, clusters, umap, tsne, markers, softwareReport, projectFile, expName, species, protocolList, confFile, rawMatrix, tpmMatrix, referenceFasta, referenceGtf, condensedSdrf
    .join(MATCHED_META_FOR_BUNDLING)                                                    // esTag, filteredMatrix, normalisedMatrix, clusters, umap, tsne, markers, softwareReport, projectFile, expName, species, protocolList, confFile, rawMatrix, tpmMatrix, referenceFasta, referenceGtf, condensedSdrf, cellMetadata 
    .set{BUNDLE_INPUTS}

process bundle {
    
    conda "${baseDir}/envs/nextflow.yml"

    publishDir "$SCXA_RESULTS/$expName/$species", mode: 'copy', overwrite: true
    
    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'ignore' }
    maxRetries 20
    
    input:
        set val(esTag), file(filteredMatrix), file(normalisedMatrix), file(clusters), file('*'), file('*'), file('*'), file(softwareReport), file(projectFile), val(expName), val(species), val(protocolList), file(confFile), file(rawMatrix), file(tpmMatrix), file(referenceFasta), file(referenceGtf), file(condensedSdrf), file(cellMetadata) from BUNDLE_INPUTS
            
    output:
        file('bundle/software.tsv')
        file('bundle/filtered_normalised/genes.tsv.gz')
        file('bundle/filtered_normalised/barcodes.tsv.gz')
        file('bundle/filtered_normalised/matrix.mtx.gz')
        file('bundle/filtered_normalised/cell_to_library.txt') optional true
        file('bundle/filtered_normalised/filtered_normalised.tsv') optional true
        file('bundle/raw_filtered/genes.tsv.gz')
        file('bundle/raw_filtered/barcodes.tsv.gz')
        file('bundle/raw_filtered/matrix.mtx.gz')
        file('bundle/raw_filtered/raw_filtered.tsv') optional true
        file('bundle/raw_filtered/cell_to_library.txt') optional true
        file('bundle/tpm/genes.tsv.gz') optional true
        file('bundle/tpm/barcodes.tsv.gz') optional true
        file('bundle/tpm/matrix.mtx.gz') optional true
        file('bundle/tpm/tpm.tsv') optional true
        file('bundle/raw/genes.tsv.gz')
        file('bundle/raw/barcodes.tsv.gz')
        file('bundle/raw/matrix.mtx.gz')
        file('bundle/raw/raw.tsv') optional true
        file('bundle/raw/cell_to_library.txt') optional true
        file('bundle/tpm_filtered/genes.tsv.gz') optional true
        file('bundle/tpm_filtered/barcodes.tsv.gz') optional true
        file('bundle/tpm_filtered/matrix.mtx.gz') optional true
        file('bundle/tpm_filtered/tpm_filtered.tsv') optional true
        file("bundle/${expName}.cell_metadata.tsv")
        file("bundle/${expName}.condensed-sdrf.tsv")
        file("bundle/reference/${referenceFasta}")
        file("bundle/reference/${referenceGtf}")
        file('bundle/tsne_perplexity_*.tsv')
        file('bundle/umap_n_neighbors_*.tsv')
        file('bundle/markers_*.tsv') optional true
        file('bundle/clusters_for_bundle.txt')
        file('bundle/MANIFEST')
        file('bundle.log')
        file('bundle/filtered_normalised_stats.csv')
        file('bundle/tpm_filtered_stats.csv') optional true
        file('bundleLines.txt') into NEW_BUNDLES
        file("bundle/${expName}.project.h5ad")
        
        """
            submitBundleWorkflow.sh "$expName" "$species" "$protocolList" "$confFile" "$referenceFasta" "$referenceGtf" "$tertiaryWorkflow" "$condensedSdrf" "$cellMetadata" "$rawMatrix" "$filteredMatrix" "$normalisedMatrix" "$tpmMatrix" "$clusters" markers tsne umap $softwareReport $projectFile
        """
}

// Bundle lines for things that haven't changed at all

process get_not_updated_bundles {

    input:
        val expName from NOT_UPDATED_EXPERIMENTS.map{r -> r[0]}

    output:
        file('bundleLines.txt') into NOT_UPDATED_BUNDLES

    """
        ls \$SCXA_RESULTS/$expName/*/bundle/MANIFEST | while read -r l; do
            species=\$(echo \$l | awk -F'/' '{print \$(NF-2)}' | tr -d \'\\n\')
            echo -e "$expName\\t\$species\\t$SCXA_RESULTS/$expName/\$species/bundle"
        done > bundleLines.txt

    """
}

// Generate bundle lines for the things updated but not actually changed for
// analysis

process get_not_changed_bundles {

    input:
        set val(expName), val(species) from NOT_CHANGED_EXPERIMENTS_FOR_BUNDLES

    output:
         file('bundleLines.txt') into NOT_CHANGED_BUNDLES

    """
        # Update the manifest time stamps to prevent re-config next time
        touch -m $SCXA_RESULTS/$expName/$species/bundle/MANIFEST

        echo -e "$expName\\t$species\\t$SCXA_RESULTS/$expName/$species/bundle" > bundleLines.txt
    """
}

// Record the completed bundles

NOT_UPDATED_BUNDLES
    .concat( NOT_CHANGED_BUNDLES )
    .concat ( NEW_BUNDLES )
    .collectFile(name: 'these.all.done.txt', sort: true)
    .set{
        BUNDLE_LINES
    }

// For each experiment and species with a completed bundle, remove the work
// dir. This can take a little while for large experiments, so background the
// process so future runs are not delayed.

process cleanup {
    
    executor 'local'

    errorStrategy 'ignore'
    
    publishDir "$SCXA_RESULTS", mode: 'copy', overwrite: true
    
    input:
        file(bundleLines) from BUNDLE_LINES
    
    output:
        file("all.done${doneSuffix}.txt") into DONEFILE

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

// Make the HTML report 

process report_table {

    cache false
    
    conda "${baseDir}/envs/dt.yml"
    
    publishDir "$SCXA_HTML_DIR", mode: 'move', overwrite: true

    input:
        file(doneFile) from DONEFILE    

    output:
        file("scxa_analysis_status.html")

    """
    makeExperimentStatusTable.R scxa_analysis_status.html 
    """
}
