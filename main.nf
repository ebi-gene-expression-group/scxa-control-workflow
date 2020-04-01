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
    .set{
        COMPILED_METADATA
    }

// Check the update status of the experiment

process checkExperimentStatus {

    input:
        set val(expName), file(idfFile), file(sdrfFile), file(cellsFile) from COMPILED_METADATA

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

// Check the update status of the experiment

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
    .map{ row-> tuple( row[0], row[2].toString().split('\\.')[1], row[1], row[2], row[3] ) }
    .into{
        META_WITH_SPECIES
        META_FOR_REF
    }    

// Locate reference files

process add_reference {

    conda 'pyyaml' 
    
    cache 'deep'
    
    errorStrategy { task.attempt<=3 ? 'retry' : 'finish' }

    input:
        set val(expName), val(species) from META_FOR_REF.map{ r -> tuple(r[0], r[1]) }
    
    output:
        set val(expName), val(species), file("*.fa.gz"), file("*.gtf.gz"), stdout into REFERENCES

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

REFERENCES.into{
    REFERENCES_FOR_QUANT
    REFERENCES_FOR_TERTIARY
    REFERENCES_FOR_BUNDLE
}

// Need to adjust the IDF to match the SDRF. This is necessary for when we're
// doing the SDRF condensation later

process adjust_idf_for_sdrf{
    
    input:
        set val(expName), val(species), file(idfFile), file(sdrfFile), file(cellsFile) from META_WITH_SPECIES

    output:
        set val(expName), val(species), file("${expName}.${species}.idf.txt"), file(sdrfFile), file(cellsFile) into META_WITH_SPECIES_IDF

    """
    outFile="${expName}.${species}.idf.txt"
    cp ${idfFile} \${outFile}.tmp
    sed -i 's/${expName}.sdrf.txt/${expName}.${species}.sdrf.txt/' \${outFile}.tmp  
    mv \${outFile}.tmp \${outFile}
    """
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
        set val(expName), val(species), file(idfFile), file(sdrfFile), file(cellsFile), file(referenceFasta), file(referenceGtf), file(contIndex) from META_WITH_SPECIES_FOR_QUANT.join(REFERENCES_FOR_QUANT, by: [0,1])

    output:
        set val(expName), val(species), file(referenceFasta), file(referenceGtf), file(contIndex), file("*.${species}.${expName}.conf") into CONF_REF_BY_EXP_SPECIES
        set val(expName), val(species), file("*.${species}.${expName}.sdrf.txt") into CONFSDRF_BY_EXP_SPECIES

    """
    mkdir -p tmp
    sdrfToNfConf.R \
        --sdrf=\$(readlink $sdrfFile) \
        --idf=$idfFile \
        --name=$expName \
        --verbose \
        --out_conf \$(pwd)
    """
}

CONF_REF_BY_EXP_SPECIES.into{
    CONF_REF_BY_EXP_SPECIES_FOR_QUANT
    CONF_REF_BY_EXP_SPECIES_FOR_TERTIARY
    CONF_REF_BY_EXP_SPECIES_FOR_BUNDLE
}

// Protocol is first part of '.' - separated file name from sdrfToNfConf, so we
// can use simpleName to extract it

CONF_REF_BY_EXP_SPECIES_FOR_QUANT
    .transpose()
    .map { r -> tuple(r[0], r[1], r[5].simpleName, r[2], r[3], r[4], r[5]) }
    .set{ CONF_REF_BY_EXP_SPECIES_PROTOCOL }

CONFSDRF_BY_EXP_SPECIES
    .transpose()
    .map { r -> tuple(r[0], r[1], r[2].simpleName, r[2]) }
    .set{ CONFSDRF_BY_EXP_SPECIES_PROTOCOL }


process extend_config_for_protocol{

    input:
        set val(expName), val(species), val(protocol),  file(referenceFasta), file(referenceGtf), file(contIndex), file(confFile), file(sdrfFile) from CONF_REF_BY_EXP_SPECIES_PROTOCOL.join(CONFSDRF_BY_EXP_SPECIES_PROTOCOL, by: [0,1,2])

    output:
        set val("${expName}-${species}-${protocol}"), val(expName), val(species), val(protocol) into TAGS 
        set val("${expName}-${species}-${protocol}"), file("${expName}.${species}.${protocol}.conf"), file(sdrfFile) into TAGGED_CONF 
        set val("${expName}-${species}-${protocol}"), file(referenceFasta), file(referenceGtf), file(contIndex) into TAGGED_REFERENCES
    
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


TAGS.into{
    TAGS_FOR_MARKING
    TAGS_FOR_REF
    TAGS_FOR_CHANGEDCHECK
    TAGS_FOR_QUANT
    TAGS_FOR_AGGR
    TAGS_FOR_IS_DROPLET
    TAGS_FOR_TERTIARY
    TAGS_FOR_BUNDLE
}

TAGGED_CONF.into{
    CONF_FOR_CHANGEDCHECK
    CONF_FOR_QUANT
    PRE_CONF_FOR_TERTIARY
}

// Mark droplet protocols

process mark_droplet {
   
    input:
        set val(tag), val(expName), val(species), val(protocol) from TAGS_FOR_MARKING

    output:
        set val(tag), stdout into IS_DROPLET
    
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
}

TAGS_FOR_IS_DROPLET
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
        set val(expName), val(species), stdout into IS_DROPLET_EXP_SPECIES

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

// Where analysis is already present, check that the config has actually changed
// in a way that impacts on analysis. This is distinct from
// checkExperimentStatus(), which just flags where the source SDRF has a newer
// timestamp than the analysis.

process check_experiment_changed{

    input:
        set val(tag), val(expName), val(species), val(protocol), file(confFile), file(sdrfFile) from TAGS_FOR_CHANGEDCHECK.join(CONF_FOR_CHANGEDCHECK)

    output:
        set val(tag), val(expName), val(species), val(protocol), file(confFile), file(sdrfFile), stdout into COMPILED_DERIVED_CONFIG_WITH_STATUS
        

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
        set val(expName), val(species) from NOT_CHANGED_EXPERIMENTS.map{r -> tuple(r[1], r[2])}

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

//[E-MTAB-6077-danio_rerio-smart-seq, E-MTAB-6077, danio_rerio, smart-seq, changed]

process reset_experiment{
    
    publishDir "$SCXA_CONF/study", mode: 'copy', overwrite: true

    input:
        set val(tag), val(expName), val(species), val(protocol), file("in/*"), file("in/*"), val(expStatus) from NEW_OR_CHANGED_EXPERIMENTS

    output:
        set val(tag), file("*.conf"), file("*.sdrf.txt") into NEW_OR_RESET_EXPERIMENTS
        
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

        cp -P in/*.sdrf.txt in/*conf .
    """
}

// Prepare a reference depending on spikes

process prepare_reference {

    conda 'pyyaml' 
    
    cache 'deep'
    
    errorStrategy { task.attempt<=3 ? 'retry' : 'finish' }

    input:
        set val(tag), file(confFile), file(sdrfFile), file(referenceFasta), file(referenceGtf), val(contaminationIndex) from NEW_OR_RESET_EXPERIMENTS.join(TAGGED_REFERENCES)
    
    output:
        set val(tag), file("out/*.fa.gz"), file("out/*.gtf.gz"), val(contaminationIndex) into PREPARED_REFERENCES

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
        set val(tag), file(referenceFasta), file(referenceGtf), val(contaminationIndex) from PREPARED_REFERENCES

    output:
        set val(tag), file('cleanedCdna.fa.gz'), file(referenceGtf), val(contaminationIndex) into CLEANED_REFERENCES
        set val(tag), file('transcript_to_gene.txt') into TRANSCRIPT_TO_GENE

    """
    gtf2featureAnnotation.R --gtf-file ${referenceGtf} --no-header --version-transcripts --filter-cdnas ${referenceFasta} \
        --filter-cdnas-field "transcript_id" --filter-cdnas-output cleanedCdna.fa.gz --feature-type "transcript" \
        --first-field "transcript_id" --output-file transcript_to_gene.txt --fields "transcript_id,gene_id"    
    """
}

TRANSCRIPT_TO_GENE.into{
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

// Compile the info we need for quantification

IS_DROPLET_PROT
    .join(TAGS_FOR_QUANT)
    .join(CONF_FOR_QUANT)
    .join(CLEANED_REFERENCES)
    .join(TRANSCRIPT_TO_GENE_QUANT)
    .set{
        QUANT_INPUTS
    }


// If we just want to run tertiary using our published quantification results
// from previous runs, we can do that

if ( skipQuantification == 'yes'){

    process spoof_quantify {
        
        executor 'local'

        input:
            set val(tag), val(expName), val(species), val(protocol) from QUANT_INPUTS.map{ r -> tuple( r[0], r[2], r[3], r[4] ) }

        output:
            set val(tag), file("results/*") into QUANT_RESULTS

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
        errorStrategy { task.attempt<=10 ? 'retry' : 'finish' }
        maxRetries 10
        
        input:
            set val(tag), val(isDroplet), val(expName), val(species), val(protocol), file(confFile), file(sdrfFile), file(referenceFasta), file(referenceGtf), val(contaminationIndex), file(transcriptToGene) from SMART_INPUTS
            val flag from INIT_DONE_SMART

        output:
            set val(tag), file ("kallisto") into SMART_KALLISTO_DIRS 
            set val(tag), file ("qc") into SMART_QUANT_QC
            file('quantification.log')    

        """
        submitQuantificationWorkflow.sh 'smart-seq' "$expName" "$species" "$protocol" "$confFile" "$sdrfFile" "$referenceFasta" "$transcriptToGene" "$contaminationIndex" "$enaSshUser"
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
            set val(tag), val(isDroplet), val(expName), val(species), val(protocol), file(confFile), file(sdrfFile), file(referenceFasta), file(referenceGtf), val(contaminationIndex), file(transcriptToGene) from DROPLET_INPUTS
            val flag from INIT_DONE_DROPLET

        output:
            set val(tag), file("alevin") into ALEVIN_DROPLET_DIRS
            file('quantification.log')    

        """
        submitQuantificationWorkflow.sh 'droplet' "$expName" "$species" "$protocol" "$confFile" "$sdrfFile" "$referenceFasta" "$transcriptToGene" "$contaminationIndex" "$enaSshUser"
        """
    }

    // Collect the smart and droplet workflow results

    SMART_KALLISTO_DIRS
        .concat(ALEVIN_DROPLET_DIRS)
        .set { 
            QUANT_RESULTS 
        }
}

// Aggregation just needs the quantifation results and the transcript/ gene
// mappings

TAGS_FOR_AGGR
    .join(QUANT_RESULTS)
    .join(TRANSCRIPT_TO_GENE_AGGR)
    .groupTuple( by: [1,2] )
    .set{ GROUPED_QUANTIFICATION_RESULTS }

// Use existing aggregations if specified

if (skipAggregation == 'yes' ){

    process spoof_aggregate {
    
        executor 'local'
        
        input:
            set val(tags), val(expName), val(species), file('quant_results/??/protocol'), file('quant_results/??/*'), file('quant_results/??/*') from GROUPED_QUANTIFICATION_RESULTS   

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
            set val(tags), val(expName), val(species), file('quant_results/??/protocol'), file('quant_results/??/*'), file('quant_results/??/*') from GROUPED_QUANTIFICATION_RESULTS   

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

COUNT_MATRICES
    .into{
        COUNT_MATRICES_FOR_CELL_TO_LIB
        COUNT_MATRICES_FOR_META_MATCHING
        COUNT_MATRICES_FOR_TERTIARY
        COUNT_MATRICES_FOR_BUNDLE
    }


// Get channel with the first config file for each exp/ species pair. For when
// we need a config but don't need the protocol-specific bits, for example a
// process needs to know the column with technical replicates etc.

TAGS_FOR_TERTIARY                                           // tag, expname, species, protocol
    .join(PRE_CONF_FOR_TERTIARY)                            // tag, expname, species, protocol, conf, confsdrf
    .groupTuple( by: [1,2] )                                // [tag], expname, species, [protocol], [conf], [confsdrf]
    .map{ r -> tuple(r[1], r[2], r[4][0], r[5][0]) }        // expname, species, conf, confsdrf
    .join(COUNT_MATRICES_FOR_CELL_TO_LIB, by: [0,1])        // expname, species, conf, confsdrf, countmatrix
    .join(IS_DROPLET_EXP_SPECIES_FOR_CELLMETA, by: [0,1])   // expname, species, conf, confsdrf, countmatrix, isdroplet
    .set{PRE_TERTIARY_INPUTS} 

// Separate Droplet from non-droplet

DROPLET_CELL_TO_LIB_INPUT = Channel.create()
SMART_CELL_TO_LIB_INPUT = Channel.create()

PRE_TERTIARY_INPUTS.map{r -> tuple(r[0], r[1], r[2], r[4], r[5], file('NO_FILE')) }.choice( DROPLET_CELL_TO_LIB_INPUT, SMART_CELL_TO_LIB_INPUT ) {a -> 
    a[4] == 'True' ? 0 : 1
}

// SMART does not need the cell/lib mapping, so pass the empty value through

SMART_CELL_TO_LIB_INPUT
    .map{ r -> tuple(r[0], r[1], r[5]) }
    .set{SMART_CELL_TO_LIB}
   
// Make a cell-run mapping, required by the SDRF condense process to 'explode'
// the SDRF annotations

process cell_run_mapping {
   
    cache 'deep'
    
    input:
        set val(expName), val(species), file(confFile), file(countMatrix), val(isDroplet), file(emptyFile) from DROPLET_CELL_TO_LIB_INPUT
 
    output:
        set val(expName), val(species), file('cell_to_library.txt') into DROPLET_CELL_TO_LIB
 
    """
    makeCellLibraryMapping.sh $countMatrix $confFile cell_to_library.txt 
    """
}

// Now make a condensed SDRF file. For this to operate correctly with droplet
// data, the cell-library mappings file must be present

SMART_CELL_TO_LIB
    .concat(DROPLET_CELL_TO_LIB)
    .join(META_WITH_SPECIES_FOR_TERTIARY, by: [0,1])
    .set{
       CONDENSE_INPUTS
   }

process condense_sdrf {
        
    cache 'deep'
    
    conda "${baseDir}/envs/atlas-experiment-metadata.yml"

    input:
        set val(expName), val(species), file(cell_to_lib), file(idfFile), file(origSdrfFile), file(cellsFile) from CONDENSE_INPUTS 

    output:
        set val(expName), val(species), file("${expName}.${species}.condensed-sdrf.tsv") into CONDENSED 

    """
    echo -e "exclusions: $ZOOMA_EXCLUSIONS"
    single_cell_condensed_sdrf.sh -e $expName -f $idfFile -o \$(pwd) -z $ZOOMA_EXCLUSIONS
    mv ${expName}.condensed-sdrf.tsv "${expName}.${species}.condensed-sdrf.tsv"
    """        
}

CONDENSED.into{
    CONDENSED_FOR_META
    CONDENSED_FOR_BUNDLE
}

// 'unmelt' the condensed SDRF to get a metadata table to pass for tertiary
// analysis

process unmelt_condensed_sdrf {
        
    conda "${baseDir}/envs/atlas-experiment-metadata.yml"
    
    cache 'deep'

    input:
        set val(expName), val(species), file(condensedSdrf) from CONDENSED

    output:
       set val(expName), val(species), file("${expName}.metadata.tsv") into UNMELTED_META 
        
    """
    unmelt_condensed.R -i $condensedSdrf -o ${expName}.metadata.tsv --retain-types --has-ontology
    """        
}

// Match the cell metadata to the expression matrix 

process match_metadata_to_cells {
    
    cache 'deep'

    input:
        set val(expName), val(species), file(cellMeta), file(countMatrix) from UNMELTED_META.join(COUNT_MATRICES_FOR_META_MATCHING, by: [0,1])

    output:
       set val(expName), val(species), file("${expName}.metadata.matched.tsv") into MATCHED_META 

    """
    matchMetadataToCells.sh $cellMeta $countMatrix ${expName}.metadata.matched.tsv 
    """
}

MATCHED_META.into{
    MATCHED_META_FOR_TERTIARY
    MATCHED_META_FOR_BUNDLE
}

// Tertiary analysis requires the count matrix, references, and droplet status

CONF_REF_BY_EXP_SPECIES_FOR_TERTIARY
    .join(COUNT_MATRICES_FOR_TERTIARY, by: [0,1])
    .join(MATCHED_META_FOR_TERTIARY, by: [0,1])
    .join(IS_DROPLET_EXP_SPECIES_FOR_TERTIARY, by: [0,1])
    .set{TERTIARY_INPUTS}

if ( tertiaryWorkflow == 'scanpy-galaxy' ) {

    // Re-use previous tertiary results

    if (skipTertiary == 'yes' ){

        process spoof_scanpy_galaxy {
        
            executor 'local'
            input:
                set val(expName), val(species), file(referenceFasta), file(referenceGtf), file(contIndex), file(confFile), file(countMatrix), file(cellMetadata), val(isDroplet) from TERTIARY_INPUTS
            
            output:
                set val(expName), val(species), file("matrices/raw_filtered.zip"), file("matrices/filtered_normalised.zip"), file("clusters_for_bundle.txt"), file("umap"), file("tsne"), file("markers"), file('clustering_software_versions.txt') into TERTIARY_RESULTS
        
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
            
            cache 'deep'
            
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
                set val(expName), val(species), file(referenceFasta), file(referenceGtf), file(contIndex), file(confFile), file(countMatrix), file(cellMetadata), val(isDroplet) from TERTIARY_INPUTS

            output:
                set val(expName), val(species), file("matrices/raw_filtered.zip"), file("matrices/filtered_normalised.zip"), file("clusters_for_bundle.txt"), file("umap"), file("tsne"), file("markers"), file('clustering_software_versions.txt') into TERTIARY_RESULTS

            script:
 
                """
                    submitTertiaryWorkflow.sh "$expName" "$species" "$countMatrix" "$referenceGtf" "$cellMetadata" "$isDroplet" "$galaxyCredentials"
                """
        }
    }
} else{

    // Don't do tertiary analysis
    
    process spoof_tertiary {
    
        executor 'local'
        
        input:
            set val(expName), val(species), file(referenceFasta), file(referenceGtf), file(contIndex), file(confFile), file(countMatrix), file(cellMetadata), val(isDroplet) from TERTIARY_INPUTS

        output:
            set val(expName), val(species), file("NOFILT"), file("NONORM"), file("NOCLUST"), file("NOUMAP"), file("NOTSNE"), file("NOMARKERS"), file('NOSOFTWARE') into TERTIARY_RESULTS

        """
            mkdir -p matrices
            cp -p ${countMatrix} matrices 
            touch NOFILT NONORM NOPCA NOCLUST NOUMAP NOTSNE NOMARKERS NOSOFTWARE
        """
    }    
}

// Bundling requires pretty much everything, so we need to gather a few channels

TAGS_FOR_BUNDLE
    .map{r -> tuple(r[1], r[2], r[3])}
    .groupTuple( by: [0,1] )
    .map{r -> tuple(r[0], r[1], r[2].join(","))}
    .set{
        PROTOCOLS_BY_EXP_SPECIES
    }

PROTOCOLS_BY_EXP_SPECIES
    .join(CONF_REF_BY_EXP_SPECIES_FOR_BUNDLE, by: [0,1])
    .join(MATCHED_META_FOR_BUNDLE, by: [0,1])
    .join(CONDENSED_FOR_BUNDLE, by: [0,1])
    .join(COUNT_MATRICES_FOR_BUNDLE,  by: [0,1])
    .join(TPM_MATRICES, remainder: true, by: [0,1])
    .join(TERTIARY_RESULTS, by: [0, 1])
    .set { BUNDLE_INPUTS } 

process bundle {
    
    conda "${baseDir}/envs/nextflow.yml"

    publishDir "$SCXA_RESULTS/$expName/$species", mode: 'copy', overwrite: true
    
    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }
    maxRetries 20
    
    input:
        set val(expName), val(species), val(protocolList), file(referenceFasta), file(referenceGtf), file(contaminationIndex), file(confFile), file(condensedSdrf), file(cellMeta), file(rawMatrix), file(tpmMatrix), file(filteredMatrix), file(normalisedMatrix), file(clusters), file('*'), file('*'), file('*'), file(softwareReport) from BUNDLE_INPUTS
        
    output:
        file('bundle/*')
        file('bundle.log')
        file('bundleLines.txt') into NEW_BUNDLES
        
        """
            submitBundleWorkflow.sh "$expName" "$species" "$protocolList" "$confFile" "$referenceFasta" "$referenceGtf" "$tertiaryWorkflow" "$consensedSdrf" "$cellMeta" "$rawMatrix" "$filteredMatrix" "$normalisedMatrix" "$tpmMatrix" "$clusters" markers tsne $softwareReport
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

