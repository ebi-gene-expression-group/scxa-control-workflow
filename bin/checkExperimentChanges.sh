#!/usr/bin/env bash

# This script prints a status that will dictate how the workflow routes
# analysis. 

scriptDir=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $scriptDir/utils.sh

expName=$1
species=$2
newConfig=$3
newMetaForQuant=$4
newMetaForTertiary=$5
skipQuantification=$6
skipAggregation=$7
skipTertiary=$8
overWrite=$9

# If there's already a MANIFEST then this is a potential re-run. Only do this
# if the config changes impact on analysis. This will be the case if there's a
# difference between new derived files and the old ones.

resultsDir=$SCXA_RESULTS/$expName/$species
manifestPath=$resultsDir/bundle/MANIFEST

# Files we use to check the presence of results for individual stages

quantExists=$(pattern_file_exists "$resultsDir/quantification/*/quantification.log")
aggExists=$(pattern_file_exists "$resultsDir/aggregation/matrices/counts_mtx.zip")
metaExists=$(pattern_file_exists "$resultsDir/metadata/${expName}-${species}.metadata.matched.tsv")
tertiaryExists=$(pattern_file_exists "$resultsDir/scanpy/clustering_software_versions.txt")
bundleExists=$(pattern_file_exists $resultsDir/bundle/MANIFEST)

# Default is 'unchanged', which will trigger no analysis. Explicity overwrites,
# changed config or missing data will trigger analysis processes

experimentStatus='unchanged'

# If overwrite is specified, then re-analyse everything except the explicit skips

if [ "$overWrite" == 'yes' ]; then
    echo "Overwriting" 1>&2

    if [ "$skipTertiary" = 'yes' ] && [ $tertiaryExists ]; then
        experimentStatus='changed_for_bundling'
    elif [ "$skipAggregation" = 'yes' ] && [ $aggExists ]; then
        experimentStatus='changed_for_tertiary'
    elif [ "$skipQuantification" = 'yes' ] && [ $quantExists ]; then
        experimentStatus='changed_for_aggregation'
    else 
        experimentStatus='changed_for_quantification'
    fi

else
    
    echo "Overwrite not set, checking existing config" 1>&2

    # Where overwrite was not specified, we'll look at the changes to see what need re-doing
   
    configDiff=0
    quantMetaDiff=0
    tertiaryMetaDiff=0 

    # Has config changed?

    while read -r nc; do
        oldConfig=$SCXA_CONF/study/$(basename $nc)
        
        diff_configs $newConfig/$nc $oldConfig
        if [ $? -eq 1 ]; then
            echo "$newConfig/$nc is different to $oldConfig" 1>&2
            configDiff=1
            break
        else
            echo "$newConfig/$nc is not different to $oldConfig" 1>&2
        fi
    done <<< "$(ls $newConfig)"

    # Has quantification-relevant metadata changed?
    
    if [ $quantExists ]; then
        while read -r nmfq; do
            oldMetaForQuant=$SCXA_CONF/study/$(basename $nmfq)
            
            diff_configs $newMetaForQuant/$nmfq $oldMetaForQuant
            if [ $? -eq 1 ]; then
                echo "$newMetaForQuant/$nmfq is different to $oldMetaForQuant" 1>&2
                quantMetaDiff=1
                break
            else
                echo "$newMetaForQuant/$nmfq is not different to $oldMetaForQuant" 1>&2
            fi
        done <<< "$(ls $newMetaForQuant)"

        # Has tertiary-relevant metadata changed?

        if [ $tertiaryExists ]; then
            while read -r nmft; do
                oldMetaForTertiary=$SCXA_CONF/study/$(basename $nmft)
                
                # Compare configs and derived SDRFS

                diff_configs $newMetaForTertiary/$nmft $oldMetaForTertiary
                if [ $? -eq 1 ]; then
                    echo "$newMetaForTertiary/$nmft is different to $oldMetaForTertiary" 1>&2
                    tertiaryMetaDiff=1
                    break
                else
                    echo "$newMetaForTertiary/$nmft is not different to $oldMetaForTertiary" 1>&2
                fi
            done <<< "$(ls $newMetaForTertiary)"
        else
            echo 'No completed tertiary results were found' 1>&2
        fi
    else
        echo 'No completed quantification results were found' 1>&2
    fi
    
    if [ $quantMetaDiff -eq 1 ] || [ ! $quantExists ]; then
        experimentStatus='changed_for_quantification'
    elif [ ! $aggExists ] || [ ! $metaExists ] || [ $tertiaryMetaDiff -eq 1 ]; then
        experimentStatus='changed_for_aggregation'
    elif [ ! $tertiaryExists ]; then
        experimentStatus='changed_for_tertiary'
    elif [ ! $bundleExists ]; then
        experimentStatus='changed_for_bundling'
    fi
fi    

echo -n "$experimentStatus"
