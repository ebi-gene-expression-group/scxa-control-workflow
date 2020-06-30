#!/usr/bin/env bash

scriptDir=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $scriptDir/utils.sh

expName=$1
species=$2
newConfig=$3
newSdrf=$4
newCells=$5

# If there's already a MANIFEST then this is a potential re-run. Only do this
# if the config changes impact on analysis. This will be the case if there's a
# difference between new derived files and the old ones.

manifestPath=$SCXA_RESULTS/$expName/$species/bundle/MANIFEST
oldConfig=$SCXA_CONF/study/$(basename $newConfig)
oldSdrf=$SCXA_CONF/study/$(basename $newSdrf)
oldCells=$SCXA_CONF/study/$(basename $newCells)

experimentStatus='unchanged'

if [ -e "$manifestPath" ]; then
    
    # Compare configs and derived SDRFS

    diff_configs $newConfig $oldConfig
    configDiff=$?

    diff_configs $newSdrf $oldSdrf
    sdrfDiff=$?
    
    diff_configs $newCells $oldCells
    cellsDiff=$?

    if [ $configDiff -eq 1 ] || [ $sdrfDiff -eq 1 ]; then
        experimentStatus='changed'
    elif [ $cellsDiff -eq 1 ]; then
        experiment_status='meta_changed'
    fi
else
    experimentStatus="changed"
fi

echo -n "$experimentStatus"
