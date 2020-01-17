#!/usr/bin/env bash

scriptDir=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $scriptDir/utils.sh

expName=$1
species=$2
newConfig=$3
newSdrf=$4

# If there's already a MANIFEST then this is a potential re-run. Only do this
# if the config changes impact on analysis. This will be the case if there's a
# difference between new derived files and the old ones.

manifestPath=$SCXA_RESULTS/$expName/$species/bundle/MANIFEST
oldConfig=$SCXA_CONF/study/$(basename $newConfig)
oldSdrf=$SCXA_CONF/study/$(basename $newSdrf)

experimentStatus='changed'

if [ -e "$manifestPath" ]; then
    
    # Compare configs and derived SDRFS

    diff_configs $newConfig $oldConfig
    configDiff=$?

    diff_configs $newSdrf $oldSdrf
    sdrfDiff=$?

    if [ $configDiff -eq 0 ] && [ $sdrfDiff -eq 0 ]; then
        experimentStatus='unchanged'
    fi
fi

echo -n "$experimentStatus"
