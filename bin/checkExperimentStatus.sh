#!/usr/bin/env bash

expName=$1
idfFile=$2
sdrfFile=$3
cellsFile=$4
overwrite=$5

# Start by assuming a new experiment

newExperiment=1
bundleManifests=$(ls $SCXA_RESULTS/$expName/*/bundle/MANIFEST 2>/dev/null || true)

cells_filesize=0
if [ -n "$cellsfile" ] && [ -e $cellsFile ]; then
    cells_filesize=$(stat --printf="%s" $(readlink ${cellsFile}))
fi

# If there are existing bundles for this experiment then just use
# those, unless the related SDRFs have been updated

if [ -n "$bundleManifests" ] && [ "$overwrite" != 'yes' ]; then
    newExperiment=0
    while read -r bundleManifest; do
        if [[ $idfFile -nt "$bundleManifest" || $sdrfFile -nt "$bundleManifest" || ( $cells_filesize -gt 0 && $cellsFile -nt "$bundleManifest" ) ]]; then
            newExperiment=1
            break
        fi
    done <<< "$(echo -e "$bundleManifests")"
fi

# Check for the existence of Atlas prod lock files indicating loading
# in progress.  Only tag new experiments when there are no loading
# locks present

bundleLocks=$(ls $SCXA_RESULTS/$expName/*/bundle/atlas_prod.loading 2>/dev/null || true)

if [ $newExperiment -eq 1 ] && [ -z "$bundleLocks" ]; then
    echo -n 'new'
else
    echo -n 'old'
fi
