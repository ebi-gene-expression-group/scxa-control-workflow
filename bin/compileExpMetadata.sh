#!/usr/bin/env bash

# Check a candidate experiment for processing. Experiments are processed if
# they're currently unknown, or have updated SDRF. Experiments are not updated
# if they were processed previously and metadata has not been updated, or if
# the bundle is locked for loading.

scriptDir=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $scriptDir/utils.sh

expName=$1
idfFile=$2
overwrite=$3
force=$4 # If a non-empty string, this will override the placement of an experiment in the excluded file

# Have we excluded this study?

set +e
grep -P "$expName\\t" $SCXA_RESULTS/excluded.txt > /dev/null
excludeStatus=$?
set -e

# This study is in the excluded file, which is why we don't have results.
# Signal to ignore it.

if [ $excludeStatus -eq 0 ] && [ -z "$force" ]; then
    exit 0
fi

if [ ! -e "$idfFile" ]; then
    die "IDF $idfFile file missing"
fi

idfDir=$(dirname $(readlink $idfFile))

sdrfLine=$(grep -iP "^SDRF File\t" $idfFile)
if [ $? -ne 0 ]; then
    die "No SDRF definition found in IDF file" 1
fi

sdrfFileName=$(echo -e "$sdrfLine" | awk -F"\t" '{print $2}')
sdrfFile=${idfDir}/$sdrfFileName

if [ ! -e "$sdrfFile" ]; then
    die "SDRF file $sdrfFile does not exist" 1
else

    # Theoretically the SDRF could have a different name to the IDF and we just
    # have to check the name of the SDRF in the IDF. But let's keep things
    # simple for ourselves and name them the same for downstream processing
    
    ln -s $sdrfFile ${expName}.sdrf.txt
fi

# If curators have provided a .cells.txt file...

cellsFile=${idfDir}/${expName}.cells.txt

if [ -e "$cellsFile" ]; then
    ln -s $cellsFile .
fi
