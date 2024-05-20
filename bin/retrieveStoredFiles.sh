#!/usr/bin/env bash

scriptDir=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $scriptDir/utils.sh

expName=$1
species=$2
storedType=$3
filenames=$4
dest=${5:-'.'}

for file in $filenames; do
    echo "Looking for $file" 1>&2
    storedFile=$SCXA_RESULTS/$expName/$species/$storedType/$file
    fileExists=$(pattern_file_exists "$storedFile")

    if [ $fileExists ]; then
        ln -s $storedFile $dest
    else
        echo "Stored file(s) matching $storedFile does not exist" 1>&2
        exit 1
    fi
done
