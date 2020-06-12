#!/usr/bin/env bash

expName=$1
confFile=$2

# Check if we're looking at a controlled access datasets, and if we have the
# necessaries

controlledAccessField=$(parseNfConfig.py --paramFile $confFile --paramKeys params,fields,controlled_access)

# This controlled access field will only have been set by the config script if
# a) a controlled access field is present in the SDRF and b) it has some 'yes'
# values

if [ "$controlledAccessField" != 'None' ]; then
    if [ -z "$CONTROLLED_ACCESS_DATASETS" ]; then
        echo "$expName is set as controlled access in $confFile, but the CONTROLLED_ACCESS_DATASETS environment variable is is not set" 1>&2
        exit 1
    elif [ ! -e "$CONTROLLED_ACCESS_DATASETS" ]; then
        echo "Specified controlled access datasets file $CONTROLLED_ACCESS_DATASETS does not exist" 1>&2
        exit 1
    else
        datasetInfo=$(grep "^$expName$(printf '\t')" $CONTROLLED_ACCESS_DATASETS)
        if [ $? -ne 0 ]; then
            echo "Can't find controlled access info for $expName in $CONTROLLED_ACCESS_DATASETS" 1>&2
            exit 1
        else
            caDir=$(echo -e "$datasetInfo" | awk '{print $2}')
            caUser=$(echo -e "$datasetInfo" | awk '{print $3}')
            if [ "$USER" != "$caUser" ]; then
                echo "I am $USER, not $caUser, I can't quantify this dataset" 1>&2
                exit 1
            elif [ ! -d "$caDir" ]; then
                echo "Controlled access directory $caDir specified in $CONTROLLED_ACCESS_DATASETS for $expName does not exist"
                exit 1
            fi 
        fi
    fi
    echo -n "yes"
else
    echo -n "no"
fi

