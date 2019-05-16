#!/usr/bin/env bash

usage() { echo "Usage: $0 [-e <experiment ID>] [-q <re-use existing quantifications where present, yes or no>] [-a <re-use existing aggregations where present, yes or no>] [-t <tertiary workflow>] [-u <re-use existing tertiary results where present, yes or no>] [-w <overwrite exising results, yes or no>]"  1>&2; }  

e=
s=no
t=none
u=no
w=no

while getopts ":e:q:a:t:u:w:" o; do
    case "${o}" in
        e)
            e=${OPTARG}
            ;;
        q)
            q=${OPTARG}
            ;;
        a)
            a=${OPTARG}
            ;;
        t)
            t=${OPTARG}
            ;;
        u)
            u=${OPTARG}
            ;;
        w)
            w=${OPTARG}
            ;;
        *)
            usage
            exit 0
            ;;
    esac
done
shift $((OPTIND-1))

# Simple script to trigger SCXA workflow

expName=$e
skipQuantification=$q
skipAggregation=$a
tertiaryWorkflow=$t
skipTertiary=$u
overwrite=$w

workflow=scxa-control-workflow

# Change working dir for experiment-specific runs

workingDir="$SCXA_WORKFLOW_ROOT/work/${workflow}"
if [ -n "$expName" ]; then
    workingDir="${workingDir}_$expName"
fi

cd $SCXA_WORKFLOW_ROOT

# Clone or update the scxa-workflows repo containing config files etc

if [ ! -d 'workflow/scxa-workflows' ]; then
    git clone --recursive https://github.com/ebi-gene-expression-group/scxa-workflows workflow/scxa-workflows
fi

pushd workflow/scxa-workflows > /dev/null
git checkout $SCXA_BRANCH > /dev/null
git pull > /dev/null
git submodule update > /dev/null
popd > /dev/null

# Add scripts from the scxa-workflows repo to the PATH

export PATH=$(pwd)/workflow/scxa-workflows/bin:$PATH

# Build the Nextflow command

expNamePart=
if [ -n "$expName" ]; then
    expNamePart="--expName $expName"
fi

skipQuantificationPart=
if [ -n "$skipQuantification" ]; then
    if [ "$skipQuantification" != 'yes' ] && [ "$skipQuantification" != 'no' ]; then
        echo "Skip quantification (-q) must be 'yes' or 'no', $skipQuantification provided." 1>&2
        exit 1
    fi

    skipQuantificationPart="--skipQuantification $skipQuantification"
fi

skipAggregationPart=
if [ -n "$skipAggregation" ]; then
    if [ "$skipAggregation" != 'yes' ] && [ "$skipAggregation" != 'no' ]; then
        echo "Skip aggregation (-a) must be 'yes' or 'no', $skipAggregation provided." 1>&2
        exit 1
    fi

    skipAggregationPart="--skipAggregation $skipAggregation"
fi

skipTertiaryPart=
if [ -n "$skipTertiary" ]; then
    if [ "$skipTertiary" != 'yes' ] && [ "$skipTertiary" != 'no' ]; then
        echo "Skip tertiary (-u) must be 'yes' or 'no', $skipTertiary provided." 1>&2
        exit 1
    fi

    skipTertiaryPart="--skipTertiary $skipTertiary"
fi

overwritePart=
if [ -n "$overwrite" ]; then
    if [ "$overwrite" != 'yes' ] && [ "$overwrite" != 'no' ]; then
        echo "Overwrite (-w) must be 'yes' or 'no', $overwrite provided" 1>&2
        exit 1
    fi

    overwritePart="--overwrite $overwrite"
fi

tertiaryWorkflowPart=
galaxyCredentialsPart=

if [ -n "$tertiaryWorkflow" ]; then
    tertiaryWorkflowPart="--tertiaryWorkflow $tertiaryWorkflow"

    if [ "$tertiaryWorkflow" == 'scanpy-galaxy' ]; then
        if [ -z "$GALAXY_CREDENTIALS" ]; then
            echo "Please set the GALAXY_CREDENTIALS environment variable"
            exit 1
        fi
        galaxyCredentialsPart="--galaxyCredentials $GALAXY_CREDENTIALS"
    fi
fi

nextflowCommand="nextflow run -N $SCXA_REPORT_EMAIL -resume workflow/${workflow}/main.nf $expNamePart $skipQuantificationPart $skipAggregationPart $tertiaryWorkflowPart $skipTertiaryPart $galaxyCredentialsPart $overwritePart --enaSshUser fg_atlas_sc --sdrfDir $SCXA_SDRF_DIR -work-dir $workingDir"

# Run the LSF submission if it's not already running

bjobs -w | grep "${SCXA_ENV}_$workflow" > /dev/null 2>&1

if [ $? -ne 0 ]; then

    # If workflow completed successfully we can clean up the work dir. If not,
    # then the caching from the work dir will be useful to resume

    successMarker="$SCXA_WORKFLOW_ROOT/work/.success"

    if [ -e "$successMarker" ]; then
        echo "Previous run succeeded, cleaning up $workingDir"
        mv $workingDir ${workingDir}_removing_$$ 
        nohup rm -rf ${workingDir}_removing_$$ &
        rm -rf $successMarker
    fi

    echo "Submitting job"
    rm -rf run.out run.err .nextflow.log*  
    bsub \
        -J ${SCXA_ENV}_$workflow \
        -M 4096 -R "rusage[mem=4096]" \
        -u $SCXA_REPORT_EMAIL \
        -o run.out \
        -e run.err \
        "$nextflowCommand" 
else
    echo "Workflow process already running" 
fi   
