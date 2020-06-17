#!/usr/bin/env bash

# This is a further wrapper layer around submitControlWorkflow.sh. That script
# was initially designed to submit a single Nextflow workflow controlling
# analysis of all or a batch of SCXA studies. However it can submit the control
# workflow for one experiment at a time, with independent working directories
# etc, and this may allow for more effective use of resources than Nextflow
# processes waiting for a batch to finish before submitting more.

usage() { echo "Usage: $0 [-q <re-use existing quantifications where present, yes or no>] [-a <re-use existing aggregations where present, yes or no>] [-t <tertiary workflow>] [-u <re-use existing tertiary results where present, yes or no>] [-w <overwrite exising results, yes or no>] [-m <maximum number of experiments to analyse at once>]"  1>&2; }  

q=no
a=no
t=scanpy-galaxy
u=no
w=no
m=5

while getopts ":q:a:t:u:w:m:" o; do
    case "${o}" in
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
        m)
            m=${OPTARG}
            ;;
        *)
            usage
            exit 0
            ;;
    esac
done
shift $((OPTIND-1))

maxExpsToRun=$m

workflow=scxa-control-workflow
controlJobSuffix=${SCXA_ENV}_$workflow

export SCXA_RESULTS=$SCXA_WORKFLOW_ROOT/results
submissionMarker="$SCXA_RESULTS/.submitting"

if [ -e  $submissionMarker ]; then
    echo "Submission process ongoing"
else
    touch $submissionMarker
fi

currentJobs=$(bjobs -w | grep $controlJobSuffix)
nJobsRunning=$(echo $currentJobs | wc -l)
maxToSubmit=$((maxExpsToRun-nJobsRunning))

if [ $nJobsRunning -ge $maxExpsToRun ]; then
    echo "$nJobsRunning already running, submitting no more\n"
    exit 0
else
    echo "Will sumbit a maximum of $maxToSubmit jobs"
fi

submitted=0

# Fetch the Git SDRFs

if [ ! -d 'metadata' ]; then
    git clone $SCXA_METADATA_REPO metadata
fi

pushd metadata > /dev/null
git pull > /dev/null
popd > /dev/null

# Submit workflows for any updatated metdata files. The Nextflow workflow
# itself contains logic for not re-analaysing where it does not need to.

export PATH=$SCXA_WORKFLOW_ROOT/workflow/scxa-control-workflow/bin:$PATH

while read -r idfFile; do
    idfFileName=$(basename $idfFile)
    expId=$(echo $idfFileName | awk -F'.' '{print $1}')
    sdrfFile=$(dirname $idfFile)/${expId}.sdrf.txt
    
    grep -P "$expId\\t" $SCXA_RESULTS/excluded.txt > /dev/null
    if [ $? -eq 0 ]; then
        experimentStatus='excluded'
    else 
        # Notes: checkExperimentStatus.sh is the same as that used by the
        # Nextflow workflow, so the logic is kept consistent   
        experimentStatus=$(checkExperimentStatus.sh $expId $sdrfFile $overwrite)
    fi

    echo "Experiment status for $expId: $experimentStatus"
    currentlyRunning='false'

    # Check if the jobs is already running

    echo -e "$currentJobs" | grep "${expId}_$controlJobSuffix" > /dev/null
    if [ $? -eq 0 ]; then
        currentlyRunning='true'
        echo "$expId is currently running"
    fi

    if [ "$currentlyRunning" = 'false' ] && [ $submitted -lt $maxToSubmit ] && [ "$experimentStatus" == 'new' ]; then

        echo "Submitting $expId for re-analysis"
        $SCXA_WORKFLOW_ROOT/workflow/scxa-control-workflow/bin/submitControlWorkflow.sh -e $expId -q $q -t $t -a $a -u $u -w $w 

        currentJobs=$(bjobs -w | grep $controlJobSuffix)
        submitted=$((submitted+1))

    fi
done <<< "$(ls $SCXA_WORKFLOW_ROOT/metadata/*/*/*.idf.txt)"

ls $SCXA_WORKFLOW_ROOT/results/*/*/bundle/MANIFEST | awk -F'/' '{OFS="\t";} {print $(NF-2),$(NF-3),$0}' > $SCXA_RESULTS/all.done.txt

rm $submissionMarker
