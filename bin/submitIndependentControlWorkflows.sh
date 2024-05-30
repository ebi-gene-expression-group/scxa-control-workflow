#!/usr/bin/env bash

# This is a further wrapper layer around submitControlWorkflow.sh. That script
# was initially designed to submit a single Nextflow workflow controlling
# analysis of all or a batch of SCXA studies. However it can submit the control
# workflow for one experiment at a time, with independent working directories
# etc, and this may allow for more effective use of resources than Nextflow
# processes waiting for a batch to finish before submitting more.

usage() { echo "Usage: $0 [-q <re-use existing quantifications where present, yes or no>] [-a <re-use existing aggregations where present, yes or no>] [-u <re-use existing tertiary results where present, yes or no>] [-w <overwrite exising results, yes or no>] [-m <maximum number of experiments to analyse at once>]"  1>&2; }  

q=no
a=no
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
submitNew='yes'

workflow=scxa-control-workflow
controlJobSuffix=${SCXA_ENV}_$workflow

export SCXA_RESULTS=$SCXA_WORKFLOW_ROOT/results
submissionMarker="$SCXA_RESULTS/.submitting"

if [ -e  $submissionMarker ]; then
    echo "Submission process ongoing"
    exit 0
else
    touch $submissionMarker
fi
# this needs migration
currentJobs=$(bjobs -w | grep $controlJobSuffix)
nJobsRunning=$(echo -e "$currentJobs" | wc -l)
maxToSubmit=$((maxExpsToRun-nJobsRunning))

if [ $nJobsRunning -ge $maxExpsToRun ]; then
    echo -e "$nJobsRunning already running, submitting no more\n"
    submitNew='no'
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

echo "Checking experiments for submission"
while read -r idfFile; do
    idfFileName=$(basename $idfFile)
    expId=$(echo $idfFileName | awk -F'.' '{print $1}')
    sdrfFile=$(dirname $idfFile)/${expId}.sdrf.txt
    cellsFile=$(dirname $idfFile)/${expId}.cells.txt
    runLog=$SCXA_WORKFLOW_ROOT/nextflow/${expId}_${SCXA_ENV}_scxa-control-workflow/run.out   

    experimentStatus='old'

    grep -P "$expId\\t" $SCXA_RESULTS/excluded.txt > /dev/null
    if [ $? -eq 0 ]; then
        experimentStatus='excluded'
    else
        if [ -e "$runLog" ]; then
            cat $runLog | grep "Exited with" > /dev/null
            if [ $? -eq 0 ]; then
                echo -e "${expId}\tUnknown error: see <workflows root>/nextflow/${expId}_${SCXA_ENV}_scxa-control-workflow/run.out and delete when resolved" >> $SCXA_RESULTS/excluded.txt
                experimentStatus='excluded'
            fi
        fi 
 
    fi

    if [ "$experimentStatus" != 'excluded' ]; then
      
        # Notes: checkExperimentStatus.sh is the same as that used by the
        # Nextflow workflow, so the logic is kept consistent   
        experimentStatus=$(checkExperimentStatus.sh $expId $idfFile $sdrfFile $cellsFile $overwrite)
    fi    

    echo "Experiment status for $expId: $experimentStatus"

    if [ "$submitNew" == 'yes' ]; then

        # Check if the jobs is already running

        currentlyRunning='false'
        echo -e "$currentJobs" | grep "${expId}_$controlJobSuffix" > /dev/null
        if [ $? -eq 0 ]; then
            currentlyRunning='true'
            echo "$expId is currently running"
        fi

        if [ "$currentlyRunning" = 'false' ] && [ $submitted -lt $maxToSubmit ] && [ "$experimentStatus" == 'new' ]; then

            echo "Submitting $expId for re-analysis"
            $SCXA_WORKFLOW_ROOT/workflow/scxa-control-workflow/bin/submitControlWorkflow.sh -e $expId -q $q -a $a -u $u -w $w 

            # this needs migration
            currentJobs=$(bjobs -w | grep $controlJobSuffix)
            submitted=$((submitted+1))

        fi
    fi

done <<< "$(ls $SCXA_WORKFLOW_ROOT/metadata/*/*/*.idf.txt | grep -v ANND)"

# List all finished bundles for loading, exluding anything that might have been
# added to the excluded set after completion

echo "Getting completed result set..."
ls $SCXA_WORKFLOW_ROOT/results/*/*/bundle/MANIFEST | while read -r l; do dirname $l; done | awk -F'/' '{OFS="\t";} {print $(NF-2),$(NF-1),$0}' | while read -r m; do
    exp=$(echo -e "$m" | awk '{print $1}');
    grep "^$exp$(printf '\t')" $SCXA_RESULTS/excluded.txt > /dev/null;
    if [ $? -ne 0 ]; then
        echo -e "$m";
    fi;
done > $SCXA_RESULTS/all.done.txt

# Now do all the cleanup for completed experiments that couldn't be done from
# within the workflow. We'll have to tweak this if we ever allow multi-species
# experiments, since the dirs are not currently tagged by species.

echo "Cleaning up..."
# this needs migration
bjobsOutput=$(bjobs -w)

cat $SCXA_RESULTS/all.done.txt | while read -r l; do
    expId=$(echo "$l" | awk '{print $1}')
    species=$(echo "$l" | awk '{print $2}')

    controlJobName="${expId}_${controlJobSuffix}"

    echo -e "$bjobsOutput" | grep " ${controlJobName}" > /dev/null 2>&1

    if [ $? -ne 0 ]; then
        for tmpWfDir in work/scxa-control-workflow_$expId work/$expId nextflow/${expId}_${SCXA_ENV}_scxa-control-workflow nextflow/${expId}; do
            if [ -d "${SCXA_WORKFLOW_ROOT}/${tmpWfDir}" ]; then
                nohup rm -rf ${SCXA_WORKFLOW_ROOT}/${tmpWfDir} &
            fi
        done
    fi
done

rm $submissionMarker
