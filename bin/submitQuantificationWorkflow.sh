#!/usr/bin/env bash

workflow=$1
expName=$2
species=$3
protocol=$4
confFile=$5
sdrfFile=$6
referenceFasta=$7
transcriptToGene=$8
contaminationIndex=$9
enaSshUser=$10

# Remove existing results downstrea of this step

for stage in aggregation scanpy bundle; do
    rm -rf $SCXA_RESULTS/$expName/$species/$stage/$protocol
done

# Create necessary output directories

RESULTS_ROOT=$PWD
SUBDIR="$expName/$species/quantification/$protocol"     

# we ran checks earlier to check that we had all we needed in the case of
# controlled access analysis- so here we just need to retrieve the directory to
# use

quantWorkDir=$SCXA_WORK/$SUBDIR

datasetInfo=$(grep "^$expName\t" $CONTROLLED_ACCESS_DATASETS)
if [ $? -ne 0 ]; then
  caDir=$(echo -e "$datasetInfo" | awk '{print $2}')
  quantWorkDir=$caDir/analysis 
fi

mkdir -p $quantWorkDir
mkdir -p $SCXA_NEXTFLOW/$SUBDIR
mkdir -p $SCXA_RESULTS/$SUBDIR/reports

# Workflow-specific parameterisation

contIndex=''
enaSshOption=

if [ "$worflow" = 'smart-seq' ]; then

    enaSshOption="--enaSshUser $enaSshUser"

    # Supply a contamination index where provided in the config

    if [ "$contaminationIndex" != 'None' ]; then
        contIndex="--contaminationIndex $contaminationIndex"
    fi 
fi

# Submit from the nextflow dir so that logs etc are collected

pushd $SCXA_NEXTFLOW/$SUBDIR > /dev/null

nextflow run \
    -config $RESULTS_ROOT/$confFile \
    --sdrf $RESULTS_ROOT/$sdrfFile \
    --protocol $protocol \
    --referenceFasta $RESULTS_ROOT/$referenceFasta \
    --resultsRoot $RESULTS_ROOT \
    --manualDownloadFolder $SCXA_DATA/ManuallyDownloaded/$expName \
    --transcriptToGene $RESULTS_ROOT/$transcriptToGene \
    $contIndex $enaSshOption \
    -resume \
    $SCXA_WORKFLOW_ROOT/workflow/scxa-workflows/w_${workflow}_quantification/main.nf \
    -work-dir $quantWorkDir \
    -with-report $SCXA_RESULTS/$SUBDIR/reports/report.html \
    -with-trace  $SCXA_RESULTS/$SUBDIR/reports/trace.txt \
    -N $SCXA_REPORT_EMAIL \
    -with-dag $SCXA_RESULTS/$SUBDIR/reports/flowchart.pdf

if [ $? -ne 0 ]; then
    echo "Workflow failed for $expName - $species - scxa-${workflow}-quantification-workflow" 1>&2
    exit 1
fi
            
popd > /dev/null
cp $SCXA_NEXTFLOW/$SUBDIR/.nextflow.log quantification.log
