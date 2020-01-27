#!/usr/bin/env bash

workflow=$1
expName=$2
species=$3
protocol=$4
confFile=$5
sdrfFile=$6
referenceFasta=$7
referenceGtf=$8
contaminationIndex=$9
enaSsshUser=$10

# Remove existing results downstrea of this step

for stage in aggregation scanpy bundle; do
    rm -rf $SCXA_RESULTS/$expName/$species/$stage/$protocol
done

# Create necessary output directories

RESULTS_ROOT=$PWD
SUBDIR="$expName/$species/quantification/$protocol"     
mkdir -p $SCXA_WORK/$SUBDIR
mkdir -p $SCXA_NEXTFLOW/$SUBDIR
mkdir -p $SCXA_RESULTS/$SUBDIR/reports

# Workflow-specific parameterisation

contIndex=''
gtfOption=
protocolOption=
enaSshOption=

if [ "$worflow" = 'smart-seq' ]; then

    enaSshOption="--enaSshUser $enaSshUser"

    # Supply a contamination index where provided in the config

    if [ "$contaminationIndex" != 'None' ]; then
        contIndex="--contaminationIndex $contaminationIndex"
    fi 

elif [ "$workflow" = 'droplet' ]; then
    protocolOption="--protocol $protocol"
    gtfOption="--referenceGtf $RESULTS_ROOT/$referenceGtf"
fi

# Submit from the nextflow dir so that logs etc are collected

pushd $SCXA_NEXTFLOW/$SUBDIR > /dev/null

nextflow run \
    -config $RESULTS_ROOT/$confFile \
    --sdrf $RESULTS_ROOT/$sdrfFile \
    --referenceFasta $RESULTS_ROOT/$referenceFasta \
    --resultsRoot $RESULTS_ROOT \
    --manualDownloadFolder $SCXA_DATA/ManuallyDownloaded/$expName \
    $contIndex $gtfOption $enaSshOption \
    -resume \
    $SCXA_WORKFLOW_ROOT/workflow/scxa-workflows/w_${workflow}_quantification/main.nf \
    -work-dir $SCXA_WORK/$SUBDIR \
    -with-report $SCXA_RESULTS/$SUBDIR/reports/report.html \
    -with-trace  $SCXA_RESULTS/$SUBDIR/reports/trace.txt \
    -N $SCXA_REPORT_EMAIL \
    -with-dag $SCXA_RESULTS/$SUBDIR/reports/flowchart.pdf

if [ $? -ne 0 ]; then
    echo "Workflow failed for $expName - $species - scxa-${workflow}-quantification-workflow" 1>&2
    exit 1
fi
            
popd > /dev/null
cp $SCXA_NEXTFLOW/\$SUBDIR/.nextflow.log quantification.log










for stage in aggregation scanpy bundle; do
    rm -rf $SCXA_RESULTS/$expName/$species/\$stage/$protocol
done
RESULTS_ROOT=\$PWD
SUBDIR="$expName/$species/quantification/$protocol"     
mkdir -p $SCXA_WORK/\$SUBDIR
mkdir -p $SCXA_NEXTFLOW/\$SUBDIR
mkdir -p $SCXA_RESULTS/\$SUBDIR/reports
pushd $SCXA_NEXTFLOW/\$SUBDIR > /dev/null
nextflow run \
    -config \$RESULTS_ROOT/$confFile \
    --resultsRoot \$RESULTS_ROOT \
    --sdrf \$RESULTS_ROOT/$sdrfFile \
    --referenceFasta \$RESULTS_ROOT/$referenceFasta \
    --referenceGtf \$RESULTS_ROOT/$referenceGtf \
    --protocol $protocol \
    --manualDownloadFolder $SCXA_DATA/ManuallyDownloaded/$expName \
    -resume \
    $SCXA_WORKFLOW_ROOT/workflow/scxa-workflows/w_droplet_quantification/main.nf \
    -work-dir $SCXA_WORK/\$SUBDIR \
    -with-report $SCXA_RESULTS/\$SUBDIR/reports/report.html \
    -with-trace  $SCXA_RESULTS/\$SUBDIR/reports/trace.txt \
    -N $SCXA_REPORT_EMAIL \
    -with-dag $SCXA_RESULTS/\$SUBDIR/reports/flowchart.pdf
if [ \$? -ne 0 ]; then
    echo "Workflow failed for $expName - $species - scxa-droplet-quantification-workflow" 1>&2
    exit 1
fi

popd > /dev/null
cp $SCXA_NEXTFLOW/\$SUBDIR/.nextflow.log quantification.log






