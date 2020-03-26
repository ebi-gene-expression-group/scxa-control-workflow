#!/usr/bin/env bash

expName=$1
species=$2
protocolList=$3
confFile=$4
referenceFasta=$5
referenceGtf=$6
tertiaryWorkflow=$7
cellMeta=$8
rawMatrix=$9
filteredMatrix=$10
normalisedMatrix=${11}
tpmMatrix=${12}
clusters=${13}
markersDir=${14}
tsneDir=${15}
softwareReport=${16}

RESULTS_ROOT=$PWD
SUBDIR="$expName/$species/bundle"     
TPM_OPTIONS=''

tpm_filesize=$(stat --printf="%s" $(readlink ${tpmMatrix}))
if [ "$tpmMatrix" != 'null' ] && [ $tpm_filesize -gt 0 ]; then
    TPM_OPTIONS="--tpmMatrix ${tpmMatrix}"
fi

TERTIARY_OPTIONS=''
if [ "$tertiaryWorkflow" == 'scanpy-galaxy' ]; then
    TERTIARY_OPTIONS="--tertiaryWorkflow $tertiaryWorkflow --rawFilteredMatrix ${filteredMatrix} --normalisedMatrix ${normalisedMatrix} --clusters ${clusters} --tsneDir $tsneDir --markersDir $markersDir --tertiarySoftwareReport ${softwareReport}"
fi 

mkdir -p $SCXA_WORK/$SUBDIR
mkdir -p $SCXA_NEXTFLOW/$SUBDIR
mkdir -p $SCXA_RESULTS/$SUBDIR/reports
pushd $SCXA_NEXTFLOW/$SUBDIR > /dev/null

nextflow run \
    -config $RESULTS_ROOT/$confFile \
    --masterWorkflow scxa-control-workflow \
    --resultsRoot $RESULTS_ROOT \
    --protocolList ${protocolList} \
    --cellMetadata ${cellMeta} \
    --rawMatrix ${rawMatrix} $TPM_OPTIONS \
    --referenceFasta $referenceFasta \
    --referenceGtf $referenceGtf $TERTIARY_OPTIONS \
    -resume \
    $SCXA_WORKFLOW_ROOT/workflow/scxa-workflows/w_bundle/main.nf \
    -work-dir $SCXA_WORK/$SUBDIR \
    -with-report $SCXA_RESULTS/$SUBDIR/reports/report.html \
    -with-trace  $SCXA_RESULTS/$SUBDIR/reports/trace.txt \
    -N $SCXA_REPORT_EMAIL \
    -with-dag $SCXA_RESULTS/$SUBDIR/reports/flowchart.pdf

if [ $? -ne 0 ]; then
    echo "Workflow failed for $expName - $species - scxa-bundle-workflow" 1>&2
    exit 1
fi

popd > /dev/null
cp $SCXA_NEXTFLOW/$SUBDIR/.nextflow.log bundle.log
chmod g+rwx $SCXA_RESULTS/$expName/$species/bundle
echo -e "$expName\t$species\t$SCXA_RESULTS/$expName/$species/bundle" > bundleLines.txt
