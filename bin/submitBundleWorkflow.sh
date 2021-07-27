#!/usr/bin/env bash

expName=$1
species=$2
protocolList=$3
confFile=$4
referenceFasta=$5
referenceGtf=$6
geneMetadata=$7
tertiaryWorkflow=$8
condensedSdrf=$9
cellMeta=${10}
rawMatrix=${11}
filteredMatrix=${12}
normalisedMatrix=${13}
tpmMatrix=${14}
clusters=${15}
markersDir=${16}
tsneDir=${17}
umapDir=${18}
softwareReport=${19}
projectFile=${20}

RESULTS_ROOT=$PWD
SUBDIR="$expName/$species/bundle"     
TPM_OPTIONS=''

tpm_filesize=$(stat --printf="%s" $(readlink ${tpmMatrix}))
if [ "$tpmMatrix" != 'null' ] && [ $tpm_filesize -gt 0 ]; then
    TPM_OPTIONS="--tpmMatrix ${tpmMatrix}"
fi

TERTIARY_OPTIONS=''
if [ "$tertiaryWorkflow" == 'scanpy-galaxy' ]; then
    TERTIARY_OPTIONS="--tertiaryWorkflow $tertiaryWorkflow --rawFilteredMatrix ${filteredMatrix} --normalisedMatrix ${normalisedMatrix} --clusters ${clusters} --tsneDir $tsneDir --umapDir $umapDir --markersDir $markersDir --tertiarySoftwareReport ${softwareReport} --projectFile ${projectFile}"
fi 

mkdir -p $SCXA_WORK/$SUBDIR
mkdir -p $SCXA_NEXTFLOW/$SUBDIR
mkdir -p $SCXA_RESULTS/$SUBDIR/reports
pushd $SCXA_NEXTFLOW/$SUBDIR > /dev/null

# Note: we could have received multiple config files, coming from multiple
# protocols- just pass the first, we don't have any protocol-specific
# behaviours in bundling

nextflow run \
    -config $RESULTS_ROOT/$(echo $confFile | awk '{print $1}') \
    --masterWorkflow scxa-control-workflow \
    --expName $expName \
    --resultsRoot $RESULTS_ROOT \
    --protocolList ${protocolList} \
    --cellMetadata ${cellMeta} \
    --condensedSdrf ${condensedSdrf} \
    --rawMatrix ${rawMatrix} $TPM_OPTIONS \
    --referenceFasta $referenceFasta \
    --geneMetadata $geneMetadata \
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
