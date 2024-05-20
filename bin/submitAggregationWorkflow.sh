#!/usr/bin/env bash

expName=$1
species=$2

# Run aggregation with
# https://github.com/ebi-gene-expression-group/scxa-aggregation-workflow

# Clean up downstream results

for stage in scanpy bundle; do
    rm -rf $SCXA_RESULTS/$expName/$species/$stage
done

# Create necessary output directories

RESULTS_ROOT=$PWD
SUBDIR="$expName/$species/aggregation"     
mkdir -p $SCXA_WORK/$SUBDIR
mkdir -p $SCXA_NEXTFLOW/$SUBDIR
mkdir -p $SCXA_RESULTS/$SUBDIR/reports

# Submit from the nextflow dir so that logs etc are collected

pushd $SCXA_NEXTFLOW/$SUBDIR > /dev/null

# If we have a species-wise config, supply to aggregation

species_conf=$SCXA_PRE_CONF/reference/${species}.conf

opt_conf=
if  [ -e $species_conf ]; then
    opt_conf=" --config $species_conf"
fi            

nextflow run $opt_conf \
    --resultsRoot $RESULTS_ROOT \
    --quantDir $RESULTS_ROOT/quant_results \
    -resume \
    $SCXA_WORKFLOW_ROOT/workflow/scxa-workflows/w_aggregation/main.nf \
    -work-dir $SCXA_WORK/$SUBDIR \
    -with-report $SCXA_RESULTS/$SUBDIR/reports/report.html \
    -with-trace  $SCXA_RESULTS/$SUBDIR/reports/trace.txt \
    -N $SCXA_REPORT_EMAIL \
    -with-dag $SCXA_RESULTS/$SUBDIR/reports/flowchart.pdf

if [ $? -ne 0 ]; then
    echo "Workflow failed for $expName - $species - scxa_aggregation_workflow" 1>&2
    exit 1
fi

popd > /dev/null
cp $SCXA_NEXTFLOW/$SUBDIR/.nextflow.log matrices/aggregation.log
