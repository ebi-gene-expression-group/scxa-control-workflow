#!/usr/bin/env bash

scriptDir=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $scriptDir/utils.sh

# Submit tertiary workflow

expName=$1
species=$2
confFile=$3
countMatrix=$4
geneMetadata=$5
cellMetadata=$6
isDroplet=$7

rm -rf $SCXA_RESULTS/$expName/$species/bundle

# Nextflow subworkflow wants the matrix components separately

zipdir=$(unzip -qql ${countMatrix} | head -n1 | tr -s ' ' | cut -d' ' -f5- | sed 's|/||')
unzip ${countMatrix} 
       
gzip ${zipdir}/matrix.mtx
gzip ${zipdir}/genes.tsv
gzip ${zipdir}/barcodes.tsv

export species=$species
export expName=$expName
export gene_meta_file=$geneMetadata
export matrix_file=${zipdir}/matrix.mtx.gz
export genes_file=${zipdir}/genes.tsv.gz
export barcodes_file=${zipdir}/barcodes.tsv.gz
export cell_meta_file=$cellMetadata
export tpm_filtering='False'
export create_conda_env=no
export SCXA_WORKDIR=$SCXA_WORK
export SCXA_OUTDIR="."
       
if [ "$isDroplet" = 'True' ]; then
    export FLAVOUR=w_droplet_clustering
else
    export FLAVOUR=w_smart-seq_clustering
fi

# Extract things we need from the conf file

cellTypeField=$(parseNfConfig.py --paramFile $confFile --paramKeys params,fields,cell_type)
if [ "$cellTypeField" != 'None' ]; then
    export cell_type_field=$(sanitise_field "$cellTypeField")
    echo "Cell type field: $cell_type_field"
fi

batchField=$(parseNfConfig.py --paramFile $confFile --paramKeys params,fields,batch)
echo "confFile $confFile"
echo "batchField $batchField"
if [ "$batchField" != 'None' ]; then
    export batch_field=$(sanitise_field "$batchField")
    echo "Batch field: $batch_field"
fi

# This script is under /bin of the scxa-workflows repo
run_tertiary_workflow.sh

pushd $SCXA_OUTDIR > /dev/null 

if [ $? -eq 0 ]; then
                    
    # Group associated matrix files
    for matrix_type in raw_filtered filtered_normalised; do
                        
        pushd matrices > /dev/null 
        zip -r ${matrix_type}.zip ${matrix_type}
        popd > /dev/null 
    done
                   
    set +e
                    
    marker_files=$(ls markers_* 2>/dev/null | grep -v markers_resolution)
    if [ $? -ne 0 ]; then
         echo "No marker files present"
         touch markers/NOMARKERS
     else
        mv $marker_files markers
    fi

    popd > /dev/null
fi 
