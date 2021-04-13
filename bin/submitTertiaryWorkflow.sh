#!/usr/bin/env bash

scriptDir=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $scriptDir/utils.sh

# Submit a workflow for Galaxy

expName=$1
species=$2
confFile=$3
countMatrix=$4
referenceGtf=$5
cellMetadata=$6
isDroplet=$7
galaxyCredentials=$8
galaxyInstance=${9:-ebi_cluster_ge_user}

rm -rf $SCXA_RESULTS/$expName/$species/bundle

# Galaxy workflow wants the matrix components separately

zipdir=$(unzip -qql ${countMatrix} | head -n1 | tr -s ' ' | cut -d' ' -f5- | sed 's|/||')
unzip ${countMatrix} 
       
gzip ${zipdir}/matrix.mtx
gzip ${zipdir}/genes.tsv
gzip ${zipdir}/barcodes.tsv

export species=$species
export expName=$expName
export gtf_file=$referenceGtf
export matrix_file=${zipdir}/matrix.mtx.gz
export genes_file=${zipdir}/genes.tsv.gz
export barcodes_file=${zipdir}/barcodes.tsv.gz
export cell_meta_file=$cellMetadata
export tpm_filtering='False'
export create_conda_env=no
export GALAXY_CRED_FILE=$galaxyCredentials
export GALAXY_INSTANCE=$galaxyInstance
       
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
if [ "$batchField" != 'None' ]; then
    export batch_field=$(sanitise_field "$batchField")
    echo "Batch field: $batch_field"
fi

# This script is under /bin of the scxa-workflows repo

export state_file=$TMPDIR/${expName}.${species}.galaxystate
run_tertiary_workflow.sh

if [ $? -eq 0 ]; then
    mkdir -p matrices
                    
    # Group associated matrix files
    for matrix_type in raw_filtered filtered_normalised; do
        mkdir -p matrices/${matrix_type} 
                        
        for file in matrix.mtx genes.tsv barcodes.tsv; do
            if [ ! -e ${matrix_type}_${file} ]; then
                echo "${matrix_type}_${file} does not exist" 1>&2
                exit 2
            else
                mv ${matrix_type}_${file} matrices/${matrix_type}/${file}
            fi
        done
                        
        pushd matrices > /dev/null 
        zip -r ${matrix_type}.zip ${matrix_type}
        popd > /dev/null 
    done
                   
    # Organise other outputs
                 
    mkdir -p tsne && mv tsne_* tsne 
    mkdir -p umap && mv umap_* umap
    mkdir -p markers
    set +e
                    
    marker_files=$(ls markers_* 2>/dev/null | grep -v markers_resolution)
    if [ $? -ne 0 ]; then
        echo "No marker files present"
        touch markers/NOMARKERS
    else
        mv $marker_files markers
    fi

    rm -f state_file
fi 
