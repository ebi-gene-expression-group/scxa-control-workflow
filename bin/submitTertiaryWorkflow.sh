#!/usr/bin/env bash

# Submit a workflow for Galaxy

expName=$1
species=$2
countMatrix=$3
referenceGtf=$4
cellMetadata=$5
isDroplet=$6
galaxyCredentials=$7

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
export GALAXY_INSTANCE=ebi_cluster_ge_user
       
if [ "$isDroplet" = 'True' ]; then
    export FLAVOUR=w_droplet_clustering
else
    export FLAVOUR=w_smart-seq_clustering
fi

run_flavour_workflows.sh
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
                    
    marker_files=$(ls markers_* 2>/dev/null | grep -v markers_clusters_resolution)
    if [ $? -ne 0 ]; then
        echo "No marker files present"
        touch markers/NOMARKERS
    else
        mv $marker_files markers
    fi
fi 
