#!/usr/bin/env bash

species=$1
referenceType=$2
outDir=${3:-'reference'}

export columns=500

function find_orig_refgenie_asset_name() {
    asset_path=$1
    basename $(grep "cp " $(dirname $asset_path)/_refgenie_build/refgenie_commands.sh | head -n 1 | awk '{print $2}')
}

function find_index_path(){
    local reference=$1
    local index_type=$2
    local reference_type=$3

    index_name=$(echo -e "$refinfo" | grep ${index_type}_index | awk -F' │ ' '{print $3}' | sed 's/,//g' | tr ' ' '\n' | grep cdna_${reference_type})
    if [ $? -ne 0 ]; then
        echo "No $index_type index available for $reference" 1>&2
    else
        asset_path=$(dirname $(refgenie seek $reference/${index_type}_index:$index_name))

        if [ "$index_type" == 'kallisto' ]; then
            asset_path=$(echo -n $(ls $asset_path/*.idx))
        fi
        echo -n "$asset_path"
    fi
}

# Get the digest for the currently used alias. We want to resolve this to a
# more specific-looking name so files are more descriptive. e.g. 
# homo_sapiens--latest to homo_sapiens--GRCh38.

digest=$(refgenie alias get -a ${species}--${referenceType}) 
if [ $? -ne 0 ]; then
    echo "Refgenie doesn't seem to have the alias ${species}--${referenceType}" 1>&2
    exit 1
else
    reference=$(refgenie alias get | grep $digest | awk -F ' │ ' '{print $2}' | awk -F',' '{print $1}')
fi

refinfo=$( env COLUMNS=500 refgenie list -g $refgenie )

# Use references from the ISL setup. Use pre-baked conversions for when
# ensembl species paths don't match organism

species=$species
if [ -n "$NONSTANDARD_SPECIES_NAMES" ] && [ -e "$NONSTANDARD_SPECIES_NAMES" ]; then
    set +e
    species_line=$(grep "^$species$(printf '\t')" $NONSTANDARD_SPECIES_NAMES)
    if [ $? -eq 0 ]; then
        species=$(echo -e "$species_line" | awk '{print $2}')
    fi
    set -e
fi

fasta=$(refgenie seek $reference/fasta_txome:cdna_$referenceType)
gtf=$(refgenie seek $reference/ensembl_gtf:$referenceType )
contamination=$(refgenie seek contamination/bowtie2_index)
kallisto_index=$(find_index_path $species kallisto $referenceType)
salmon_index=$(find_index_path $species salmon $referenceType)

mkdir -p $outDir
ln -s $fasta $outDir/$(find_orig_refgenie_asset_name $fasta)
ln -s $gtf $outDir/$(find_orig_refgenie_asset_name $gtf)
echo -en "$contamination" > $outDir/contamination.txt
ln -s $salmon_index $outDir/salmon_index
ln -s $kallisto_index $outDir/$(find_orig_refgenie_asset_name $fasta).index
