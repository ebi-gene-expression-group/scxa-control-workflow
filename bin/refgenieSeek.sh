#!/usr/bin/env bash

set -e

species=$1
referenceType=$2
spikes=$3
outDir=${4:-'reference'}

export COLUMNS=500

function find_orig_refgenie_asset_name() {
    asset_path=$1
    basename $(grep "cp " $(dirname $asset_path)/_refgenie_build/refgenie_commands.sh | head -n 1 | awk '{print $2}')
}

# Get the digest for the currently used alias. We want to resolve this to a
# more specific-looking name so files are more descriptive. e.g. 
# homo_sapiens--latest to homo_sapiens--GRCh38.

spikesPart=''
if [ -n "$spikes" ]; then
    spikesPart="--spikes_$spikes"
fi

digest=$(refgenie alias get -a ${species}--${referenceType}${spikesPart}) 
if [ $? -ne 0 ] || [ -z "$digest" ]; then
    echo "Refgenie doesn't seem to have the alias ${species}--${referenceType}" 1>&2
    exit 1
else
    echo "Got digest $digest for ${species}--${referenceType}${spikesPart}" 1>&2
    reference=$(refgenie alias get | grep $digest | awk -F ' │ ' '{print $2}' | awk -F',' '{print $1}')
fi
refinfo=$( refgenie list -g $reference )

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

kallisto_version=$(grep "kallisto=" ${SCXA_WORKFLOW_ROOT}/workflow/scxa-workflows/w_smart-seq_quantification/envs/kallisto.yml | awk -F'=' '{print $2}' | tr -d '\n')
kallisto_index=$(refgenie seek $reference/kallisto_index:cdna_${referenceType}--kallisto_${kallisto_version})

salmon_version=$(grep "salmon=" ${SCXA_WORKFLOW_ROOT}/workflow/scxa-workflows/w_droplet_quantification/envs/alevin.yml | awk -F'=' '{print $2}' | tr -d '\n')
salmon_index=$(refgenie seek $reference/salmon_index:cdna_${referenceType}--salmon_${salmon_version})

mkdir -p $outDir
ln -s $fasta $outDir/$(find_orig_refgenie_asset_name $fasta)
ln -s $gtf $outDir/$(find_orig_refgenie_asset_name $gtf)
ln -s $(dirname $salmon_index) $outDir/salmon_index
ln -s $(ls $(dirname $kallisto_index)/*.idx) $outDir/$(find_orig_refgenie_asset_name $fasta).index
