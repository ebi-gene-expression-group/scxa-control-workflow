#!/usr/bin/env bash

scriptDir=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $scriptDir/utils.sh

mtxArchive=$1
confFile=$2
outFile=${3:-'cell_to_library.txt'}

for inFile in "$mtxArchive" "$confFile"; do
    if [ -z "$inFile" ]; then
        die "Empty input file supplied" 1 
    elif [ ! -e "$inFile" ]; then
        die "Input file $inFile does not exist" 2
    fi
done

# Uncompress MTX to extract matrix, genes, barcodes

zipdir=$(unzip -qql ${mtxArchive} | head -n1 | tr -s ' ' | cut -d' ' -f5- | sed 's|/||')

if [ $? -ne 0 ]; then
    die "Unzip of $mtxArchive failed" 3
fi

sampleField=$(parseNfConfig.py --paramFile $confFile --paramKeys params,fields,run)
techrepField=$(parseNfConfig.py --paramFile $confFile --paramKeys params,fields,techrep)

if [ $techrepField != 'None' ]; then
    sampleField=$techrepField
fi

echo "# $sampleField" > $outFile
unzip -p $mtxArchive ${zipdir}/barcodes.tsv | while read -r b; do 
    barcode=${b##*-} 
    run=${b/-$barcode/''}
    echo -e "$b\t$run" 
done >> $outFile
