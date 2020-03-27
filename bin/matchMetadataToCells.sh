#!/usr/bin/env bash

cellMeta=$1
countMatrix=$2
outFile=$3

# Make a santised header that won't be an issue for downstream tools

head -n 1 $cellMeta | \
    sed -e "s/\t\(Characteristic\|Factor\)/\t/g" | \
    sed -e "s/\t\[/\t/g" | \
    sed -e "s/\]\t/\t/g" | \
    sed -e "s/\] ontology/ ontology/g" |
    sed "s/ /_/g" > ${outFile}.tmp
    
# Find the matching lines from the meatadata in the correct order wrt the
# matrix

zipdir=$(unzip -qql ${countMatrix} | head -n1 | tr -s ' ' | cut -d' ' -f5- | sed 's|/||')
unzip -p $countMatrix ${zipdir}/barcodes.tsv | while read -r l; do 
    grep -P "^$l\t" $cellMeta
     if [ $? -ne 0 ]; then 
        echo "Missing metadata for $l" 1>&2
        exit 1
     fi 
done >> ${outFile}.tmp

mv ${outFile}.tmp ${outFile}
