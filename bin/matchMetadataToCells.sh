#!/usr/bin/env bash

cellMeta=$1
countMatrix=$2
outFile=$3

# Borrowed this incantation from here:
# https://unix.stackexchange.com/questions/349363/compare-two-files-and-print-matches-large-files.
# We use 'join' to grab cell names from the metadata file, and to preserve
# ordering we add an index which we-re-sort by at the end

function match_barcodes_to_meta {
    <$1 nl -s $'\t' |
    sort -t$'\t' -k 2 |
    join -t$'\t' -1 2 - <(sort -k 1,1 $2 ) |
    sort -t$'\t' -k 2,2n |
    cut -d$'\t' -f2 --complement
}

# Make a santised header that won't be an issue for downstream tools

head -n 1 $cellMeta | \
    sed -e "s/\t\(Characteristic\|Factor\)/\t/g" | \
    sed -e "s/\t\[/\t/g" | \
    sed -e "s/\]\t/\t/g" | \
    sed -e "s/\] ontology/ ontology/g" |
    sed "s/ /_/g" > header.tsv
    
# Find the matching lines from the meatadata in the correct order wrt the
# matrix

zipdir=$(unzip -qql ${countMatrix} | head -n1 | tr -s ' ' | cut -d' ' -f5- | sed 's|/||')
unzip -p $countMatrix ${zipdir}/barcodes.tsv > barcodes.tsv

# Run the match

match_barcodes_to_meta barcodes.tsv $cellMeta > ${outFile}.tmp

# Just check that we've actually got the match

cat ${outFile}.tmp | awk '{print $1}' > test_barcodes.tsv

diff barcodes.tsv test_barcodes.tsv > /dev/null 2>&1
if [ $? -ne 0 ]; then
    echo "Unable to match all barcodes in cell metadata, trying a fix to remove odd characters from cell meta" 1>&2
    
    # Try again with a fixed version of the metadata- the SMART pipeline
    # removes odd characters from techrep group names
    
    cat $cellMeta | awk -F $'\t' 'BEGIN {OFS = FS} {gsub(/[^a-zA-Z0-9]/,"_",$1); print}' > ${cellMeta}.fixed

    match_barcodes_to_meta barcodes.tsv ${cellMeta}.fixed > ${outFile}.tmp
    cat ${outFile}.tmp | awk '{print $1}' > test_barcodes.tsv
    diff barcodes.tsv test_barcodes.tsv > /dev/null 2>&1
    
    if [ $? -ne 0 ]; then
        echo "Nope, still failed" 1>&2
        exit 1
    else
        echo "Fix succeeded!"
    fi
fi 

echo "All cells matched successfully!"    

cat header.tsv ${outFile}.tmp > ${outFile}
