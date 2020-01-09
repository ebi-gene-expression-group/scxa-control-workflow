#!/usr/bin/sh

expName=$1
isDroplet=$2
barcodesFile=$3
sampleMetadataFile=${4}
dropletMetadataFile=${5:-''}

if [ "$isDroplet" = 'False' ]; then
    
    # For SMART experiments, grep out the metadata lines we need for each run.
    # If matching metadata can't be found then there's something wrong

    cellMetaFile=$metaFile

    cat $barcodesFile | while read -r l; do 
        grep -P "^$l\t" $sampleMetadataFile
         if [ $? -ne 0 ]; then 
            echo "Missing metadata for $l" 1>&2
            exit 1
         fi 
    done > cell_metadata.tsv

else

    # For droplet experiments life is complicated. Cell IDs are of the form
    # SAMPLE-BARCODE, so we split the sample ID out and use that to to
    # duplicate sample info per row. Where a .cells.tsv file is available, this
    # can be used to add cell-wise metadata.

    # 1. Duplicate library-wise metadata across component cells, looking for
    # the sample identifier of each cell ID in the first column of the sample
    # metadata

    echo "Copying library-wise metadata across cells..."

    head -n 1 $sampleMetadataFile > sample_metadata.tsv && cat $barcodesFile | awk -F'-' '{print $1}' | while read -r l; do 
        grep -P "^$l\t" $sampleMetadataFile
         if [ $? -ne 0 ]; then 
            echo "Missing metadata for \$l" 1>&2
            exit 1
         fi 
    done >> sample_metadata.tsv

    # 2. Derive a cell list that will be the first column of the final output,
    # and paste in front of the sample metadata

    echo "cell" > cells.tsv
    cat $barcodesFile >> cells.tsv

    if [ "$(cat cells.tsv | wc -l)" = "$(cat sample_metadata.tsv | wc -l)" ]; then 
        paste -d "\t" cells.tsv sample_metadata.tsv > cell_metadata.tsv.tmp
    else
        echo "Number of cells not equal to number of lines in cell-expanded sample metadata" 1>&2
        exit 1
    fi

    # 3. Now take cell-wise metadata found in the *.cells.tsv file (where
    # present), find the matching lines for each cell identifier, and add the
    # resulting columns to the output

    type=$(echo $expName | awk -F'-' '{print $2}')
    cells_file_name="$SCXA_WORKFLOW_ROOT/metadata/$type/${expName}/${expName}.cells.txt"
    
    if [ -e "$cells_file_name" ]; then

        echo "Found cell metadata at $cells_file_name, adding to output metadata table..."
    
        # Barcodes without entries in the annotation, print the right number of
        # delimiters such that we get empty fields

        emptystring=$(head -n 1 $cells_file_name | sed s/[^\\t]//g)

        head -n 1 $cells_file_name > droplet_cell_metadata.tsv && cat $barcodesFile | while read -r l; do 
            grep -P "^$l\t" $cells_file_name
            if [ $? -ne 0 ]; then  
                echo -e "$emtpyString"
            fi
        done >> droplet_cell_metadata.tsv

        if [ "$(cat cell_metadata.tsv.tmp | wc -l)" = "$(cat droplet_cell_metadata.tsv | wc -l)" ]; then 
            paste -d "\t" cell_metadata.tsv.tmp droplet_cell_metadata.tsv > cell_metadata.tsv    
        else
            echo "Inconsistent number of lines derived from cell metadata" 1>&2
            exit 1
        fi

    else
        "No cells file present at $cells_file_name"
        cp cell_metadata.tsv.tmp cell_metadata.tsv
    fi
fi






