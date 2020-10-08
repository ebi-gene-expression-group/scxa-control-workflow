#!/usr/bin/env bash

cellMeta=$1
confFile=$2
outFile=$3

# Extract things we need from the conf file

cellTypeField=$(parseNfConfig.py --paramFile $confFile --paramKeys params,fields,cell_type)

if [ "$cellTypeField" != 'None' ]; then
    cp ${cellMeta} ${outFile}.tmp

    echo "Cell type field found, searching for singlets to remove"
    head -n 1 ${cellMeta} | awk -F"\t" -v fieldName=$cellTypeField '{for(i=1;i<=NF;i++) if($i==fieldName) print i}' | head -n 1 | while read -r fieldNum; do
        tail -n +2 ${cellMeta} | awk -F'\t' -v fieldNum=$fieldNum '{print $fieldNum}' | sort | uniq -c | sed -e 's/^[ ]*//' | while read -r l; do 
            wordFreq=$(echo -e "$l" | awk '{print $1}'); 
            if [ $wordFreq -lt 2 ]; then
                singleValue=$(echo -e "$l" | cut -d" " -f2- | sed 's/\//\\\//g' | sed 's/(/\\\(/g' | sed 's/)/\\\)/g') 
                echo -e "$singleValue is a single value, would take it out"
                awk -v fieldNum=$fieldNum 'BEGIN{FS=OFS="\t"} {gsub(/'"^${singleValue}$"'/, "", $fieldNum)} 1' ${outFile}.tmp > ${outFile}.tmp.tmp && mv ${outFile}.tmp.tmp ${outFile}.tmp      
            fi
        done
    done
    
    mv ${outFile}.tmp $outFile
else
    cp -P $cellMeta $outFile
fi
    

