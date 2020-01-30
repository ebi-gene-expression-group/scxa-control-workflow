#!/usr/bin/env bash

scriptDir=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $scriptDir/utils.sh

sdrfFile=$1
outDir=$2

if [ ! -e "$sdrfFile" ]; then
    die "$sdrfFile does not exist"    
fi 

# Extract the field containing organism information

organismField=$(sed -n $'1s/\\\t/\\\n/gp' $sdrfFile | grep -nxP 'Characteristics ?\[organism\]' | cut -d: -f1)
if [ $? -ne 0 ]; then
    die "No SDRF field available containing organism information"
fi

mkdir -p $outDir

# Now split the file into parts by species

expName=$(echo $(basename $sdrfFile) | awk -F'.' '{print $1}')
headerLine=$(head -n 1 $sdrfFile)

tail -n +2 "$sdrfFile" | awk -v ind=$organismField -F'\t' '{gsub(/ /, "_", $ind); print tolower($ind)}' | sort | uniq | while read -r organism; do
    echo -e "$headerLine" > $outDir/${expName}.${organism}.sdrf.txt
done

awk -v ind=$organismField -v dir="$outDir" -v expName=$expName -F"\t" 'BEGIN { OFS = FS } {print >> dir"/"expName "." tolower(gensub(/ /, "_", "g", $ind)) ".sdrf.txt"}' <(tail -n +2 $sdrfFile)
