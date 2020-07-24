#!/usr/bin/env bash

sdrf_file=$1
cells_file=$2
cell_type_fields=$3

# Determine if there is a valid cell types file

has_cells=0
if [ -n "$cells_file" ] && [ -e "$cells_file" ];then
    cells_filesize=$(stat --printf="%s" $(readlink ${cells_file}))
    if [ $cells_filesize -gt 0 ]; then
        has_cells=1
    fi
fi

# Loop over the cell types fields, checking in the SDRF and cells files

IFS=',' read -ra ctfs <<< "$cell_type_fields"
for ctf in "${ctfs[@]}"; do
    
    echo "Checking for $ctf" 1>&2

    # Allow for bracketed field in cells file, preceded by any run of non-tab
    # characters (Factor value etc plus spaces)
    
    sdrf_cell_type_field=$(head -n 1 $sdrf_file | grep -io "[^$(echo -e "\t")]*$(echo "\\[$ctf\\]") *")    
    if [ $? -eq 0 ]; then
        echo -n "$ctf"
        break
    elif [ $has_cells -eq 1 ]; then

        # Allow for bracketed field in cells file, preceded by any run of
        # non-tab characters (Factor value etc plus spaces)

        echo "No cell type in SDRF, checking .cells file" 1>&2
        cells_cell_type_field=$(head -n 1 $cells_file | grep -io "[^$(echo -e "\t")]*$(echo "\\[$ctf\\]") *")    
        if [ $? -eq 0 ]; then
            echo -n "$ctf"
            break
        else
            
            # Allow for bare field in cells file, accounting for leading or
            # trailing spaces
        
            cells_cell_type_field=$(head -n 1 $cells_file | grep -io "[$(echo -e "\t")] *$(echo "$ctf") *")    
            if [ $? -eq 0 ]; then
                echo -n "$ctf"
                break
            fi
        fi
    fi
done
