#!/usr/bin/env bash

species=$1
field=$2

filename=$(grep $field "$ISL_CONFIG_DIR/$species.conf" | awk -F'=' '{print $2}')
if [ "$filename" == '' ];then
    echo -n 'None'
else
    echo -n "$filename"
fi
