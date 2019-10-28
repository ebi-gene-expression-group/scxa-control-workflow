#!/usr/bin/env bash

conf_file=$1
field=$2

if [ ! -e "$conf_file" ]; then
    echo "Supplied iRAP conf file $conf_file does not exist" 1>&2
    exit 1
fi

filename=$(grep $field "$conf_file" | awk -F'=' '{print $2}')
if [ "$filename" == '' ];then
    echo -n 'None'
else
    echo -n "$filename"
fi
