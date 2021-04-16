#!/usr/bin/env python

import argparse
import yaml
import sys
from functools import reduce

parser = argparse.ArgumentParser(description='Extract parameter values from a Nextflow configuration')

parser.add_argument('--paramFile', dest='param_file', help='a Nextflow configuration file')
parser.add_argument('--paramKeys', dest='param_keys', help='a Nextflow configuration file')

args = parser.parse_args()
nf_content = ''

# Read and make YAML-like for parsing

with open(args.param_file, 'r') as nf_file:
    for line in nf_file:
        if (not 'includeConfig' in line) and (not '}' in line) and (not line.startswith('//')):
            line = line.replace('{', ':')
            line = line.replace('=', ':')
            line = line.replace('//', '#')
            nf_content += line

# Parse content as YAML

nf_config = yaml.load(nf_content, Loader=yaml.FullLoader)

# Get the value out of the dict by list keys

result = reduce(dict.get, args.param_keys.split(','), nf_config)

if result is not None:
    sys.stdout.write(str(result))
else:
    sys.stdout.write('None')
