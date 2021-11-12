#!/usr/bin/env Rscript

# Script to make sure genen annotation fields necessary for analysis are all
# present, even if we have to spoof them

cl <- commandArgs(trailingOnly=TRUE)

infile <- cl[1]
outfile <- cl[2]

anno <- read.delim(infile, stringsAsFactors=FALSE)

if (! 'gene_name' %in% colnames(anno)){
    anno$gene_name <- anno$gene_id
}

write.table(anno, file = outfile, sep = "\t", quote=FALSE, row.names = FALSE)

