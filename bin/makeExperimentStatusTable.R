#!/usr/bin/env Rscript

library(DT)

cl <- commandArgs(trailingOnly = TRUE)

outfile <- cl[1]

workflow_root=Sys.getenv('SCXA_WORKFLOW_ROOT')

idf_files=Sys.glob(file.path(workflow_root, 'metadata/*/*/*.idf.txt'))
exclusions <- read.delim(file.path(workflow_root, 'results', 'excluded.txt'), header = F)

status_table <- do.call(rbind, lapply(idf_files, function(idf){

    sdrf <- sub('\\.idf', '.sdrf', idf)
    cells <- sub('\\.idf', '.cells', idf)

    metadata_files <- c(idf, sdrf)
    if (file.exists(cells)){
        metadata_files <- c(metadata_files, cells)
    }

    metamods <- unlist(lapply(metadata_files, function(x) as.character(file.info(x)$mtime)))
    metadata_order <- order(as.Date(metamods), decreasing = TRUE)
    metadata_last_modified <- metamods[metadata_order][1]
    last_modified_metadata_file <- metadata_files[metadata_order[1]]

    exp_id=unlist(strsplit(basename(idf), '\\.'))[1] 
    exp_dir=file.path(workflow_root, 'results', exp_id)

    # The fields to be reported for each experiment

    species <- 'unknown'
    protocols <- 'unknown'
    excluded <- 'no'
    up_to_date <- 'no'

    file_tests <- list(
        Configured = 'no',
        Quantified = 'no',
        Aggregated = 'no',
        'Tertiary analysis complete' = 'no',
        Bundled = 'no'
    )
    
    # Check for exclusion status

    excluded <- ifelse(exp_id %in% exclusions$V1, 'yes', 'no')
    excluded_reason <- 'N/A'    
    
    if (excluded == 'no'){
        # Check for species
        
        if ( dir.exists(exp_dir)){
            species <- list.files(exp_dir) 
        }
    }else{
        excluded_reason <- exclusions$V2[exclusions$V1 == exp_id]
    }

    # If we have the species dir, then something has happened, proceed to check
    # for config step (which will separate things by protocol)

    analysis_completed <- 'N/A'

    if ( species[1] != 'unknown' ){
        for (spec in species){
            conf_glob <- file.path(workflow_root, 'conf/study', paste(exp_id, spec, '*', 'conf', sep = '.'))
            conf_files <- Sys.glob(conf_glob)
            if (length(conf_files) > 0){
                protocols <- paste(unlist(lapply(conf_files, function(cf) unlist(strsplit(basename(cf), '\\.'))[3])), collapse=',')
            }

            test_files <- list(
                Configured = conf_glob,
                Quantified = file.path(workflow_root, 'results', exp_id, spec, 'quantification/*/quantification.log'),
                Aggregated = file.path(workflow_root, 'results', exp_id, spec, 'aggregation/matrices/counts_mtx.zip'),
                'Tertiary analysis complete' = file.path(workflow_root, 'results', exp_id, spec, 'scanpy/clustering_software_versions.txt'),
                Bundled = file.path(workflow_root, 'results', exp_id, spec, 'bundle/MANIFEST')
            )

            # Generate statuses from file presence/ absence

            file_tests <- lapply(test_files, function(x){
              ifelse(length(Sys.glob(x) > 0), 'yes', 'no')  
            })

            if ( file.exists(test_files$Bundled)){
                analysis_completed <- as.character(file.info(test_files$Bundled)$mtime)
            }

            # Check updates

            if(file_test('-nt', test_files$Bundled, last_modified_metadata_file)){
                up_to_date <- 'yes'
            }
        } 
    }

    cbind(
        file_tests,
        data.frame(
            Experiment=exp_id,
            'Metadata last modified' = metadata_last_modified,
            Species=species,
            Protocols=protocols,
            Excluded=excluded,
            'Excluded reason' = excluded_reason,
            'Up to date' = up_to_date,
            'Analysis completed' = analysis_completed,
            check.names = FALSE
        )
    )[, c('Experiment', 'Metadata last modified', 'Species', 'Protocols', 'Excluded', 'Excluded reason', 'Configured', 'Quantified', 'Aggregated',  'Tertiary analysis complete', 'Bundled', 'Analysis completed', 'Up to date')] 
}))

write.table(status_table, file = 'status_table.txt')

completed <- which(status_table[['Analysis completed']] != 'N/A')
completed <- completed[order(status_table[completed, 'Analysis completed'], decreasing = TRUE)]
not_completed <- which(status_table[['Analysis completed']] == 'N/A' & status_table[['Excluded']] == 'no')
not_completed <- not_completed[order(status_table[not_completed, 'Metadata last modified'], decreasing = TRUE)]
excluded <- which(status_table[['Excluded']] == 'yes')

htmlwidgets::saveWidget(
    datatable(status_table[c(not_completed, completed, excluded),], extensions = c('Buttons', 'FixedHeader'), rownames=FALSE, filter = 'top', options=list(dom = 'rtiBp', pageLength="500", fixedHeader = TRUE)) 
        %>% formatStyle(c('Quantified', 'Configured', 'Aggregated', 'Bundled', 'Tertiary analysis complete', 'Up to date'), backgroundColor = styleEqual(c('no', 'yes'), c('red', 'lightgreen'))) 
        %>% formatStyle(c('Excluded'), backgroundColor = styleEqual(c('yes', 'no'), c('red', 'lightgreen')))
        %>% htmlwidgets::prependContent(htmltools::tags$h1("SCXA analysis pipeline status"))
        %>% htmlwidgets::prependContent(htmltools::p(paste0('Last updated: ', date()))),
    outfile,
    selfcontained = FALSE,
    title = "SCXA analysis pipeline status"
)

