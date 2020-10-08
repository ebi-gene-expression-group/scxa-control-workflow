#!/usr/bin/env Rscript

source(file.path(Sys.getenv(c("SCXA_BIN")), "utils.R"))

# Load optparse we need to check inputs

suppressPackageStartupMessages(require(optparse))

# Load common functions

suppressPackageStartupMessages(require(workflowscriptscommon))

################################################################################
# Argument parsing
################################################################################

usage <- "sdrfToNfConf.R --name exp_name --species species --sdrf file --out_conf  [options]"
filenames <- c("sdrf_file","idf_file")

option_list <- list(
  make_option(c("-n", "--name"), type="character", dest="name", default=NULL, help="(short) name of the experiment (no spaces)"),
  make_option(c("--species"), type="character", dest="species", default=NULL, help="Species name."),
  make_option(c("--nc"), type="numeric", dest="nc", default=0, help="estimated number of clusters [default %default]"),
  make_option(c("--nc_window"), type="numeric", dest="nc_window", default=10, help="Clusters window [default %default]"),  
  make_option(c("-s", "--sdrf"), type="character", dest="sdrf_file", default=NA, help="SDRF file name"),
  make_option(c("-i", "--idf"), type="character", dest="idf_file", default=NULL, help="IDF file name"),
  make_option(c("-e", "--cells"), type="character", dest="cells_file", default=NULL, help="Cells file name (droplet experiments only)"),
  make_option(c("-f", "--cell_meta_fields"), type="character", dest="cell_meta_fields", default=NULL, help="Comma-separated list of cell metadata fields important for analysis"),
  make_option(c("-t", "--cell_type_fields"), type="character", dest="cell_type_fields", default=NULL, help="Comma-separated list of possible fields for cell type"),
  make_option(c("--sc"), action="store_true",default=FALSE,dest="is_sc",help="Single cell experiment?. Default is FALSE."),
  make_option(c("-c","--check_only"),action="store_true",dest="check_only",default=FALSE,help="Checks the sdrf/idf but skips the generation of the configuration file."),
  make_option(c("-v","--verbose"),action="store_true",dest="verbose",default=FALSE,help="Produce verbose output."),
  make_option(c("--debug"),action="store_true",dest="debug",default=FALSE,help="Debug mode"),
  make_option(c("-o", "--out_conf"), type="character",default=NULL,help="Output directory for configuration files.")
)

multiple.options = list()

mandatory <- c("sdrf_file","name", "idf_file")

opt <- wsc_parse_args(option_list, mandatory = mandatory)

################################################################################
# Function definitions
################################################################################

pinfo <- function(...) {
    cat(paste0("[INFO ",format(Sys.time(), format="%d/%m-%H:%M"),"] ",...,"\n"))
}

# Record warnings as we run checks

warnings.lst <- c()
addWarning <- function(...) {
  warnings.lst <- append(warnings.lst,paste0(...))
  assign("warnings.lst",warnings.lst,envir = .GlobalEnv)
}

printAllWarnings <- function() {
  n <- length(warnings.lst)
  if ( n > 0 ) return
  pinfo(n," warnings. Unique entries:")
  prints <- sapply(unique(warnings.lst),pwarning)
}

# Print info only in the verbose case

print.info <- function(...){
  if (opt$verbose){
    pinfo(paste(...))
  }
}

# Simplify column names

normalize.cols <- function(vals,strip.brackets=TRUE) {
  ## try to workaround typos & inconsistencies
  x <- vals
  if (strip.brackets)
    x <- gsub("\\].*","",gsub(".*\\[","",x))
  else {
    x <- gsub(" +\\]","]",gsub(" +\\[","[",x))
    x <- gsub("] +","]",x)
  }
  # a space or underscore seems to be the same thing
  x <- gsub("_+"," ",x) ## no underscores
  return(tolower(gsub(" +"," ",x))) ##lowercase and single space
}

# Return the actual col name for query name

getActualColnames <- function(cols, sdrf){
   
    sdrf.cols <- colnames(sdrf)

    cols.norm <- normalize.cols(cols)
    sdrf.cols.norm <- normalize.cols(sdrf.cols)

    cols.norm.comments <- paste0('comment[', cols.norm, ']')
    cols.norm.factors <- paste0('factorvalue[', cols.norm, ']')
    cols.norm.characteristics <- paste0('characteristics[', cols.norm, ']')

    actual.colnames <- unlist(lapply(1:length(cols), function(i){
        if (cols.norm[i] %in% sdrf.cols.norm){
            sdrf.cols[match(cols.norm[i], sdrf.cols.norm)]
        }else if (cols.norm.factors[i] %in% sdrf.cols.norm){
            sdrf.cols[match(cols.norm.factors[i], sdrf.cols.norm)]
        }else if (cols.norm.characteristics[i] %in% sdrf.cols.norm){
            sdrf.cols[match(cols.norm.characteristics[i], sdrf.cols.norm)]
        }else if (cols.norm.comments[i] %in% sdrf.cols.norm){
            sdrf.cols[match(cols.norm.comments[i], sdrf.cols.norm)]
        }else{
            NULL
        }
    }))

    # Deal with possible ambiguities

    if ( length(cols) == 1 && is.null(actual.colnames)){
        if ( cols.norm == 'ena run' ) {
            getActualColnames('run', sdrf)
        }else if (cols.norm == 'library construction'){
            getActualColnames('library preparation', sdrf)
        }else{
            actual.colnames
        }
    }else if (length(cols) != length(actual.colnames)){
        perror(print(paste('Do not have actual cols', paste(actual.colnames, collapse=',') ,'of same length as query', paste(cols, collapse=','))))
        q(status=1)
    }else{
        actual.colnames
    }
}

# Check that a value is filled

populated <- function(y){
  unlist(lapply(y, function(x){
    x != '' & x != -1 & gsub(' ', '', x) != '' & ! is.na(x)
  }))
}

# Check if a method is a droplet protocol

is.droplet.protocol <- function(x){
  is.singlecell && tolower(x) %in% tolower(sc.droplet.protocols)
}

################################################################################
# Basic initalisation 
################################################################################

cat("\n")
pinfo(">>> Commencing SDRF validation")

# Do some introspection on arguments

print.info("### Argument values")

for ( n in names(opt) ) {
  print.info(paste0(n,"=",opt[[n]]))
}

################################################################################
# First round SDRF field checks: before we know about field content 
################################################################################

# Load the SDRF

print.info("Loading SDRF ",opt$sdrf_file," ...")
sdrf <- read.tsv(opt$sdrf_file,header=TRUE,comment.char="",fill=TRUE)
print.info("Loading SDRF...done.")

if ( ! is.null(opt$cells_file)){
    print.info("Loading cells ",opt$cells_file," ...")
    cells <- read.tsv(opt$cells_file,header=TRUE,comment.char="",fill=TRUE)
    print.info("Loading cell-wise metadata...done.")
}else{
    cells <- sdrf
}

# Set some variable names we'll be using a lot. getActualColnames() will set
# things to null when they're not present

run.col <- getActualColnames('Comment [ENA_RUN]', sdrf)
if ( is.null(run.col) ) {
  perror("SDRF: expected RUN or ENA_RUN - none found")
  q(status=1)
}else{
  ## fix the run names (should not include # )
  sdrf[[run.col]] <- gsub("#","_",sdrf[[run.col]])
}
techrep.col <- getActualColnames('technical replicate group', sdrf)
strand.col <- getActualColnames('library strand', sdrf)
library.layout.col <- getActualColnames('library layout', sdrf)
sc.identifier.col <- getActualColnames('single cell identifier', sdrf)
protocol.col <- getActualColnames('library construction', sdrf)
sc.quality.col <- getActualColnames('single cell quality', sdrf)
libary.layout.col <- getActualColnames('library layout', sdrf)
spike.in.col <- getActualColnames('spike_in', sdrf)
fastq.col <- getActualColnames('fastq uri', sdrf)
cell.count.col <- getActualColnames('cell count', sdrf)
hca.bundle.uuid.col <- getActualColnames('HCA bundle uuid', sdrf)
hca.bundle.version.col <- getActualColnames('HCA bundle version', sdrf)
controlled.access.col <- getActualColnames('controlled access', sdrf)
batch.col <- NULL

# We may have been supplied with some cell-wise fields to extract

cell.meta.cols <- NULL
cell.type.col <- NULL

if (! is.null(opt$cell_meta_fields)){

    if ( ! is.null(opt$cells_file)){
        cell.meta.source <- cells
        cell.id.col <- getActualColnames('cell id', cells)
    }else{
        cell.meta.source <- sdrf
        cell.id.col <- run.col
    }
    
    # See if cell type fields specifically have been supplied

    if (! is.null(opt$cell_type_fields)){
        cell_type_fields <- unlist(strsplit(opt$cell_type_fields, ','))
        
        for (ctf in cell_type_fields){
            actf <- getActualColnames(ctf, cell.meta.source)
            if (! is.null(actf)){
                cell.type.col <- actf
                break
            }
        }
    }
    
    cell_meta_fields <- unlist(strsplit(opt$cell_meta_fields, ','))
    cell_meta_fields <- unlist(lapply(cell_meta_fields, function(x) getActualColnames(x, cell.meta.source)))
    cell.meta.cols <- unique(c(cell.id.col, cell.type.col, cell_meta_fields[ ! is.null(cell_meta_fields)]))
}
ena.sample.col <- getActualColnames('ena_sample', sdrf)
organism.col <- getActualColnames('organism', sdrf)

# Define protocols in case of single cell  

sc.opt.cols <- c("single cell quality","input molecule","end bias","single cell library method","read1 file","read2 file","index1 file", "index2 file","index3 file")
supported.single.cell.protocols <- c("smart-seq", "smart-seq2","smarter","smart-like","10xv2","10xv3","drop-seq","seq-well","10x5prime")
sc.droplet.protocols <- c('10xv1', '10xv1a', '10xv1i', '10xv2', '10xv3', 'drop-seq','seq-well', '10x5prime')
  
# Check the protocol and use to determine single-cell

is.singlecell <- FALSE
is.hca <- ! is.null(hca.bundle.uuid.col)

if ( ! is.null(protocol.col)){

    protocol.corrections <- list(
        'smart-like' = 'smart-seq',
        'Smart-seq2' = 'smart-seq2'
    )

    protocols <- tolower(sdrf[[protocol.col]])

    for (prot in names(protocol.corrections)){
        protocols[protocols == prot] <- protocol.corrections[[prot]]
    }

    sdrf[[protocol.col]] <- protocols 

    if (all(protocols %in% supported.single.cell.protocols)){
        print.info("Found single-cell experiment")
        is.singlecell <- TRUE

    }else if (any(protocols %in% supported.single.cell.protocols)){
        perror("Some but not all rows have single-cell protocols")
        q(status=1)
    }else{
        pwarning(paste('None of', paste(unique(protocols), collapse=','), 'indicate single-cell' ))
    }
}else{
    stop(paste("No protocol column in ", paste(colnames(sdrf), collapse=',')))
}

## Factors should not be simultaneously comments, niether should characteristics 

factors <- normalize.cols(colnames(sdrf)[grepl("^factor ?value",colnames(sdrf),ignore.case=TRUE)])
comments <- normalize.cols(colnames(sdrf)[grepl("^comment",colnames(sdrf),ignore.case=TRUE)])
characteristics <- normalize.cols(colnames(sdrf)[grepl("^characteristics",colnames(sdrf),ignore.case=TRUE)])

compare <- list( factors = factors, characteristics = characteristics  )
for (comp in names(compare)){
  common <- intersect( compare[[comp]], comments )
  if (length(common)>0) {
    err_msg <- paste0("SDRF: common ", comp," and comments - ",paste(common,collapse=","))
    perror(err_msg)
    q(status=1)
  }
}

## Verify presence of mandatory columns (SDRF) and regularise names

cols.for.download <- ifelse(is.hca, c(hca.bundle.uuid.col, hca.bundle.version.col) , fastq.col)
expected.cols <- c("Source Name")
expected.comment.cols <- c("LIBRARY_STRATEGY","LIBRARY_SOURCE","LIBRARY_SELECTION","LIBRARY_LAYOUT",cols.for.download)
expected.characteristic.cols <- c()
expected.factor.cols <- c()
opt.cols <- c("ORGANISM","organism part","sex","spike in","molecule","technical replicate group","ENA_RUN","ENA_SAMPLE","Scan Name")

################################################################################
# Load and check the IDF where provided, and cross-reference with the SDRF
################################################################################

idf <- NULL
idf.attrs <- idf.factors <- c()

if ( !is.null(opt$idf_file) ) {
  print.info("Loading IDF...")
  idf <- read.tsv(opt$idf_file,header=FALSE,comment.char="",fill=TRUE,quote="")
  if (is.null(idf) ) {
    perror("Error while loading IDF file ",opt$idf_file)
    q(status=1)
  }
  print.info("Loading IDF...done.")
  
  ## exclude blank lines
  
  idf <- idf[idf$V1!="",,drop=FALSE]
  
  ## SecondaryAccession and SequenceDataURI may be non-unique
  
  allowed_dups <- c('comment[relatedexperiment]', 'comment[additionalfile:txt]', 'comment[sequencedatauri]', 'comment[secondaryaccession]')
  idffieldchecks <- gsub(' ','',tolower(idf[,1])) 

  for (ad in allowed_dups){
    if ( sum(idffieldchecks == ad) > 1 ) {
      addWarning("IDF: Multiple '", ad,"' entries")
    }
  }
  
  idf <- idf[! idffieldchecks %in% allowed_dups, ]
  
  ## basic idf check (non unique rownames). Some duplicate fields allowed.
  
  dups <- unique(idf[duplicated(idf[,1]), 1])
  
  if ( length(dups) > 0 ) {
    dupss <- paste(dups, sep=",", collapse=",")
    perror("Duplicated entries: ", dupss)
    q(status=1)
  }
  
  ## error handling should be improved...
  rownames(idf) <- tolower(gsub("\\s*\\]\\s*","",gsub("^\\s*COMMENT\\s*\\[\\s*","",as.character(idf[,1]),ignore.case=TRUE)))
  idf <- idf[,-1,drop=FALSE]
  
  ## be tolerant about the EA comments. If both EAExperimentType and ExperimentType are present, remove Experimenttype

  if ( 'experimenttype' %in% rownames(idf) && 'eaexperimenttype' %in% rownames(idf)){
    idf <- idf[rownames(idf) != 'experimenttype',]
  }

  rownames(idf) <- gsub("^ea","",rownames(idf))
  
  ## initial check
  if (! "sdrf file" %in% rownames(idf) ) {
    perror("SDRF file missing from IDF file")
    q(status=1)
  }
  
  if ( idf["sdrf file",1]!=basename(opt$sdrf_file) ) {      
    perror("SDRF file in IDF (",idf["sdrf file",1],") differs from given sdrf file name ",basename(opt$sdrf_file))
    q(status=1)
  }
 
  ##expectedclusters
  ##
  ## ExperimentType       differential|baseline (mandatory)
  ## ExpectedClusters     may be empty
  ## AdditionalAttributes - characteristics (optional) exist in the SDRF
  ## AEExperimentType RNA-seq of coding RNA from single cells (single cell)
  ## Protocol Name should match SDRF
  ## Experimental Factor Name (mandatory) exist in the SDRF
  ## note: factors and characteristics should always be lower case!! @Laura
  
  expected.in.idf <- c()
  
  # Required fields for single-cell
  
  #if ( is.singlecell ) {
  #  expected.in.idf <- c("expectedclusters")
  #  names(expected.in.idf) <- c("EAExpectedClusters")
  #}
  
  # Required fields for Atlas
  
  expected.in.idf2 <- c("experimenttype","experimental factor name","protocol name")
  names(expected.in.idf2) <- c("EAExperimentType","Experimental Factor Name","Protocol Name")
  expected.in.idf <- c(expected.in.idf,expected.in.idf2)
  
  # Error where required IDF fields not present
  
  not.present <- (!expected.in.idf %in% rownames(idf))
  if ( sum(not.present)>0 ) {
    perror("IDF incomplete. Missing ",paste(names(expected.in.idf)[not.present],sep=",",collapse=","))
    q(status=1)
  }
  
  # Additional IDF checks for Atlas
  
  if ( ! tolower(idf["experimenttype",1]) %in% c("baseline","differential") ){
    perror("IDF error in EAExperimentType: Invalid value (",idf["experimenttype",1],") expected baseline or differential")
    q(status=1)
  }
   
  ## Check the experiment type
    
  exp.type <- tolower(idf["aeexperimenttype",1])
  valid.exp.types = c( single_cell = "RNA-seq of coding RNA from single cells", bulk = "RNA-seq of coding RNA" )
    
  if ( ! exp.type %in% tolower(valid.exp.types) ){
    perror("IDF error in AEExperimentType: Invalid value (", exp.type,") expected ",paste(valid.exp.types, sep=" or ", collapse=" or "))
    q(status=1)
  }
    
  # Override argument if we find single-cell experiment type

  exp.type.name <- names(valid.exp.types)[tolower(valid.exp.types) == exp.type]

  if (is.singlecell && exp.type.name == 'bulk'){
    addWarning("IDF shows bulk experiment ('RNA-seq of coding RNA'), but SDRF indicated single-cell (so should be 'RNA-seq of coding RNA from single cells')")
  } else if ( exp.type.name == 'single_cell' && ! is.singlecell ){
    addWarning("IDF shows single cell ('RNA-seq of coding RNA from single cells'), but SDRF indicated bulk (so should be 'RNA-seq of coding RNA')")
  } 
 
  ## - Experimental Factor Name (mandatory) exist in the SDRF
  
  idf.factors <- idf["experimental factor name",][idf["experimental factor name",] != '' & ! is.na(idf["experimental factor name",]) ]
  
  ## and match the sdrf. Ignore single cell identifier when validating against SDRF- it's not a real factor 
 
  idf.factors.comp <- sort(setdiff(normalize.cols(idf.factors), c(sc.identifier.col, normalize.cols(sc.identifier.col))))
  sdrf.factors.comp <- sort(setdiff(normalize.cols(factors), c(sc.identifier.col, normalize.cols(sc.identifier.col))))
 
  if ( length(idf.factors.comp) != length(sdrf.factors.comp) || any( idf.factors.comp != sdrf.factors.comp )){
    perror(paste("SDRF/IDF inconsistency: IDF factors (", paste(idf.factors.comp, collapse=','),") do not match (SDRF factors", paste(sdrf.factors.comp, collapse=','), ')'))
    q(status=1)
  }
  
  # Check the ExpectedClusters specification
  
  if ( "expectedclusters" %in% rownames(idf) && idf["expectedclusters",1] != "" ) {
    x <- as.numeric(idf["expectedclusters",1])
    if (is.na(x) || x < 1 ) {
      perror("IDF error in AEExpectedClusters: Invalid value  (",idf["expectedclusters",1],")")
      q(status=1)
    }
    opt$nc <- x
    print.info(paste('Expected cluster number from IDF:', opt$nc))
  }
  
  # AdditionalAttributes - characteristics (optional) exist in the SDRF
  
  if ( "additionalattributes" %in% rownames(idf) && idf["additionalattributes",1] != "" ) {
    idf.attrs <- idf["additionalattributes",]
    
    # exclude "" and ignore case
    idf.attrs <- sapply(idf.attrs[which(idf.attrs != "" && ! is.na(idf.attrs))],tolower)
   
    # Check for fields absent in the SDRF
    not.present <- (! idf.attrs %in% normalize.cols(colnames(sdrf)))
    
    if (sum(not.present)) {
      perror("IDF error in EAAdditionalAttributes:",paste(idf.attrs[not.present],sep=",")," not found in SDRF")
      q(status=1)
    }
  }

  # Retrieve batch field if present

  if ( "batcheffect" %in% rownames(idf)){
    batch.col <- getActualColnames(idf["batcheffect",1], sdrf)
    cell.meta.cols <- c(cell.meta.cols, batch.col)
  }
}

# Having derived possible single-cell status from the IDF we can add
# single-cell specific things to the checks.

if( is.singlecell ) {
  
  expected.cols <- append(expected.cols, c("Material Type","library_construction","single cell isolation"))
  expected.comment.cols <- append(expected.comment.cols, c("LIBRARY_STRAND"))
  opt.cols <- append(opt.cols,sc.opt.cols)

  if ( (! is.null(hca.bundle.uuid.col) ) || ( ! is.null(hca.bundle.version.col) )){
    print.info('Looks like an HCA experiment')
    if ( is.null(hca.bundle.uuid.col)){
      perror("... but no UUID field")
      q(status=1)
    }else if ( is.null(hca.bundle.version.col)){
      perror("... but no version field")
      q(status=1)
    }
  }
}

# Now do the checks

print.info("Checking presence of columns...")

all.expected.fields <- c(expected.cols, expected.characteristic.cols, expected.comment.cols, expected.factor.cols)
all.actual.cols <- getActualColnames( all.expected.fields, sdrf )

# Exit if any of the fields come back null from getActualColnames

if (any(is.null(all.actual.cols))){
    missing.error <- paste0("SDRF: Column(s) ", paste(all.expected.fields[which(is.null(all.actual.cols))], collapse=","), " not found")
    
    perror(missing.error)
    q(status=1)
}
print.info("Checking the presence of columns... done.")

################################################################################
# Further SDRF checks now we know the right fields are there and we can look at
# the content
################################################################################

## Check if there are any controlled access runs

controlled.access='no'
if ( ! is.null(controlled.access.col) && any(sdrf[[controlled.access.col]] == 'yes')){
   controlled.access='yes'
}

## User can specify species, in which case we can discard irrelevant rows

sdrf[[organism.col]] <- gsub(" ","_",tolower(sdrf[[organism.col]]))

if ( !is.null(opt$species) ) {
  if ( !opt$species %in% sdrf[[organism.col]] ) {
    perror("SDRF: Species ", opt$species, " not found")
    q(status=1)
  }
  sdrf <- sdrf[sdrf[[organism.col]] == opt$species, ]
}

print.info("Species: ",paste(unique(sdrf[[organism.col]])),"\n")

## SDRF-wide single-cell checks. This can result in filtering out rows, so do it
## early. In the process we also work out what fastq files are required

if ( is.singlecell ) {
  
  if (! is.null(sc.quality.col)){
    
    sdrf[[sc.quality.col]] <- gsub(" +"," ", tolower(sdrf[[sc.quality.col]]))
    single.cell.quality.counts <- table(sdrf[[sc.quality.col]])
    
    ## Pipeline will discard cells from analysis that are low quality
    ## OK, OK filtered or not OK
    
    print.info(single.cell.quality.counts["not ok"]," runs/cells will be discarded due to 'not ok' quality status")
   
    if ( "not ok" %in% names(single.cell.quality.counts)){ 
      if (single.cell.quality.counts["not ok"] == nrow(sdrf)){
        perror(paste("All rows will be removed by the", sc.quality.col, "criterion - stopping"))
        q(status=1)
      }
    }
  }else{
    print("NO QUALITY COL")
  }
  
  # Check that all sc protocols are valid
  
  protocols <- tolower(protocols)
  
  not.valid.sc.protocol <- unique(protocols[! protocols %in% supported.single.cell.protocols])
  
  if ( length(not.valid.sc.protocol) ) {
    perror("SDRF: Unsupported single cell protocol ", not.valid.sc.protocol)
    q(status=1)
  }
  
  # Ignore single-cell identifier for droplet-only experiments
  
  if (all(protocols %in% sc.droplet.protocols) && (! is.null(sc.identifier.col)) && sc.identifier.col %in% factors) {
    addWarning("You have supplied 'single cell identifier' for an exclusively droplet-based experiment. This column is meaningless in this context and will be ignored")
    factors <- factors[factors != sc.identifier.col]
    sdrf <- sdrf[, colnames(sdrf) != sc.identifier.col]
  }
  
  # Correct layout for droplet libraries
  
  library.layout <- sdrf[[library.layout.col]]
  library.layout[ library.layout == 'PAIRED' & is.droplet.protocol(protocols) ] <- 'SINGLE'
  sdrf[[library.layout.col]] <- library.layout  

  ## For every run, derive required file fields
  
  # Rules for droplet protocols
  # 10Xv1: ?
  # 10Xv1a: read 1 = cDNA, read 2 = UMI barcode, index 1 = cell barcode, index 2 = sample barcode (not required)
  # 10Xv1i: read 1 = interleaved cDNA + UMI barcode, read 2 = interleaved cDNA + UMI barcode, index 1 = cell barcode, index 2 = sample barcode (not required)
  # 10Xv2: read 1 = cell + UMI barcodes, read 2 = cDNA, index 1 = sample barcode (not required)
  # Drop-seq: read 1 = cell + UMI barcodes, read 2 = cDNA
  
  sc.required.fastq.fields = list(
    '10xv1' = c('read1 file', 'read2 file'),    
    '10xv1a' = c('read1 file', 'read2 file', 'index1 file'),         
    '10xv1i' = c('read1 file', 'read2 file', 'index1 file'),
    '10xv2' = c(cols.for.download, 'read1 file', 'read2 file', 'cDNA read', 'umi barcode read', 'cell barcode read'),
    '10xv3' = c(cols.for.download, 'read1 file', 'read2 file', 'cDNA read', 'umi barcode read', 'cell barcode read'),
    'drop-seq' = c(cols.for.download, 'read1 file', 'read2 file', 'cDNA read', 'umi barcode read', 'cell barcode read'), 
    'seq-well' = c(cols.for.download, 'read1 file', 'read2 file', 'cDNA read', 'umi barcode read', 'cell barcode read'), 
    '10x5prime' = c(cols.for.download, 'read1 file', 'read2 file', 'cDNA read', 'umi barcode read', 'cell barcode read'), 
    "smart-seq" = cols.for.download,
    "smart-seq2" = cols.for.download,
    "smarter" = cols.for.download,
    "smart-like" = cols.for.download
  )
  
  row.required.fastq.fields <- lapply(sc.required.fastq.fields[match(protocols, tolower(names(sc.required.fastq.fields)))], getActualColnames, sdrf)
  
  sc.optional.fastq.fields = list(
    '10xv1' = c('index1 file', 'index2 file'),    
    '10xv1a' = c('index2 file'),         
    '10xv1i' = c('index2 file'),
    '10xv2' = c('index1 file'),
    '10xv3' = c('index1 file'),
    'drop-seq' = c(),
    'seq-well' = c(),
    '10x5prime' = c(),
    "smart-seq2" = c(),
    "smarter" = c(),
    "smart-like" = c()
  )

  droplet.protocol.defaults = list(
    '10xv2' = list(
      'umi_barcode_offset' = 16,  
      'umi_barcode_size' = 10,  
      'cell_barcode_size' = 16,  
      'cell_barcode_offset' = 0,  
      'cdna_read_offset' = 0,
      'end' =  '5'
    ),
    '10x5prime' = list(
      'umi_barcode_offset' = 16,  
      'umi_barcode_size' = 10,  
      'cell_barcode_size' = 16,  
      'cell_barcode_offset' = 0,  
      'cdna_read_offset' = 0,
      'end' =  '5'
    ),
    '10xv3' = list(
      'umi_barcode_offset' = 16,  
      'umi_barcode_size' = 12,  
      'cell_barcode_size' = 16,  
      'cell_barcode_offset' = 0,  
      'cdna_read_offset' = 0,
      'end' =  '5'
    ),
    'drop-seq' = list(
      'umi_barcode_offset' = 12,  
      'umi_barcode_size' = 8,  
      'cell_barcode_size' = 12,  
      'cell_barcode_offset' = 0,  
      'cdna_read_offset' = 0,
      'end' =  '5'
    ),
    'seq-well' = list(
      'umi_barcode_offset' = 12,  
      'umi_barcode_size' = 8,  
      'cell_barcode_size' = 12,  
      'cell_barcode_offset' = 0,  
      'cdna_read_offset' = 0,
      'end' =  '5'
    )
  )
  
  # Consider any of the optional file fields that occur in the SDRF
 
  row.optional.fastq.fields <- lapply(sc.optional.fastq.fields[match(protocols, tolower(names(sc.optional.fastq.fields)))], function(x){
    if (is.null(x)){
        NULL
    }else{
        actual.opt.fields <- getActualColnames(x, sdrf)
        actual.opt.fields[!is.null(actual.opt.fields)]    
    }
  })
  
}else{
 
  # For bulk, just check rows for fastq_uri
  
  row.required.fastq.fields <- rep(list(cols.for.download), nrow(sdrf))
  row.optional.fastq.fields <- rep(list(c()), nrow(sdrf))
}

################################################################################
# Check that:
# 1 - all the necessary file fields exist
# 2 - they have defined values for compulsory field for their technology
# 3 - the indicated files actually exist
#
# SDRF files have a FASTQ_URI column for specifying links to FASTQ files, and
# for droplet techs we also have fields like 'read1 file', 'index1 file' etc,
# which are just file names.
#
# The following assumes:
#
#  - FASTQ files can be found in the specified directories
#  - Single-cell droplet protocols have one run per line with multiple files. We ignore the multiple FASTQ URIs
#  - Bulk and SMART-like protocols have one line per file with a single FASTQ URI in each
################################################################################

all.required.fastq.fields <- unique(unlist(row.required.fastq.fields))
fastq.fields.in.sdrf <- getActualColnames(all.required.fastq.fields, sdrf)

## Check the actual fields exist

if (any(is.null(fastq.fields.in.sdrf))){
  perror(paste("FASTQ file field (s)", paste(all.required.fastq.fields[is.null(fastq.fields.in.sdrf)], collapse=','), 'missing'))
  q(status=1)
}

## Determine required file names on each row using the derived field names

fastq.fields <- paste(Reduce(union, row.required.fastq.fields), collapse = ',') 

row.files <- lapply(1:nrow(sdrf), function(row.no){
  req.fields <- row.required.fastq.fields[[row.no]]
  structure(basename(as.character(sdrf[row.no, req.fields, drop = FALSE])), names = req.fields)
})
## Check that the required fields are populated

rows.with.undefined.files <- which(unlist(lapply(row.files, function(x) any(! populated(x)))))

if (length(rows.with.undefined.files) > 0){
  undefined.row.files <- paste(unlist(lapply(rows.with.undefined.files, function(x) paste(paste0('Row ', x, ':'), paste(names(row.files[[x]])[ ! populated(row.files[[x]])], collapse=',')))), collapse="\n")
  
  perror("SDRF: some rows have missing values for required files: \n", undefined.row.files)
  q(status=1)
}

## Add in any optional files with defined values

row.files <- lapply(1:nrow(sdrf), function(row.no){
  opt.fields <- row.optional.fastq.fields[[row.no]]
  if (length(opt.fields) > 0 ){
    opt.files <- basename(sdrf[row.no, opt.fields])
    names(opt.files) <- opt.fields
    c(row.files[[row.no]], opt.files[populated(opt.files)])
  }else{
    row.files[[row.no]]
  }
})

## Check there are the required number of rows (and therefore files) for
## paired-endedness

# single-end runs
se.runs <- unique(sdrf[[run.col]][!grepl("PAIRED", sdrf[[library.layout.col]],ignore.case=TRUE)])

# paired-end runs
pe.runs <- unique(sdrf[[run.col]][grepl("PAIRED", sdrf[[library.layout.col]],ignore.case=TRUE)])

n.paired.run.files <- table(sdrf[[run.col]][sdrf[[run.col]] %in% pe.runs])

if(any(n.paired.run.files != 2)){
  perror(paste('Paired runs', paste(names(n.paired.run.files[which(n.paired.run.files != 2)]), collapse=', '), 'did not have two run files'))
  q(status=1)
}

## Finally, collapse to files per run

run.fastq.files <- lapply(unique(sdrf[[run.col]]), function(run){
  
  # If there are multiple rows, then this is paired-end. Sort so that first
  # read file comes first. Otherwise maintain ordering

  fq.files <- unlist(row.files[sdrf[[run.col]] == run])
  
  if ( length(which(sdrf[[run.col]] == run)) > 1 ){
    sort(fq.files)  
  }else{
    fq.files
  }
})
names(run.fastq.files) <- unique(sdrf[[run.col]])

################################################################################
# We'll be making a config per species and protocol, and these sub-experiments
# will then be analysed separately. So checks that with an experimental context
# (presence of technical replication, consistency of spikeins etc) need to be
# done within that species.
################################################################################

sdrf.by.species.protocol <- lapply(split(sdrf, sdrf[[organism.col]]), function(x) split(x, x[[protocol.col]]) )

# Filter any cells file to match the protocol-wise data
cells.by.species.protocol <- list()

# Run the checks first

species.protocol.properties <- list()

for (species in names(sdrf.by.species.protocol)){

  for (protocol in names(sdrf.by.species.protocol[[species]])){

    species.protocol.sdrf <- sdrf.by.species.protocol[[species]][[protocol]]
  
    # Some default propertiesa
      
    properties <- list(
      protocol = tolower(protocol),
      has.spikes = FALSE,
      has.techreps = FALSE,
      has.strandedness = FALSE,
      has.controlled.access = FALSE,
      has.cell.meta = FALSE
    )
  
    # Filtering could have removed all rows for a species
  
    if (nrow(species.protocol.sdrf) == 0){
      perror(paste0('No rows remaining for ', species, ', check filtering'))
      q(status=1)
    }
  
    ## Spikes? Assert that there are no spike-ins unless the field is populated with
    ## non-empty values of a single type.
  
    if (! is.null(spike.in.col)){
      print.info("Found 'spike in' column")
      spikein <- unique(species.protocol.sdrf[[spike.in.col]])
      if ( length(spikein) > 1 ) {
        perror("SDRF error: Expected one or no spike ins (found ",length(spikein),").")
        q(status=1)
      }else if (populated(spikein)){
        properties$has.spikes <- TRUE
      }
    }
  
    ## Check for technical replicates . Consider that there are no tech rep groups
    ## if:
    ##   1) all values are empty ('' or NA) 
    ##   or 2) there is a single value
  
    nruns <- length(unique(species.protocol.sdrf[[run.col]]))
 
    if (! is.null(techrep.col)){ 
      tr <- species.protocol.sdrf[[techrep.col]]
      n.techrep.groups <- length(unique(tr))
    
      is.empty <- tr == "" | is.na(tr) | tr == 'NA' | tr == 'not applicable'
    
      if ( (! all(is.empty)) && any(is.empty) ){
        perror("SDRF: Found ",sum(is.empty)," entries in technical replicate group without values where some values are set")
        q(status=1)
      }else if(n.techrep.groups == 1 ||  n.techrep.groups == nruns){
      
        # If the technica replicate field has been set, but the group numbers don't
        # indicate technical replication (i.e. every run has its own group, or
        # there's only only 1), we would normally say that there is in fact no
        # technical replication. But since droplet protocols have multiple cells per
        # run we'll give these entries the benefit of the doubt.
      
        if (is.singlecell && all(is.droplet.protocol(unique(species.protocol.sdrf[[protocol.col]])))){
          properties$has.techrep <- TRUE
        }
      }else{
        properties$has.techrep <- TRUE
      }
    }
  
    # Check for strandedness

    if ( ! is.null(strand.col)){
      strandeness <- species.protocol.sdrf[[strand.col]]

      if (any(tolower(strandeness) != 'not applicable')){
          properties$has.strandedness <- TRUE
      }  
    }
  
    # ENA sample. Check that the sample names don't suggest we should have technical
    # replicates when we don't
  
    if (! is.null(ena.sample.col)){
      nsamples <- length(unique(species.protocol.sdrf[[ena.sample.col]]))
      has.techrep.from.samples <- nsamples != nruns
      if ( has.techrep.from.samples != properties$has.techrep  ) {
        addWarning(paste("Technical replicate group from", techrep.col, paste0('(',properties$has.techrep,')'), " inconsistent with ena_sample", paste0('(',has.techrep.from.samples,')')))
      }
    }
  
    ## Check if there are any controlled access runs

    if ( ! is.null(controlled.access.col) && any(sdrf[[controlled.access.col]] == 'yes')){
       properties$has.controlled.access=TRUE
    }

    # Check if there are inferred cell types in the SDRF

    if ( (! is.null(cell.meta.cols))){
        properties$has.cell.meta=TRUE
    }

    species.protocol.properties[[species]][[protocol]] <- properties 
  }
}

# With all checks done, we can stop if no further output is required

pinfo("SDRF/IDF looks OK!")    

if ( opt$check_only ) {
  printAllWarnings()
  pinfo("<<< All done: SDRF validation complete, exiting before config creation")
  q(status=0)
}

# Now generate config and metadata

# Silly device to make sure the list output in the subsetquent step is named correctly
species_list <- as.list(names(sdrf.by.species.protocol))
names(species_list) <- unlist(species_list)

configs <- lapply(species_list, function(species){

  protocol_list <- as.list(names(sdrf.by.species.protocol[[species]]))
  names(protocol_list) <- unlist(protocol_list)
 
  lapply( protocol_list, function(protocol){
 
    species.protocol.sdrf <- sdrf.by.species.protocol[[species]][[protocol]]
    properties <- species.protocol.properties[[species]][[protocol]]
 
    # layout is the SDRF with the duplicated lanes for paired end (and possibly technical replication) removed
    if(properties$has.techrep) {
      species.protocol.layout <- species.protocol.sdrf[!duplicated(species.protocol.sdrf[[techrep.col]]),,drop=FALSE]
      rownames(species.protocol.layout) <- species.protocol.layout[[techrep.col]]
    }else{   
      species.protocol.layout <- species.protocol.sdrf[!duplicated(species.protocol.sdrf[[run.col]]),,drop=FALSE]
      rownames(species.protocol.layout) <- species.protocol.layout[[run.col]]
    }

    # Generate starting config file content

    config <- c(
      "\nparams{",
      paste0("    name = '", opt$name, "'"),
      paste0("    organism = '", species, "'"),
      paste0("    protocol = '", protocol, "'")
    )

    config_fields <- c(run = run.col, layout = library.layout.col)
    
    if (!  is.droplet.protocol(protocol)){
      if ( length(fastq.fields) > 1 ){
        perror('Multiple fastq fields on non-droplet experiment')
        q(status=1)
      }

      # For non-droplet HCA experiments...
      if (is.hca){
        species.protocol.sdrf[['hca_uri']] <- paste('hca', species.protocol.sdrf[[hca.bundle.uuid.col]], species.protocol.sdrf[[hca.bundle.version.col]], species.protocol.sdrf[[fastq.fields]], sep='/')
        config_fields['fastq'] <- 'hca_uri'
      }else{
        config_fields['fastq'] <- fastq.fields
      }
    }

    ## Field to use for quality filtering

    if ( ! is.null(sc.quality.col)){
      config_fields['quality'] <- sc.quality.col 
    }

    ## Anything other than ERCC will be ignored

    spikes <- c()

    if (properties$has.spikes){
      spikein <- tolower(unique(species.protocol.sdrf[[spike.in.col]]))
      if ( grepl("ercc.*",spikein) ) {
        config_fields['spike'] <- spike.in.col
        spikes <- paste0("    spikes = 'ercc'"  )
      }else{
        print(paste("ignoring", spikein, 'for spikein'))
      }
    }
  
    ## Set up for technical replication (or not)
  
    if(properties$has.techrep) {
      config_fields['techrep'] <- techrep.col
    }

    # Do any rows need stranded analysis

    if (properties$has.strandedness){
      config_fields['strand'] <- strand.col
    }

    # Do any rows need controlled access analysis?

    if (properties$has.controlled.access){
      config_fields['controlled_access'] <- controlled.access.col
    }
    
    # For droplet techs, add colums with strighforward statements of the URIs
    # that contain barcodes and cDNAs, and record which columns to use (unless
    # HCA in which case the UUID added above will work
        
    if ( is.droplet.protocol(protocol)){
 
      cdna_field <- getActualColnames('cdna read', sdrf)
      umi_field <- getActualColnames('umi barcode read', sdrf)
      cb_field <- getActualColnames('cell barcode read', sdrf)
      
      if (! is.hca){
        uri_cols <- which(colnames(sdrf) == fastq.col)
        if (length(uri_cols) < 2 ){
          perror('Less than 2 FASTQ URI fields supplied for droplet experiment- expect at least two, probably one for barcode/UMI, one for cDNA.')
          q(status=1)
        }    
      }

      # Right now we need umi and cell barcodes to be in the same file. Might be
      # different when we start to handle 10xv1

      if ( length(umi_field) == 0 || length(cb_field) == 0 ){
        perror('UMI barcode read or cell barcode read fields not supplied')
        q(status=1)
      }

      if ( any(sdrf[[umi_field]] != sdrf[[cb_field]]) ){
        perror("Cell barcodes and UMIs must be in the same file for currently enabled droplet protocols")
        q(status=1)
      }

      # Work out which URI holds the cDNA reads. We have to work out which read
      # number has the cdna reads, then check the 'readN file' column for the
      # file. For the full URI we then have to work out which of the FASTQ URI
      # fields contains that file. Things are simpler for HCA, where we
      # just use the file field (downstream workflows then use that in
      # combination with the UUIDs to get the files.

      for (field_type in c('cdna', 'cell barcode', 'umi barcode')){
        read_field <- getActualColnames(paste(field_type, 'read'), species.protocol.sdrf)
        uri_field <- gsub(' ', '_', paste(field_type, 'uri'))
        file_field_name <- paste(sub(' ', '', species.protocol.sdrf[[read_field]]), 'file')
        file_fields <- getActualColnames(file_field_name, species.protocol.sdrf)

        if (is.null(file_fields)){
          perror(paste(file_field_name, 'field not found in SDRF'))
          q(status=1)
        }
            
        files <- unlist(lapply(1:nrow(species.protocol.sdrf), function(x) species.protocol.sdrf[x, file_fields[x]]))

        if (is.hca){
          species.protocol.sdrf[[uri_field]] <- paste('hca', species.protocol.sdrf[[hca.bundle.uuid.col]], species.protocol.sdrf[[hca.bundle.version.col]], files, sep='/')
        }else{

          nlibs <- nrow(species.protocol.sdrf)
          uri_select <- apply(species.protocol.sdrf[,uri_cols], 2, function(x) basename(x) == files)
          
          if (nlibs > 1){
            uri_fields <- uri_cols[apply(uri_select, 1, function(x) which(x))]
          }else{
            uri_fields <- uri_cols[uri_select]
          }
          species.protocol.sdrf[[uri_field]] <- unlist(lapply(1:nrow(species.protocol.sdrf), function(x) species.protocol.sdrf[x, uri_fields[x]]))      
        }  
        config_fields[uri_field] <- uri_field
      }

      # Record the barcode position and offset fields
     
      droplet_fields_for_analysis <- c(apply(expand.grid(c('cdna read', 'cell barcode', 'umi barcode'), c('offset', 'size')), 1, function(x) paste(x, collapse=' ')), 'end')

      for (dffa in droplet_fields_for_analysis){
        field_label = gsub(' ', '_', dffa)
        field_name = getActualColnames(dffa, sdrf)

        if (is.null(field_name)){

          if ( protocol %in% names(droplet.protocol.defaults) && field_label %in% names(droplet.protocol.defaults[[protocol]]) ){

              # We can populate a field with the default value for the protocol if necessary

              species.protocol.sdrf[[field_label]] <- droplet.protocol.defaults[[protocol]][[field_label]]
              config_fields[field_label] <- field_label
          }else{
            perror(paste(field_label, 'field not present, and default not known for protocol', protocol))
            q(status=1)
          }

        }else{
          config_fields[field_label] <- field_name
        }
      }
      
    }

    # Choose the info we'll be putting in the cells metadata file. This file
    # will mostly be used to determine when analysis-relevant metdata has
    # changed between runs. In the case of a droplet experiment without a cells
    # file, the file of run IDs will not be cell identifiers per se, but will
    # serve to differentiate this run from a future one where that file is
    # present.
    
    if ( is.droplet.protocol(protocol) && ! is.null(opt$cells_file)){
      cell_run_techrep_ids <- unlist(lapply(strsplit(cells[[cell.id.col]], '-'), function(x) x[1]  ))
      if(properties$has.techrep) {
        cell.relate.col <- techrep.col    
      }else{
        cell.relate.col <- run.col
      }

      # This is a droplet protocol, expect cell type info to be in the cells
      # file. We relate this back to our protocol-wise SDRF dataframe to match
      # cells to protocol
   
      if ( ! is.null(batch.col)){
        cells[[batch.col]] <- sdrf[match(cell_run_techrep_ids, species.protocol.sdrf[[cell.relate.col]]), batch.col ]
      }

      cells.by.species.protocol[[species]][[protocol]] <<- cells[cell_run_techrep_ids %in% species.protocol.sdrf[[cell.relate.col]], cell.meta.cols, drop = FALSE]
    }else{
      # This is not a droplet protocol, or does not have a cells file, expect
      # cell type info to be in the SDRF file
      cells.by.species.protocol[[species]][[protocol]] <<- species.protocol.sdrf[, cell.meta.cols, drop = FALSE]        
    }    

    # Record field containing cell counts, where present
    
    if ( is.null(cell.count.col)){
      cell.count.col <- 'cell count'
      species.protocol.sdrf[[cell.count.col]] <- NA
    }
    config_fields['cell_count'] <- cell.count.col

    # Re-save the tweaked SDRF for output
    sdrf.by.species.protocol[[species]][[protocol]] <<- species.protocol.sdrf[, config_fields]        

    # Add cell metadata fields we need to know about, but which didn't need to
    # go in the tweaked SDRF (i.e. they're relevant for tertiary analysis)

    if ( ! is.null(batch.col)){
        config_fields['batch'] <- batch.col
    }

    if ( ! is.null(cell.type.col)){
        config_fields['cell_type'] <- cell.type.col
    }
   
    # Create the config fields section

    config <- c(
      config,
      paste0("\n    fields {"),
      unlist(lapply(names(config_fields), function(x) paste0("        ", x," = '", config_fields[x], "'"))),
      '    }\n'
    )
    
    # Put spikes in if present

    config <- c (config, spikes)
 
    # SC prototol
  
    sc.prot.conf = c()
    sc.clusters.conf = c()
    platform.config = c() 
 
    if (is.singlecell && opt$nc > 0){

      ## Decide what range of K values to try. Use a window around the provided
      ## central k
    
      # Guess at a central K when not provided
    
      if ( opt$nc < 2 ) {
        opt$nc <- round(log2(nrow(species.protocol.layout)))
      }
    
      num.min.c <- max(round(opt$nc - opt$nc_window/2,0),2)
      num.max.c <- num.min.c + opt$nc_window
    
      config <- c(config, 
          paste0('\n    clusters {'),
          paste0('        min = ', num.min.c),
          paste0('        max = ', num.max.c),
          paste0('    }')
      )
    }
  
    config <- c(paste("// Generated from", opt$sdrf), config, paste0("\n    extra_metadata = '", paste0(opt$name, ".metadata.tsv'")), "}\n")

    # Generate metadata file content to save
  
    metadata <- NULL

    metadata_cols <- unique(c(factors, idf.factors, idf.attrs, techrep.col))
    sdfcols2save2tsv <- getActualColnames(metadata_cols, species.protocol.sdrf)
  
    if (length(sdfcols2save2tsv) > 0){
      metadata <- species.protocol.layout[,sdfcols2save2tsv, drop = FALSE]
      metadata <- cbind(run=rownames(metadata),metadata)
      colnames(metadata) <- c('run', metadata_cols)
    }

    list(config=config, metadata=metadata)
  })
})

if ( length(configs)==0 ) {
  perror("Unable to create configuration file, see above errors")
  q(status=1)
}

## Write final outputs

cat("\n")
pinfo("### SDRF summary info\n")
cat("Columns: ",ncol(sdrf),"\n")
cat("Entries: ",nrow(sdrf),"\n")
pinfo("Columns found in SDRF:")
cat(paste(colnames(sdrf),collapse=","),"\n")

## print the values
pinfo("Most frequent values per column...")
for (col in tolower(c(expected.cols,expected.comment.cols,expected.factor.cols,expected.characteristic.cols)) ) {
  actualCol <- getActualColnames(col, sdrf)  
  if (! is.null(actualCol)){
    tt <- table(sdrf[,actualCol])
    if (length(tt)==nrow(sdrf) ) {
      pinfo(col,": single value per row")
    } else {
      pinfo(col,": ",length(tt)," unique values")
    }
  }
}

cat("\n")
pinfo("### Warnings\n")
printAllWarnings()

cat("\n")
pinfo("### Outputs\n")

if ( is.null(opt$out_conf) ){
    opt$out_conf <- getwd()
}

for (species in names(configs)){

  protocol_configs <- configs[[species]]

  for ( protocol in names(protocol_configs)){

    # If there's only one protocol, don't add a protocol suffix to the the outputs

    file_prefix <- paste(protocol, species, sep='.')

    conf.file <- file.path(opt$out_conf, paste0(file_prefix,'.',  opt$name, ".conf"))
    meta.file <- file.path(opt$out_conf, paste0(file_prefix, '.', opt$name, ".metadata.tsv"))

    # These two files subset the metadata to that which inpacts on
    # quantification and on tertiary analysis. Changes in these outputs between
    # runs can be used to detect if changes to metadata require re-analysis
    
    sdrf.file <- file.path(opt$out_conf, paste0(file_prefix, '.', opt$name, ".meta_for_quant.txt"))
    cells.file <- file.path(opt$out_conf, paste0(file_prefix, '.', opt$name, ".meta_for_tertiary.txt"))
  
    config <- configs[[species]][[protocol]]$config
    metadata <- configs[[species]][[protocol]]$metadata

    dir.create(opt$out_conf, showWarnings = FALSE)
 
    # If there is metadata then point to it from the config file
  
    if ( ! is.null(metadata)){
      write.tsv(metadata, file=meta.file)
      pinfo("Created ", meta.file)
    }
    write.tsv(sdrf.by.species.protocol[[species]][[protocol]], file=sdrf.file)
    write.tsv(cells.by.species.protocol[[species]][[protocol]], file=cells.file)
    writeLines(config, con = conf.file)
    pinfo("Created ", conf.file)
    pinfo("Created ", cells.file)
  }
}

cat("\n")
pinfo("<<< All done: SDRF validation complete")
q(status=0)
