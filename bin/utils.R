pinfo <- function(...) {
    cat(paste0("[INFO ",format(Sys.time(), format="%d/%m-%H:%M"),"] ",...,"\n"))
}

perror <- function(...) {
  cat(paste0("[ERROR] ",...,"\n"),file=stderr())
}

pwarning <- function(...) {
  cat(paste0("[WARNING] ",...,"\n"),file=stderr())
}

read.tsv <- function(f,header=TRUE, comment.char="", nrows=-1L,fill=FALSE, quote="\"", colClasses=NULL, drop=NULL) {
    tsv.data <- NULL

    suppressPackageStartupMessages(require(data.table))
    
    fread(input=f,sep = "\t", nrows=nrows, header=header,check.names=FALSE,data.table=FALSE,fill=fill,drop=drop, quote=quote,colClasses=colClasses)
}

write.tsv <- function(x, file, header=TRUE){
    fwrite(x,file=file,sep="\t",row.names=FALSE,col.names=header,quote=FALSE,verbose=FALSE)
}
