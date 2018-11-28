#' Loop a function over a set of results directories
#' @param fun the function to apply
#' @param maindir the parent directory
#' @export
applyFunctionToResultsDirectories <- function(fun,maindir="./",pattern="",...)
{
  subdirs <- grep(pattern,list.dirs(path=maindir),value=T)
  subdirs <- as.list(grep("results_from_time_",subdirs,value=T))
  names(subdirs) <- gsub("\\.|/|results_from_time_[0-9]+","",subdirs)
  return( lapply(subdirs,fun,...) )
}

#' Calculate the counts of a list of time step data frames
#' @param data a list of the chaste results data at each time step
#' @param column the column of the data frame to count over
#' @return a data frame containing the counts
#' export
getCellCounts <- function(data,column)
{
  d <- lapply(data,function(x){table(as.factor(data[,column]))})
  d <- lapply(d,as.data.frame)
  d <- rbindlist(d,use.names=TRUE,idcol=TRUE)
  d$time <- as.numeric(gsub("TimeStep_","",d$.id))
  return( d )
}
