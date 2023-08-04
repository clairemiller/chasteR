#' Efficient loading of a Chaste results file
#'
#' Reading in of the generic results.viz[data] file with timesteps
#' @param filename the name or filepath of the results file
#' @param filepath the path to the file
#' @param columns (optional) a vector of names for the columns of the outputs for each time step
#' @return A list of data frames, one for each time step. If columns are specified data is split into these columns.
#' @export
readChasteResultsFile <- function(filename, filepath=".", columns="v")
{
  # New read method using readChar:
  fullfilename <- file.path(filepath,filename)
  size <- file.info(fullfilename)$size
  dat <- readChar(fullfilename,size,useBytes=T)
  dat <- unlist(strsplit(dat,"\r\n|\n|\r"))

  # Split the string into the desired columns
  dat <- strsplit(dat,"\\s+")
  # Remove any empty rows
  i <- which(sapply(dat,length) > 1)
  dat <- dat[i]
  dat <- lapply(dat, function(x) {
    t <- as.numeric(x[1]) # Extract the time
    d <- as.numeric(x[-1]) # Convert data to numeric
    d <- matrix(d,ncol=length(columns),byrow=TRUE,dimnames = list(rows=c(),columns))
    return( data.frame(time=t,d) )
  })
  dat <- do.call(rbind,dat)
  return(dat)
}

#' Load cell locations file
#'
#' Read in the results.viznodes file and split into locations
#' @param filepath the path to the containing results folder
#' @param num_dims the number of dimensions of the simulation
#' @return A list of data frames, one for each time step, containing the cell locations.
#' @export
readCellLocations<-function(filepath,num_dims=3)
{
  columns <- c("x","y","z")[1:num_dims]
  cell_locations <- readChasteResultsFile("results.viznodes",filepath=filepath,columns=columns)
  return(cell_locations)
}


#' Read a cell data file
#'
#' Loads in data from a cell data output
#' @param var_name the variable name of the cell data from chaste
#' @return A data frame containing the cell data output
#' @export
readCellData<-function(filepath,var_name, dimensions = 3)
{
  filename <- paste0("celldata_",var_name,".dat")
  dim_names <- c("x","y","z")[1:dimensions]
  columns <- c("loc_idx","id",dim_names,var_name)
  return(readChasteResultsFile(filepath=filepath,filename=filename,columns=columns))
}

#' Read a cell velocities file
#'
#' Loads in data from NodeVelocityWriter output
#' Note, has got my modification where the cell id is not repeated
#' @param filepath the folder path to where the simulation output is located
#' @param dimensions the number of dimensions of the simulation
#' @export
readCellVelocities<-function(filepath, dimensions = 3)
{
  filename <- "nodevelocities.dat"
  dim_seq <- 1:dimensions
  columns <- c("id",c("x","y","z")[dim_seq],c("vx","vy","vz")[dim_seq])
  return(readChasteResultsFile(filepath=filepath,filename=filename,columns=columns))
}
