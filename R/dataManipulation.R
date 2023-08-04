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

#' Approximate connectivity by voronoi tesselation in 2d
#' @param x the x locations
#' @param y the y locations
#' @param width the periodic width
#' @param b_dist the boundary distance
#' @export
getVoronoiConnectivity2d <- function(x,y,width,b_dist=2.0)
{
  ghost_nodes = list( x=c(x[x<=b_dist]+width,x[x >= (width-b_dist)]-width),
                      y=c(y[x<=b_dist],y[x >= (width-b_dist)]) )
  tess_data <- deldir(x,y,dpl=ghost_nodes)
  output <- list(connectivity = tess_data[["summary"]][tess_data[["summary"]]$pt.type=="data",c("x","y","n.tside")],
                 voronoi_lines = tess_data[["dirsgs"]][,c("x1","y1","x2","y2")],
                 neighbours = tess_data[["delsgs"]],
                 ghost_nodes=do.call(cbind,ghost_nodes))
  return(output)
}


#' Approximate connectivity by voronoi tesselation
#' @param x the x locations
#' @param y the y locations
#' @param z the z locations
#' @param width the periodic width
#' @param b_dist the boundary distance
#' @export
getVoronoiConnectivity3d <- function(x,y,z,width,b_dist=2.0)
{
  # Create an appropriate matrix
  d <- as.matrix( cbind(x,y,z) )
  output <- delaunayn(d)
  return(output)
}

#' Convert a named vector from time steps to data frame with days
#' @param v the named vector
#' @param value_name the string part of the names to gsub out
#' @export
convertNamedTimeStepVector <- function(v, value_name="v")
{
  d <- data.frame(Day=as.numeric(gsub("TimeStep_","",names(v)))/24)
  d[,value_name] <- as.vector(v)
  return(d)
}

#' Load in a results.vizancestors file
#' Reading in of the cell ancestory output
#' @param filepath the path to the containing results folder
#' @export
getCellAncestory<-function(filepath)
{
  # Read in ancestors and determines the colonies
  ancestors<-readChasteResultsFile(paste(filepath,"results.vizancestors",sep="/"))
  n_colonies <- lapply(ancestors,table)
  n_colonies <- melt(n_colonies)

  # Tidy up the dataset
  colnames(n_colonies) <- c("Ancestor","NumberCells","TimeStep")
  n_colonies$Ancestor = factor(paste("Cell",n_colonies$Ancestor))#,levels=paste("Cell",as.character(5:0)))
  n_colonies$TimeStep = gsub("TimeStep_","",n_colonies$TimeStep)
  n_colonies$Day = as.numeric(n_colonies$TimeStep)/24

  # Get the number of colonies over time
  counts_colonies <- lapply(ancestors,unique)
  counts_colonies <- sapply(counts_colonies,nrow)
  counts_colonies <- data.frame( Day=as.numeric(gsub("TimeStep_","",names(counts_colonies)))/24, Colonies=as.numeric(counts_colonies)  )

  return(list(raw=ancestors,colony_size=n_colonies,num_colonies=counts_colonies))
}

#' Generate the hybrid force functions given peak, location of peak, and cutoff
#' @param x the x locations
#' @param max_F the peak force
#' @param peak_location the peak location
#' @export
calcHybridForce<-function(
    x = seq(-0.3,2.5,0.01),
    max_F = 0.2, peak_location = 1.1)
{
  alpha = 1/peak_location
  y <- data.frame(x,force=0,fn="Repulsion")
  # Calculate repulsion
  spring_stiffness = max_F * alpha * exp(1.0)
  y$force[x < 0] = spring_stiffness * log(1.0 + x[x < 0])
  # Calculate adhesion/attraction
  y$force[x > 0] = max_F * alpha * x[x>0] * exp(1.0 - alpha * x[x>0])
  # Add an extra 0 row for visualisation purposes
  y$fn[x > 0] = "Adhesion"
  if (sum(x==0) > 0) {
    y <- rbind(y,data.frame(x=0,force=0,fn="Adhesion"))
  }
  y <- y[order(y$x),]
  return(y)
}



