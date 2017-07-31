#'get a descriptor (fingerprints and/or CDK physical descriptors) from SMILES strings
#' @description Get a descriptor (fingerprints and/or CDK physical descriptors) from SMILES strings with possbility
#' to request the scaling (for continuous descriptors, e.g. physical) or re-casting (for binary descriptors, e.g. fingerprints) of
#' the output descriptors.
#' @param smis is a SMILES strings vector ("C1=CC=C(C=C1)O", canonical SMILES of a phenol by default).
#' @param desctypes is a vector of characters defining the fingerprints and/or physical descriptors types to compute ("standard" by default).
#' The actual entire list of available fingerprints: "standard", "extended", "graph", "hybridization", "maccs", "estate", "pubchem", "kr",
#' "shortestpath" and "circular", and physical descriptors: "constitutional","topological","electronic" can be computed.
#' @param scale sets to TRUE (FALSE by default) for scaling the physical descriptors only (i.e. continuous features) - mean = 0, s.d. = 1.
#' @param scale_init sets to TRUE (FALSE by default) to keep in memory the means and s.d. related to each descriptor after a first scaling. Indeed,
#' after the descriptors on a training set have been first computed, the mean and s.d. have to be kept fixed for future descriptors
#' computation on test and/or validation sets. In this last case, the scale_init variable is set to FALSE.
#' @param mdesc is a scalar (0 by default) or vector of means for a post-scaling of physical descriptors.
#' @param sddesc is a scalar (1 by default) or vector of standard deviations for a post-scaling of physical descriptors.
#' @param quiet keeps the console's outputs quiet if sets to TRUE (FALSE by default).
#'
#' @examples \dontrun{
#'
#' descriptors <- get_descriptor(smis = "C1=CC=C(C=C1)O", desctypes = c("standard","topological"))
#'
#' }
#'
#' @return the descriptor(s) with the associated means and standard deviations for scaling.
#'
#' @import fingerprint
#' @import rcdk
#' @import stats
#'
#' @export get_descriptor

get_descriptor <- function(smis = c("C1=CC=C(C=C1)O"),
                           desctypes = c("standard"),
                           scale = F, scale_init = F,
                           mdesc = 0, sddesc = 1,
                           quiet = F){

  if(is.null(smis)){
    cat("The vector of input SMILES you entered is empty\n\n")
    return(NULL)
  }

  fp <- NULL
  pdesc <- NULL
  pdesc_types <- c("constitutional","topological","electronic")

  fpdesc <- NULL

  if(mdesc==0 && sddesc==1 && scale==T && scale_init==F){
    cat("Please perform the same job with scale_init = T\n")
    cat("You have requested a re-scaling without pre-computing the mean and standard deviation.\n\n")
  }

  if(!quiet)cat("Parsing of the SMILES strings in progress...\n")
  mols <- parse.smiles(as.character(smis), kekulise=F) # kekulise on FALSE, the check should be done a priori
  mols_op <- lapply(mols,do.aromaticity)
  mols_op <- lapply(mols,do.typing)
  mols_op <- lapply(mols,do.isotopes)
  rm(mols_op)
  if(!quiet)cat("Parsing of the SMILES strings done.\n")

  fpdescname <- c()
  physdescname <- c()
  for(descname in desctypes){
    if(!descname%in%pdesc_types){
      fpdescname <- c(fpdescname,descname)
    } else {
      physdescname <- c(physdescname,descname)
    }
  }

  # Trick to check the fingerprints sizes to manage
  # Useful to name the columns and retrieve relevant features later on
  mol_checker <- parse.smiles("CC", kekulise=F)
  ##

  if(length(fpdescname)!=0){
    if(!quiet){
      cat("get the", fpdescname ,"fingerprint(s) from SMILES strings\n")
    }
    fp.obj <- list()
    fp_temp <- c()
    fp_temp_colnames <- c()
    for(fpdescname_i in fpdescname){
      # Compute the fingerprint of type fpdescname
      fp.obj <- lapply(mols, get.fingerprint, type=fpdescname_i)
      fp_temp <- fp.to.matrix(fp.obj)
      fp <- cbind(fp,fp_temp)
      ##
      fp_temp_dim_tmp <- dim(fp.to.matrix(lapply(mol_checker, get.fingerprint, type=fpdescname_i)))
      fp_temp_colnames <- c(fp_temp_colnames,paste(fpdescname_i,"_",seq(1,fp_temp_dim_tmp[2],1), sep=""))
    }
    colnames(fp) <- fp_temp_colnames
    rm(fp.obj,fp_temp, fp_temp_colnames)
  }
  rm(mol_checker)

  if(length(physdescname)!=0){
    if(!quiet){
      cat("get the", physdescname,"physical descriptors from SMILES strings\n")
    }
    descriptors <- c()
    for(desc_i in physdescname){
      descriptors <- c(descriptors,get.desc.names(type=desc_i)) # access to the different descriptors of type physdescname
    }

    # Computation of the physical descriptors for the type physdescname
    desc.list = list() # initialize a list
    for (i in 1:length(descriptors)) {  # loop through descriptors
      desc.list[[i]] = eval.desc(mols, descriptors[i]) # calculate descriptor - each one returns a vector
    }
    desc.data = as.data.frame(do.call("cbind", desc.list)) # convert to data frame
    desc.nannames.veto <- which(colnames(desc.data) %in% "geomShape")
    if(length(desc.nannames.veto)!=0){
      desc.data <- desc.data[,-desc.nannames.veto]
    }
    # handle the non-availability or divergence of a descriptor for a given molecule
    desc.data[is.na(desc.data)] <- 0
    #desc.data[is.infinite(desc.data)] <- 0
    ##
    pdesc <- as.matrix(desc.data)
    rm(desc.data)
    rm(desc.list)

    # pdesc = matrix(0,nrow = length(mols), ncol = length(descriptors)) # initialize a matrix
    # for (i in 1:length(descriptors)) {  # loop through descriptors
    #   pdesc[,i] <- eval.desc(mols, descriptors[i]) # calculate descriptor - each one returns a vector
    # }
    ## pdesc <- as.matrix(eval.desc(mols, descriptors)) # calculate descriptor - each one returns a vector
    ##
    # handle the non-availability or divergence of a descriptor for a given molecule
    # pdesc[is.na(pdesc)] <- 0
    ##
  }
  rm(mols)

  if(scale){
    if(!is.null(pdesc)){
      if(scale_init){
        mdesc <- apply(pdesc,2,mean, na.rm = TRUE)
        sddesc <- apply(pdesc,2,sd, na.rm = TRUE)
        scaledesc <- scale(pdesc)
      }else{
        scaledesc <- t((t(pdesc)-mdesc)/sddesc)
      }
      scaledesc[is.na(scaledesc)] <- 0         # prevent NA from NULL columns
      scaledesc[is.infinite(scaledesc)] <- 0   # prevent Inf if sddesc == 0
      pdesc <- scaledesc
    } else {
      cat("You have requested a scaling on the physical descriptors but their matrix is NULL.\n\n")
    }
  }

  fpdesc <- cbind(fp, pdesc)
  intercept <- matrix(1,nrow=dim(fpdesc)[1],ncol=1)
  colnames(intercept) <- c("Intercept")
  fpdesc <- cbind(intercept,fpdesc)

  return(list(fpdesc,mdesc,sddesc))
}

#' get a list of the available descriptors types
#' @description get a list of the available descriptors types
#'
#' @examples
#' list_descriptors <- get_descriptors()
#' @return a list of the available descriptors types
#'
#' @export get_descriptors
get_descriptors <- function(){
  descriptors_list = c("standard",
                  "extended",
                  "graph",
                  "hybridization",
                  "maccs",
                  "estate",
                  "pubchem",
                  "kr",
                  "shortestpath",
                  "signature",
                  "circular",
                  "topological",
                  "constitutional",
                  "electronic")

  return(descriptors_list)
}

