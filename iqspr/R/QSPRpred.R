#'QSPRpred class
#' @description Quantitative Structure-Properties Relationship (QSPR) model construction. This class contains all the required
#' functions to train linear and non-linear models, to produce bootstrap datasets for variance estimation, and to provide prediction
#' capabilities over a matrix or vector of studied properties.
#'
#' @field propndim is the number of properties received as input data.
#' @field propmin is a vector representing the expected minimal value for each targeted property.
#' @field propmax is a vector representing the expected maximal value for each targeted property.
#' @field filtermin is a vector representing the expected minimal value for each filtered property.
#' @field filtermax is a vector representing the expected maximal value for each filtered property.
#' @field filterfunc is a function to compute the properties to filter.
#' @field X is the nxd matrix, with d features for n input SMILES, returned by \code{\link{get_descriptor}}.
#' @field Y is a nxp matrix of p properties for n input SMILES.
#' @field fnames is a list of vectors of fingerprints and/or physical descriptors types used as features in each regression model by
#' \code{\link{get_descriptor}}.
#' @field mdesc is a scalar or vector of means used for physical descriptors scaling, returned by \code{\link{get_descriptor}}.
#' @field sddesc is a scalar or vector of standard deviations used for physical descriptors scaling, returned by \code{\link{get_descriptor}}.
#' @field scale tags the scaling statement (TRUE or FALSE) of the physical descriptors only (i.e. continuous features) - mean = 0, standard deviation = 1.
#' @field func defines the analytic function to use in the computation of a subsequent property.
#' @field func_args is a vector of integers that tags the used columns in the property array prop for the computation of a subsequent property.
#' @field trmodel is the name of the used regression model for training and predictions.
#' @field trnboot is the number of bootstrap dataset used for the training.
#' @field trndf is the number of input SMILES, i.e. the number of degrees of freedom, available in the training of the regression process.
#'
#' @param smis is a list of vectors of SMILES from which a regression model will be trained, or for which targeted properties will be predicted.
#' @param prop is a list of vectors/matrices of available targeted physico-chemical properties for the training dataset.
#' @param v_filterfunc defines the filtering function (NULL by default) to use in the computation of properties to filter.
#' @param v_filtermin is a vector representing the expected minimal value for each filtered property.
#' @param v_filtermax is a vector representing the expected maximal value for each filtered property.
#' @param v_fnames is a vector, or a list of vectors, of fingerprints and/or physical descriptors types used as features for each regression model
#' (see \code{\link{get_descriptor}} for an exhaustive list of available descriptors).
#' @param v_scale sets (FALSE by default) the scaling of physical descriptors only (i.e. continuous features) - mean = 0, standard deviation = 1.
#' @param v_func defines the analytic function (NULL by default), or a list of analytic functions, to use in the computation of a subsequent property, or properties respectively.
#' A given function will return a new property computed analytically via a list of known properties in prop. This is particularly
#' useful when data and regression models can be stated for some properties (e.g. A and B), but not for a targeted property of
#' interest (e.g. A+B, A/B, etc.) for which constrains are defined via the \emph{set_target} method.
#' @param v_func_args is a vector, or a list of vectors, of integers that tags the used properties in prop for the computation of a subsequent property.
#' For example, v_func=list(func1,func2), where func1 and func2 are \emph{a priori} defined functions, and prop=list(V1,M23), where
#' V1 is a numerical vector and M23 is a two columns matrix. In this case, v_func_args=list(c(1,3),c(2)), i.e. the function func1
#' uses the 1st and 3rd output properties located in prop, and func2 uses the 2nd only. Therefore, the defined empirical functions
#' know where to find their inputs.
#' @param kekulise enables (FALSE by default) electron checking and allows for parsing of incorrect SMILES (see \code{\link[rcdk]{parse.smiles}}).
#'
#' @param model is the name of a regression model to be used (see \code{\link{get_Models}} for an exhaustive list).
#' @param params is a list of parameters to submit to a given regression model (see \code{\link{get_Model_params}} for
#' examples).
#' @param n_boot is the number of requested bootstrap datasets (1 by default) in the training process. This is used for
#' an estimation of the means and standard deviations of subsequent non-Bayesian predictions. A higher number of bootstrap
#' datasets will allow more accuracy in this estimation. However, it exists a trade-off between accuracy and computation
#' time that the user has to figure out. Consequently, in order to ease the bootstrap analysis, a parallelization
#' capability is implemented.
#' @param s_boot is the proportion of input data (0.85 by default), defined in ]0,1], used to construct bootstrap datasets.
#' @param r_boot allows (FALSE by default) the sampling in a bootstrap analysis to be performed with replacement.
#' @param parallelize allows (FALSE by default) to use the full computational capability of a user's machine for a bootstrap analysis.
#' Indeed, N-1 cores, with N the total number of cores available on the machine, will be used.
#'
#' @param v_propmin is a vector representing the expected minimal value for each targeted property.
#' @param v_propmax is a vector representing the expected maximal value for each targeted property.
#'
#' @param temp is a vector/matrix of numerical values which sets the initial temperatures in the annealing process for the
#' sequential Monte-Carlo sampler (see \code{vignette("tutorial", package = "iqspr")} for details).
#'
#' @examples \dontrun{
#'
#' # Load pre-existing data
#' data(qspr.data)
#' # Define input SMILES
#' smis <- paste(qspr.data[,1])
#' # Define associated properties
#' prop <- qspr.data[,c(2,5)]
#' # Define training set
#' trainidx <- sample(1:nrow(qspr.data), 5000)
#' # Initialize the prediction environment
#' # and compute fingerprints/descriptors associated to input SMILES
#' qsprpred_env <- QSPRpred()
#' qsprpred_env$initenv(smis=smis[trainidx], prop=as.matrix(prop[trainidx,]), v_fnames="graph")
#' # Train a regression model with associated parameters,
#' # number of bootstrapped datasets without CPUs parallelization
#' qsprpred_env$model_training(model="elasticnet",params=list("alpha" = 0.5),n_boot=10,parallelize=F)
#'
#' # Predict properties for a test set
#' predictions <- qsprpred_env$qspr_predict(smis[-trainidx])
#' # Plot the results
#' par(mfrow=c(1,2))
#' plot(predictions[[1]][1,], prop[-trainidx,1], xlab="prediction", ylab="true")
#' segments(-100,-100,1000,1000,col=2,lwd=2)
#' plot(predictions[[1]][2,], prop[-trainidx,2], xlab="prediction", ylab="true")
#' segments(-100,-100,1000,1000,col=2,lwd=2)
#'
#' # Set a targeted properties space
#' qsprpred_env$set_target(c(8,100),c(9,200))
#' # Predict properties for any input SMILES
#' # and their probability to be close to the targeted properties space
#' inv_pred <- qsprpred_env$qspr_predict(smis = smis[-trainidx], temp=c(3,3))
#'
#' See \code{vignette("tutorial", package = "iqspr")} for further options and details.
#'
#' }
#'
#' @import methods
#' @importFrom methods new
#' @import doParallel
#' @import tictoc
#'
#' @export QSPRpred
#' @exportClass QSPRpred

QSPRpred <- setRefClass("QSPRpred",
  fields = list(
    propndim = "numeric",
    propmin = "numeric",
    propmax = "numeric",
    filtermin = "numeric",
    filtermax = "numeric",
    filterfunc = "function",
    X = "list",
    Y = "list",
    fnames = "list",
    mdesc = "list",
    sddesc = "list",
    scale = "logical",
    func = "list",
    func_args = "list",
    filter = "list",
    filter_func = "list",
    trmodel = "character",
    trnboot = "numeric",
    trndf = "numeric"
  ),
  methods = list(
              ################################################
              ## Initialization of the QSPRpred environment ##
              ################################################
              init_env = function(smis=NULL, prop=matrix(0),
                                  v_filterfunc=NULL, v_filtermin=NULL, v_filtermax=NULL,
                                  v_fnames=NULL,
                                  v_scale=FALSE,
                                  v_func=NULL, v_func_args = NULL,
                                  kekulise = F){

                "initialize the QSPR predictor: implicitly called via the QSPRpred$new() method"

                if(is.null(smis)){
                  cat("Please, provide a vector (or a list of vectors) of characters as input SMILES strings smis.\n\n")
                  return(NULL)
                }

                if(class(prop)=="matrix" && length(prop)==1){
                  cat("Please, provide a numerical matrix (or a list of numerical matrices) as input properties prop.\n\n")
                  return(NULL)
                }

                list_descriptors <- get_descriptors()
                v_fnames_check <- unique(unlist(v_fnames))
                v_fnames_veto <- which(v_fnames_check%in%list_descriptors==F)
                if(length(v_fnames_veto)!=0){
                  cat("Caution:",v_fnames_check[v_fnames_veto],"descriptors type(s) is(are) not available.\n")
                  cat("Please, have a look to the full list of available descriptors types below.\n")
                  print(list_descriptors)
                  cat("\n\n")
                  return(NULL)
                }

                smis_err <- F
                if(class(smis)=="character"){
                  smis <- list(smis)
                } else if(class(smis)=="list"){
                  smis_elts_class = unique(as.character(lapply(smis,class)))
                  if(length(smis_elts_class)!=1){
                    cat("Input SMILES list smis contains objects from more than one class:\n",smis_elts_class,",\n")
                    cat("when this list should contain objects of character class only.\n\n")
                    smis_err <- T
                  } else{
                    if(smis_elts_class!="character"){
                      cat("Input SMILES list smis contains objects of",smis_elts_class,"class, not character.\n")
                      cat("Please, check it and/or cast them as character class: e.g. as.character().\n\n")
                      smis_err <- T
                    }
                  }
                } else {
                  cat("Input SMILES smis is neither a characters vector nor a list of characters vectors:\n",class(smis),"\n\n")
                  smis_err <- T
                }

                if(smis_err){
                  return(NULL)
                }

                prop_err <- F
                if(class(prop)=="matrix"){
                  prop <- list(prop)
                  prop_matrix_class <- class(unlist(prop))
                  if(prop_matrix_class!="numeric"){
                    cat("Output properties prop is not a matrix of numerical objects:\n",prop_matrix_class,"\n\n")
                    prop_err <- T
                  }
                } else if(class(prop)=="numeric"){
                  prop <- list(as.matrix(prop))
                } else if(class(prop)=="list"){
                  prop_elts_class = unique(as.character(lapply(prop,class)))
                  if(length(prop_elts_class)!=1){
                    if(length(prop_elts_class)==2){
                      if("matrix"%in%prop_elts_class&&"numeric"%in%prop_elts_class){
                        prop <- lapply(prop,as.matrix)
                      } else {
                        cat("Output properties list prop contains objects from class:\n",prop_elts_class,",\n")
                        cat("when this list should contain objects of matrix and/or numeric class only.\n\n")
                        prop_err <- T
                      }
                    } else {
                      cat("Output properties list prop contains objects from more than two classes:\n",prop_elts_class,",\n")
                      cat("when this list should contain objects of matrix and/or numeric class only.\n\n")
                      prop_err <- T
                    }
                  } else{
                    if(prop_elts_class=="numeric"){
                      prop <- lapply(prop,as.matrix)
                    } else if(prop_elts_class=="matrix"){
                      prop_matrix_class <- class(unlist(prop))
                      if(prop_matrix_class!="numeric"){
                        cat("Output properties list prop contains matrices with objects other than numerical:\n",prop_matrix_class,"\n\n")
                        prop_err <- T
                      }
                    } else {
                      cat("Output properties list prop contains objects of",prop_elts_class,"class, neither numeric nor matrix.\n\n")
                      prop_err <- T
                    }
                  }
                } else {
                  cat("Output properties prop is neither a numerical matrix nor a list of numerical matrices:\n",class(prop),"\n\n")
                  prop_err <- T
                }

                if(prop_err){
                  return(NULL)
                }

                if(length(smis)!=length(prop)){
                  cat("Caution: lists of input SMILES smis and output properties prop of different length are provided:\n")
                  cat(length(smis),"and",length(prop),"respectively.\n\n")
                  return(NULL)
                }

                prop_tmp <- list()
                smis_tmp <- list()
                na_null_tmp <- 0
                for(yi in 1:length(prop)){
                  prop_mat_n <- dim(prop[[yi]])[2]
                  if(prop_mat_n>1){
                    for(prop_mat_i in 1:prop_mat_n){
                      prop_mat_i_na <- which(is.na(prop[[yi]][,prop_mat_i]))
                      smis_null <- which(smis[[yi]]=="")
                      prop_mat_smis_veto <- unique(c(prop_mat_i_na,smis_null))
                      if(length(prop_mat_smis_veto)!=0){
                        prop_tmp <- c(prop_tmp,list(as.matrix(prop[[yi]][-prop_mat_smis_veto,prop_mat_i])))
                        smis_tmp <- c(smis_tmp,list(smis[[yi]][-prop_mat_smis_veto]))
                        na_null_tmp <- na_null_tmp + length(prop_mat_smis_veto)
                      } else {
                        prop_tmp <- c(prop_tmp,list(as.matrix(prop[[yi]][,prop_mat_i])))
                        smis_tmp <- c(smis_tmp,list(smis[[yi]]))
                      }
                      prop_mat_smis_veto <- NULL
                      smis_null <- NULL
                      prop_mat_i_na <- NULL
                    }
                  } else {
                    prop_mat_i_na <- which(is.na(prop[[yi]]))
                    smis_null <- which(smis[[yi]]=="")
                    prop_mat_smis_veto <- unique(c(prop_mat_i_na,smis_null))
                    if(length(prop_mat_smis_veto)!=0){
                      prop_tmp <- c(prop_tmp,list(as.matrix(prop[[yi]][-prop_mat_smis_veto])))
                      smis_tmp <- c(smis_tmp,list(smis[[yi]][-prop_mat_smis_veto]))
                      na_null_tmp <- na_null_tmp + length(prop_mat_smis_veto)
                    } else {
                      prop_tmp <- c(prop_tmp,list(prop[[yi]]))
                      smis_tmp <- c(smis_tmp,list(smis[[yi]]))
                    }
                    prop_mat_smis_veto <- NULL
                    smis_null <- NULL
                    prop_mat_i_na <- NULL
                  }
                }
                prop <- prop_tmp
                smis <- smis_tmp
                rm(prop_tmp)
                rm(smis_tmp)
                if(na_null_tmp!=0){
                  cat("Caution:",na_null_tmp,"input SMILES with their related output properties have been erased.\n")
                  cat("These input SMILES and/or output properties were respectively \"\" or NA (or NaN) in the submitted dataset(s).\n\n")
                }

                propndim <<- length(prop)
                cat(length(unique(unlist(smis))),"unique SMILES strings are submitted with",propndim,"related properties.\n\n")

                fnames <<- list()
                if(is.null(v_fnames)){
                  for(mi in 1:propndim){
                    fnames[[mi]] <<- c("standard")
                  }
                  cat("Standard fingerprint bits will be used by default as features for all the properties predictions.\n\n")
                }else{
                  v_fnames_class <- class(v_fnames)
                  if(v_fnames_class=="list"){
                    if(length(v_fnames)==1){
                      for(mi in 1:propndim){
                        fnames[[mi]] <<- v_fnames[[1]]
                      }
                      cat(v_fnames[[1]],"fingerprint/physical descriptors will be used as features for all the properties predictions.\n\n")
                    } else {
                      if(length(v_fnames)!=propndim){
                        cat("Caution: The number of the descriptors types vectors in v_fnames:",length(v_fnames),"does not match the number of output properties in prop:",propndim,".\n")
                        cat("A single vector of descriptors types, e.g. c(\"standard\",\"extended\"), can be provided or a list of vectors of descriptors types of length that matches the number of output properties in prop.\n\n")
                        return(NULL)
                      } else {
                        fnames <<- v_fnames
                      }
                    }
                  } else if(v_fnames_class=="character"){
                    for(mi in 1:propndim){
                      fnames[[mi]] <<- v_fnames
                    }
                    cat(v_fnames,"fingerprint/physical descriptors will be used as features for all the properties predictions.\n\n")
                  }else {
                    cat("Caution: A single vector of descriptors types, e.g. c(\"standard\",\"extended\"), can be provided or a list of vectors of descriptors types of length that matches the number of output properties in prop.\n\n")
                    return(NULL)
                  }
                }
                cat(length(unique(unlist(fnames))),"unique type(s) of descriptors is(are) used:\n")
                cat(unique(unlist(fnames)),"\n\n")

                if(length(v_scale)==1){
                  cat("A features scaling v_scale puts to",v_scale,"will be used for all the properties predictions.\n\n")
                  scale <<- rep(v_scale,propndim)
                } else {
                  if(length(v_scale)!=propndim){
                    cat("Please submit a logical vector v_scale of length that equals the number of output properties in prop:",propndim,".\n")
                    cat("A single logical v_scale value, e.g. TRUE, can also be provided.\n\n")
                    return(NULL)
                  } else {
                    scale <<- v_scale
                    cat("A scaling of the features for output properties predictions will be done or not following the logical vector v_scale:\n")
                    cat(scale,"\n\n")
                  }
                }

                if(is.null(v_filterfunc)){
                  filterfunc <<- function()
                  filtermin <<- numeric(0)
                  filtermax <<- numeric(0)
                  cat("No function for filtering has been submitted.\n\n")
                } else {
                  filterfunc_class <- class(v_filterfunc)
                  if(filterfunc_class=="function"){
                    filterfunc <<- v_filterfunc
                    if(class(v_filtermin)=="numeric"&&class(v_filtermax)=="numeric"){
                      filtermin <<- v_filtermin
                      filtermax <<- v_filtermax
                    } else {
                      cat("Caution: v_filtermin and v_filtermax of classes:",class(v_filtermin),"and",class(v_filtermax),"should be both of class numeric.\n\n")
                      return(NULL)
                    }
                  } else {
                    cat("Caution: v_filterfunc of class:",filterfunc_class,"has been submitted.\n")
                    cat("Please, submit a function v_func of class function.\n\n")
                    return(NULL)
                  }
                }

                if(is.null(v_func)){
                  cat("No function for empirical calculations has been submitted.\n\n")
                  func <<- list()
                } else {
                  func_class <- class(v_func)
                  if(func_class=="function"){
                    func <<- list(v_func)
                  } else if(func_class=="list"){
                    func_list_class <- unique(unlist(lapply(v_func,class)))
                    if(length(func_list_class)!=1){
                      cat("Caution: a list of functions v_func that contains more than one class:",func_list_class,"has been submitted.\n")
                      cat("Please, submit a new function or list of functions v_func containing only objects from the class function.\n\n")
                      return(NULL)
                    } else {
                      if(func_list_class!="function"){
                        cat("Caution: a list of functions v_func that doesn't contain function objects has been submitted.\n")
                        cat("Please, submit a new function or list of functions v_func containing only objects from the class function.\n\n")
                        return(NULL)
                      } else {
                        func <<- v_func
                      }
                    }
                  }
                  cat("A list of",length(func),"function(s) for empirical calculations has been submitted.\n\n")
                }

                if(is.null(v_func_args)){
                  if(!is.null(v_func)){
                    cat("Caution: a function or list of functions v_func for empirical calculations has been submitted without arguments location v_func_args in prop.\n")
                    cat("A vector, or list of vectors, of integers as arguments location v_func_args in the output properties prop (columns location in prop) should be submitted,\n")
                    cat("e.g. for a single function v_func using the output properties in the second and third columns of prop as arguments, v_func_args=c(2,3),\n")
                    cat("e.g. for a list of two functions v_func using the output properties in prop, v_func_args=list(c(2,3),c(1)),\n")
                    cat("where the first function v_func[[1]] uses the second and third, and the second function v_func[[2]] uses only the first, output properties in prop.\n\n")
                    return(NULL)
                  } else {
                    func_args <<- list()
                  }
                } else {
                  func_args_class <- class(v_func_args)
                  if(func_args_class=="numeric"){
                    if(length(v_func)!=1){
                      cat("Caution: a single vector of integers v_func_args is provided for a list of",length(v_func),"functions v_func for empirical calculations.\n")
                      cat("Please, equalize the v_func and v_func_args lists length,\n")
                      cat("e.g. for a list of two functions v_func using the output properties in prop, v_func_args=list(c(2,3),c(1)),\n")
                      cat("where the first function v_func[[1]] uses the second and third, and the second function v_func[[2]] uses only the first, output properties in prop.\n\n")
                      return(NULL)
                    } else {
                      func_args <<- list(v_func_args)
                      cat("A single vector of integers v_func_args has been submitted as arguments location for a single function v_func for empirical calculations.\n\n")
                    }
                  } else if(func_args_class=="list"){
                    if(length(v_func_args)!=length(v_func)){
                      cat("Caution: a list of",length(v_func_args),"vectors of integers v_func_args is provided for a list of",length(v_func),"functions v_func for empirical calculations.\n")
                      cat("Please, equalize the v_func and v_func_args lists length,\n")
                      cat("e.g. for a list of two functions v_func using the output properties in prop, v_func_args=list(c(2,3),c(1)),\n")
                      cat("where the first function v_func[[1]] uses the second and third, and the second function v_func[[2]] uses only the first, output properties in prop.\n\n")
                      return(NULL)
                    } else {
                      func_args_list_class <- unique(unlist(lapply(v_func_args,class)))
                      if(length(func_args_list_class)!=1){
                        cat("Caution: a list of",func_args_list_class,"class objects has been submitted in v_func_args, when it should only contain numerical vectors.\n")
                        cat("A vector, or list of vectors, of integers as arguments location v_func_args in the output properties prop (columns location in prop) should be submitted,\n")
                        cat("e.g. for a list of two functions v_func using the output properties in prop, v_func_args=list(c(2,3),c(1)),\n")
                        cat("where the first function v_func[[1]] uses the second and third, and the second function v_func[[2]] uses only the first, output properties in prop.\n\n")
                        return(NULL)
                      } else {
                        if(func_args_list_class=="numeric"){
                          func_args <<- v_func_args
                          cat("A list of",length(v_func_args),"vectors of integers has been submitted as arguments location for a list of",length(v_func),"functions v_func for empirical calculations.\n\n")
                        } else {
                          cat("Caution: the list of submitted vectors in v_func_args aren't numerical.\n")
                          cat("Please check it and/or cast them via as.numeric().\n\n")
                          return(NULL)
                        }
                      }
                    }
                  }
                }

                cat("Check of the SMILES strings in progress...\n")
                for(sli in 1:length(smis)){
                  mols <- parse.smiles(smis[[sli]], kekulise=kekulise)
                  mols.na <- which(is.na(mols)==T)
                  lmols.na <- length(mols.na)
                  if(lmols.na!=0){
                    prop[[sli]] <- as.matrix(prop[[sli]][-mols.na])
                    write(smis[[sli]][mols.na], paste("incorrectSMILES_",sli,".smi",sep = ""), append=F, sep="\n")
                    smis[[sli]] <- smis[[sli]][-mols.na]
                    cat(lmols.na,"input SMILES didn't pass the parse.smiles() function of rcdk package and aren't kept in the process.\n")
                    cat("Incorrect SMILES are listed in the ",getwd(),"/incorrectSMILES_",sli,".smi file.\n\n", sep = "")
                  }
                }
                cat("Check of the SMILES strings done.\n\n")

                cat("Computing the descriptors for the submitted SMILES strings...\n")
                X <<- list()
                Y <<- list()
                for(sli in 1:length(smis)){
                  scale_init <- ifelse(scale[sli], T, F)
                  Xlist <- get_descriptor(smis = smis[[sli]], desctypes = fnames[[sli]], scale = scale[sli], scale_init = scale_init)
                  X[[sli]] <<- Xlist[[1]] # features
                  Y[[sli]] <<- prop[[sli]] # properties
                  mdesc[[sli]] <<- as.numeric(Xlist[[2]])
                  sddesc[[sli]] <<- as.numeric(Xlist[[3]])
                }
                cat("Computation of the descriptors for the submitted SMILES strings done.\n\n")
              },

              ################################################
              ####  Regression models training function  #####
              ################################################
              model_training = function(model = "linear_Bayes",
                                        params = NA,
                                        n_boot = 10,
                                        s_boot = 0.85,
                                        r_boot = F,
                                        parallelize = F){
                "allows to train regression models, define their parameters, request bootstrap approach and CPU parallelization"
                Models_list = c("linear_Bayes",
                                "elasticnet",
                                "svm",
                                "ranger",
                                "bagging",
                                "gradientboost",
                                "deeplearning")

                model_length <- length(model)

                if(model_length==1){
                  cat(model,"regression model will be trained and used for predictions for the whole set of properties.\n\n")
                  model <- rep(model,propndim)
                  model_length <- length(model)
                } else {
                  if(model_length!=propndim){
                    cat("Please submit a character vector model of length that equals the output properties number in prop:",propndim,",\n")
                    cat("or a single character object is sufficient, e.g. model = \"linear_Bayes\"\n\n")
                    return(NULL)
                  }
                }

                model_check <- model%in%Models_list
                model_tag <- which(model_check==F)
                if(length(model_tag)!=0){
                  cat("Caution:",model[model_tag],"model(s) is(are) not implemented.\n")
                  cat("Please check the model(s).\n")
                  cat("Type get_Models() for an exhaustive list of implemented models.\n\n")
                  return(NULL)
                }

                if(class(params)=="logical"){
                  params <- as.list(rep(NA,model_length))
                  cat("A list of parameters params hasn't been assigned, then the regression model(s) default parameters are assigned.\n")
                  cat("The default parameters are available from get_Model_params(), e.g. get_Model_params(\"elasticnet\")\n\n")
                } else {
                  if(class(params)=="list"){
                    params_class <- unique(unlist(lapply(params,class)))
                    if(!"list"%in%params_class){
                      if(model_length!=1){
                        model_uniqueness <- length(unique(model))
                        if(model_uniqueness==1){
                          params_tmp <- params
                          params_length <- length(params)
                          params <- list()
                          for(mi in 1:model_length){
                            params[[mi]] <- params_tmp
                          }
                          cat("The same list of",params_length,"parameter(s) will be used for each of the",unique(model),"regression model.\n\n")
                        } else {
                          cat("Caution: a list of parameters has been provided for one model only,\n")
                          cat("but",model_length,"regression models are requested.\n")
                          cat("Please, submit a list of",model_length,"parameters lists respectively assigned to each model.\n\n")
                          return(NULL)
                        }
                      } else {
                        params <- list(params)
                        params_length <- length(params)
                        cat(" A list of",params_length,"parameters has been provided for a",model,"regression model.\n\n")
                      }
                    } else {
                      if(length(params_class)!=1){
                        cat("Caution: a list of different class objects has been provided as params.\n")
                        cat("Please, submit a list of",model_length,"parameters lists respectively assigned to each model.\n\n")
                        return(NULL)
                      } else {
                        params_length <- length(params)
                        if(params_length!=model_length){
                          cat("Caution: a list of parameters lists has been provided for",params_length,"model(s),\n")
                          cat("but",model_length,"regression models are requested.\n")
                          cat("Please, submit a list of",model_length,"parameters lists respectively assigned to each model.\n\n")
                          return(NULL)
                        } else {
                          cat(" A list of",params_length,"parameters lists has been provided for",model_length,"regression model(s):\n")
                          cat(model,"\n\n")
                        }
                      }
                    }
                  } else {
                    cat("A list of parameters or a list of lists of parameters should be provided in params, but params is of class",class(params),"\n\n")
                    return(NULL)
                  }
                }

                if(n_boot <= 0){
                  cat("Caution: The number of requested datasets in the bootstrap approach should be > 0.\n\n")
                  return(NULL)
                } else {
                  cat(n_boot,"datasets are requested for a bootstrap analysis.\n\n")
                }

                if(s_boot <= 0 || s_boot > 1){
                  cat("Caution: The proportion of input data, used for training, should be defined in ]0,1].\n\n")
                  return(NULL)
                } else {
                  cat(s_boot*100,"% of the whole data will be used for training.\n")
                  replace_b <- ifelse(r_boot,"with","without")
                  cat("Note that the training dataset is re-sampled",replace_b,"replacement from the whole data for each new iteration of a bootstrap analysis.\n\n")
                }

                paravailable <- F
                no_cores <- 0
                cl <- NULL
                if(parallelize){
                  paravailable <- ifelse("doParallel"%in%rownames(installed.packages()), T, F)
                  if(paravailable){
                    require(package = doParallel,quietly = T)
                    no_cores <- detectCores() - 1
                    if(no_cores>0){
                      cat("Parallelization requested. Initialization of a cluster of",no_cores,"cores...\n")
                      cat("(Note that the total of available cores - 1 will be used)\n")
                      cl <- makeCluster(no_cores)
                      registerDoParallel(cl)
                      cat("Cluster registered and ready.\n")
                      cat("The training will use the multi-core capabilities of the machine.\n\n")
                    } else {
                      cat("The number of available cores (",no_cores+1,") does not allow a parallelization.\n",sep="")
                      cat("The training will be processed sequentially.\n\n")
                    }
                  } else {
                    cat("Parallelization is unavailable because the package doParallel is missing.\n")
                    cat("Please, consider to install the package doParallel via install.packages('doParallel').\n")
                    cat("The training will be processed sequentially.\n\n")
                  }
                } else {
                  cat("Parallelization isn't requested. The training will be processed sequentially.\n\n")
                }

                mxnetavailable <- "mxnet"%in%rownames(installed.packages())
                if("deeplearning"%in%model&&!mxnetavailable){
                  cat("Caution: You have requested deeplearning as a regression model, but the package mxnet is not installed.\n")
                  cat("Please, go to http://mxnet.io/get_started/setup.html for installing mxnet for your system.\n\n")
                  return(NULL)
                }

                trmodel <<- model
                trnboot <<- as.numeric(n_boot)
                trndf <<- as.numeric(unlist(lapply(X,function(x)dim(x)[1])))

                for(sli in 1:propndim){

                  tmp_model <- model[sli]

                  cat("Processing the training of the regression model ",tmp_model," #",sli," of the models list...\n",sep="")
                  tic()
                  if(tmp_model == "linear_Bayes"){
                    linearBayes_model <- get_linearBayes(X[[sli]], Y[[sli]], params[[sli]])
                    save(linearBayes_model, file=paste(getwd(),"/",tmp_model,"_params_",sli,".Rda",sep=""))
                  } else {
                    nrX <- dim(X[[sli]])[1]
                    if(tmp_model == "deeplearning")parallelize <- F
                    if(parallelize&&paravailable){
                      packages_set <- c('glmnet','kernlab','xgboost','ranger')
                      #if(mxnetavailable)packages_set <- c(packages_set,'mxnet')
                      Boot_model_prop <- foreach(nboot=1:n_boot, .combine = rbind, .packages = packages_set, .export = c('X','Y','get_Model'))%dopar%{
                        sboot <- sample(1:nrX, floor(s_boot*nrX),replace = r_boot)
                        yboot = as.matrix(Y[[sli]][sboot])
                        xboot = as.matrix(X[[sli]][sboot,])
                        return(list(get_Model(X = xboot, Y = yboot, model_name = tmp_model, params = params[[sli]])))
                      }
                    } else {
                      Boot_model_prop <- matrix(rep(list(), n_boot),nrow = n_boot, ncol = 1)
                      for(nboot in 1:n_boot){
                        sboot <- sample(1:nrX, floor(s_boot*nrX),replace = r_boot)
                        yboot = as.matrix(Y[[sli]][sboot])
                        xboot = as.matrix(X[[sli]][sboot,])
                        Boot_model_prop[[nboot]] <- get_Model(X = xboot, Y = yboot, model_name = tmp_model, params = params[[sli]])
                      }
                    }
                    # save the n_boot*propndim matrix of regression tmp_models (rgr)
                    # each cell containing a related c(list(rgr),list(rmse),list(mae),list(max_ae),list(corr))
                    # see in AIModels.R for details
                    save(Boot_model_prop, file=paste(getwd(),"/",tmp_model,"_params_",sli,".Rda",sep=""))
                  }
                  cat("Regression model ",tmp_model," #",sli," saved to:\n",getwd(),"/",tmp_model,"_params_",sli,".Rda\n",sep="")
                  exec_time <- toc(quiet=T)
                  exec_time <- exec_time$toc - exec_time$tic
                  cat("Training of the regression model ",tmp_model," #",sli," done in ",exec_time," s.\n",sep="")
                }
                cat("The full training is processed.\n\n")
                if(parallelize&paravailable){
                  stopCluster(cl)
                  cat("Cluster of",no_cores,"cores has been properly stopped.\n\n")
                }
                ##
              },

              ################################################
              ####        Get the full features list      ####
              ################################################
              get_features = function(){
                "returns a list of nxd matrix X with d features for n input SMILES"
                return(X)
              },

              ################################################
              ###        Get the full properties list      ###
              ################################################
              get_props = function(){
                "returns a list of nxp matrix Y of p properties for n input SMILES"
                return(Y)
              },

              ################################################
              ####    Set the targeted properties space   ####
              ################################################
              set_target = function(v_propmin, v_propmax){
                "sets the targeted properties space in vectors propmin and propmax"
                if(class(v_propmin)!="numeric"&class(v_propmax)!="numeric"){
                  cat("Caution: the boundaries vectors v_propmin and v_propmax aren't numerical vectors.\n\n")
                  return(NULL)
                }

                if(length(v_propmin)!=length(v_propmax)){
                  cat("Caution: the boundaries vectors v_propmin and v_propmax have different length.\n\n")
                  return(NULL)
                }

                if(!is.null(func)){
                  if(length(v_propmin)!=(propndim+length(func))){
                    cat("Empirical calculations via a list of",length(func),"functions func are expected, but the length of the boundaries\n")
                    cat("vectors v_propmin and v_propmax is different. Please, check it.\n")
                    cat("It is exected that for a newly empirically calculated output property, minimum and maximum target boundaries are defined,")
                    cat("i.e. if 2 output properties exists in prop and 1 supplementary property is calculated in fine, v_propmin and v_propmax should have a length of 3.\n")
                    cat("For example, v_propmin=c(2,12,0.9) or v_propmin=c(NA,NA,0.9) if the boundaries for the two first output properties are unknown or set free.\n\n")
                    return(NULL)
                  } else {
                    propmin <<- v_propmin
                    propmax <<- v_propmax
                    cat("Targeted properties space boundaries:\n")
                    cat("Min:",propmin,"\n")
                    cat("Max:",propmax,"\n")
                  }
                } else {
                  if(length(v_propmin)!=propndim){
                    cat("Caution: The length of the boundaries numerical vectors v_propmin and v_propmax, ",length(v_propmin),", should match the number of output properties in prop, ",length(Y),".\n\n",sep="")
                    return(NULL)
                  } else {
                    propmin <<- v_propmin
                    propmax <<- v_propmax
                    cat("Targeted properties space boundaries:\n")
                    cat("Min:",propmin,"\n")
                    cat("Max:",propmax,"\n")
                  }
                }
              },

              ################################################
              ####     Get the predictions from SMILES    ####
              ################################################
              qspr_predict = function(smis=NULL){
                "predicts properties for input SMILES from a given regression model"
                predy <- c()
                predvar <- c()
                for(sli in 1:length(X)){
                  newx <- get_descriptor(smis = smis, desctypes = fnames[[sli]],
                                         scale = scale[sli], scale_init = FALSE,
                                         mdesc = mdesc[[sli]], sddesc = sddesc[[sli]])[[1]]

                  tmp_model <- trmodel[sli]

                  filenametmp <- paste(getwd(),"/",tmp_model,"_params_",sli,".Rda",sep="")
                  if(file.exists(filenametmp)){
                    load(filenametmp)
                  } else {
                    cat("The file containing the trained forward model does not exist.\n")
                    cat("Check your working directory via getwd(), please.\n")
                    cat("The working directory should be the same than the directory containing the trained regression models.\n")
                    cat("If the trained regression models aren't yet created, train the forward model via model_training() and try again, please.\n\n")
                    return(NULL)
                  }

                  if(tmp_model == "linear_Bayes"){
                    pred_tmp <- use_linearBayes(newx = newx, model = linearBayes_model)
                    predy <- rbind(predy, pred_tmp[[1]])
                    predvar <- rbind(predvar, pred_tmp[[2]])
                  } else {
                    predy_boot <- c()
                    for(nboot in 1:trnboot){
                      pred_tmp <- use_Model(newx = newx, model_name = tmp_model, model = Boot_model_prop[[nboot]])
                      predy_boot <- rbind(predy_boot, pred_tmp)
                    }
                    predy_tmp <- t(apply(predy_boot,2,mean, na.rm =TRUE))
                    predvar_tmp <- t(apply(predy_boot,2,var, na.rm =TRUE))
                    rm(predy_boot)
                    predy <- rbind(predy, predy_tmp)
                    predvar <- rbind(predvar, predvar_tmp)
                  }
                }
                rm(newx)
                rm(pred_tmp)

                if(length(func)!=0&length(func_args)!=0){
                  for(fli in 1:length(func)){
                    func_argslist <- list()
                    func_argslist[[1]] <- as.matrix(predy[func_args[[fli]],])
                    func_argslist[[2]] <- as.matrix(predvar[func_args[[fli]],])
                    analytic_pred <- do.call(func[[fli]],func_argslist)

                    predy <- rbind(predy,t(analytic_pred[[1]]))
                    predvar <- rbind(predvar,t(analytic_pred[[2]]))
                  }
                  rm(func_argslist)
                  rm(analytic_pred)
                }

                return(list(predy,predvar))
              },

              ################################################
              ####     Get the predictions from SMILES    ####
              ####   and likelihood to reach the target   ####
              ################################################
              iqspr_predict = function(smis=NULL, temp=c(1,1)){
                "predicts properties for input SMILES from a given regression model and evaluates the probability to
                reach a targeted properties space"
                forward_pred <- qspr_predict(smis)
                predy <- forward_pred[[1]]
                if(length(temp)==1){
                  predvar <- (temp^2)*forward_pred[[2]]
                } else {
                  predvar <- diag(temp)^2 %*% forward_pred[[2]]
                }
                rm(forward_pred)

                res <- rep(1, length(smis))
                for(i in 1:length(propmin)){
                  if(!is.na(propmin[i])&&!is.na(propmax[i])){
                    lx <- (propmin[i]-predy[i,])/sqrt(predvar[i,])
                    ux <- (propmax[i]-predy[i,])/sqrt(predvar[i,])
                    res_tmp <- (sapply(ux, function(x) pt(x, trndf[i])) - sapply(lx, function(x) pt(x, trndf[i])))
                  } else {
                    res_tmp <- 1
                  }
                  res <- res * res_tmp
                }

                if(length(filtermin)!=0&&length(filtermax)!=0){
                  filterval <- filterfunc(smis)
                  filter_veto <- t(t(filterval)>filtermin&t(filterval)<filtermax)
                  for(i in 1:dim(filter_veto)[2]){
                    res <- res * as.numeric(filter_veto[,i])
                  }
                  if(sum(res)==0)res <- rep(1, length(smis))
                }

                return(list(res,predy,predvar))
              }
  )
)
