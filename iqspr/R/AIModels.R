#' get a list of the implemented regression algorithms
#' @description Get a list of the implemented regression algorithms.
#'
#' @examples
#' list_models <- get_Models()
#' @return a list of the implemented regression algorithms with their nickname, as implemented in iqspr, and their common
#' name from literature.
#'
#' @export get_Models
get_Models <- function(){
  Models_list = c("linear_Bayes",
                  "elasticnet",
                  "svm",
                  "ranger",
                  "bagging",
                  "gradientboost",
                  "deeplearning")
  Models_list <- cbind(Models_list, c("Bayesian linear",
                                      "Elastic Net",
                                      "Support Vector Machines",
                                      "Random Forests",
                                      "Bagging",
                                      "Gradient Boosting",
                                      "Deep Learning"))
  Models_list <- cbind(Models_list, c("custom implementation",
                                      "glmnet",
                                      "kernlab",
                                      "ranger",
                                      "ranger",
                                      "xgboost",
                                      "mxnet"))

  Models_list <- cbind(Models_list, c("custom implementation",
                                      "cv.glmnet",
                                      "ksvm",
                                      "ranger",
                                      "ranger",
                                      "xgboost",
                                      "mx.model.FeedForward.create"))

  colnames(Models_list) <- c("nickname","name","package", "call")

  print(Models_list)
}

#' get a list of default parameters for a given regression algorithm
#' @description Get a list of default parameters for a given regression algorithm.
#'
#' @param model_name is the model's name for which you want to display the default parameters ("elasticnet" by default).
#' For the "svm" model, a.k.a Support Vector Machine, Gaussian RBF or Laplace kernel ("rbfdot" or "laplacedot" respectively)
#' are supported (limitation from the \code{\link[kernlab]{sigest}} function used as optimizer on the sigma parameter).
#'
#' @examples
#' list_params <- get_Model_params(model_name = "elasticnet")
#' @return a list of the default parameters for a given regression algorithm.
#'
#' @export get_Model_params
get_Model_params <- function(model_name = "elasticnet"){

  cat("\n")
  cat("Default parameters for the ", model_name, " regression model\n\n")
  if(model_name=="elasticnet"){
    params <- list("alpha" = 0.5)
  } else if(model_name=="svm"){
    params <- list("kernel" = "rbfdot",
                   "type" = "nu-svr",
                   "scaled" = T,
                   "C" = 50,
                   "cross" = 3)
  } else if(model_name=="ranger"){
    params <- list("num.trees" = 200)
  } else if(model_name=="bagging"){
    params <- list("num.trees" = 200)
  } else if(model_name=="gradientboost"){
    params <- list("params"= list("objective" = "reg:linear",
                   "eta" = 0.3,
                   "max_depth" = 6,
                   "min_child_weight" = 1,
                   "subsample" = 1,
                   "colsample_bytree" = 1),
                   "nrounds" = 10,
                   "nfold" = 5)
  } else if(model_name=="deeplearning"){
    params <- list("symbol" = "(see tutorial)",
                   "ctx" = "mx.cpu()",
                   "num.round" = 50,
                   "array.batch.size" = 20,
                   "learning.rate" = 2e-6,
                   "momentum" = 0.9,
                   "eval.metric" = "mx.metric.rmse",
                   "array.layout" = "rowmajor")

    # params <- list("hidden_node" = rep(500,5),
    #                "activation" = "relu",
    #                "dropout" = 0.5,
    #                "out_node" = 1,
    #                "out_activation" = "rmse",
    #                "optimizer" = "sgd",
    #                "num.round" = 50,
    #                "learning.rate" = 0.1,
    #                "eval.metric" = mx.metric.rmse,
    #                "verbose" = F,
    #                "array.layout" = "rowmajor")

  } else if(model_name=="linear_Bayes"){
    params <- list("w0" = "matrix(0, ncol(X), ncol(Y))",
                   "V0_inv" = "diag(c(0,rep(1, ncol(X)-1)))",
                   "b0" = "numeric(ncol(Y))",
                   "a0" = "0")
  } else {
    params <- list("Caution" = "This model is not implemented. Please check the name, or type get_Models() for an exhaustive list of available models.")
  }

  print(params)
}


#' get the computation of parameters issued from the training of a Bayesian linear regression
#' @description Get the computation of parameters issued from the training of a Bayesian linear regression.
#' See the iqspr package paper for details concerning the definition of the cited parameters.
#' @param X is a matrix of features (e.g. fingerprints and/or physical descriptors) issued from a training set.
#' @param Y is a one-dimensional vector, or a matrix, of observables (e.g. targeted chemical properties) issued from a training
#' set.
#' @param params corresponds to the initialization of the parameters (see \code{\link{get_Model_params}} for a parameters
#' list) if known \emph{a priori} (NA by default).
#'
#' @examples
#' \dontrun{linearBayes_model <- get_linearBayes(X = X[training,], Y = Y[training,],
#' params = c(list(w0),list(V0_inv),list(a0),list(b0)))}
#' @return the latest w0, V0_inv, a0 and b0 parameters after training over the data X and Y.
#'
#' @export get_linearBayes
get_linearBayes <- function(X = NULL, Y = NULL, params = NA){
  if(length(params)==4){
    w0=params[["w0"]]
    V0_inv=params[["V0_inv"]]
    a0=params[["a0"]]
    b0=params[["b0"]]
    cat("Bayesian linear regression model is initialized\n")
  } else {
    w0 <- matrix(0, ncol(X), ncol(Y))
    V0_inv <- diag(c(0,rep(1, ncol(X)-1)))  # uninformative prior I (id matrix) with a 0 in [1,1] to avoid
    # regularization of the intercept
    b0 <- numeric(ncol(Y))
    a0 <- 0
    cat("Bayesian linear regression model sets by default\n")
  }

  cat("getting posterior distribution of parameters...\n")
  XX <- t(X)%*%X

  VnMat <- (V0_inv+XX)
  cholStatus <- try(cVnMat <- chol(VnMat), silent = TRUE)
  cholError <- ifelse(class(cholStatus) == "try-error", TRUE, FALSE)
  if(!cholError){
    TVn <- chol2inv(cVnMat) # twice faster than solve(V0_inv+XX)
  }else{
    TVn <- solve(VnMat) # return the inverse of (V0_inv+XX)
  }

  nw0 <- all(abs(w0) < .Machine$double.eps) # check if w0 is a zero matrix
  if(nw0){
    Twn <- TVn%*%(t(X)%*%Y)
  }else{
    Twn <- TVn%*%(V0_inv%*%w0+t(X)%*%Y)
  }

  TaN <- a0 + nrow(X)/2

  if(nw0){
    TbN <- b0 + diag(0.5*( (t(Y)%*%Y) - (t(Twn)%*%(V0_inv+XX)%*%Twn) ))
  }else{
    TbN <- b0 + diag(0.5*((t(w0) %*% V0_inv %*% w0) + (t(Y) %*% Y) - (t(Twn) %*% (V0_inv + XX) %*% Twn) ))
  }

  return(list(Twn,TVn,TaN,TbN))
}

#' construct a given regression model thanks to a training set
#' @description Construct a given regression model thanks to a training set.
#'
#' @param X is a matrix of features (e.g. fingerprints and/or physical descriptors) issued from a training set.
#' @param Y is a one-dimensional vector, or a matrix, of observables (e.g. targeted chemical properties) issued from a training
#' set.
#' @param model_name is the model's name to be built ("elasticnet" by default).
#' @param params is a list of fixed parameters applied to the model (see \code{\link{get_Model_params}} for a detailed view of the
#' default paramaters per model). For the "svm" model, a.k.a Support Vector Machine, Gaussian RBF or Laplace kernel
#' ("rbfdot" or "laplacedot" respectively) are supported (limitation from the \code{\link[kernlab]{sigest}} function used as optimizer on the sigma parameter).
#'
#' @examples
#' \dontrun{elasticNet_model <- get_Model(X = X[training,], Y = Y[training,],
#' model = "elasticnet", params = list("alpha" = 0.5))}
#' @return the desired regression model optimized.
#'
#' @import glmnet
#' @importFrom kernlab ksvm sigest
#' @import xgboost
#' @import ranger
#'
#' @export get_Model
get_Model <- function(X = NULL, Y = NULL, model_name = "elasticnet", params = NA){

  if(unique(is.na(params))){
    if(model_name=="elasticnet"){
      params <- list("alpha" = 0.5)
    } else if(model_name=="svm"){
      params <- list("kernel" = "rbfdot",
                     "type" = "nu-svr",
                     "scaled" = T,
                     "C" = 50,
                     "cross" = 3)
    } else if(model_name=="ranger"){
      params <- list("num.trees" = 200)
    } else if(model_name=="bagging"){
      params <- list("num.trees" = 200)
    } else if(model_name=="gradientboost"){
      params <- list("params"= list("objective" = "reg:linear",
                    "eta" = 0.3,
                    "max_depth" = 6,
                    "min_child_weight" = 1,
                    "subsample" = 1,
                    "colsample_bytree" = 1),
                    "nrounds" = 10,
                    "nfold" = 5)
    } else if(model_name=="deeplearning"){
      data <- mx.symbol.Variable("data")
      fc1 <- mx.symbol.FullyConnected(data, num.hidden = 1)
      act <- mx.symbol.Activation(fc1, act.type = 'tanh')
      fc <- mx.symbol.FullyConnected(act, num.hidden = 1)
      lro <- mx.symbol.LinearRegressionOutput(fc)
      params <- list("symbol" = lro,
                     "ctx" = mx.cpu(),
                     "num.round" = 50,
                     "array.batch.size" = 20,
                     "learning.rate" = 2e-6,
                     "momentum" = 0.9,
                     "eval.metric" = mx.metric.rmse,
                     "array.layout" = "rowmajor")

      #params <- list("hidden_node" = rep(500,5),
      #                  "activation" = "relu",
      #                  "dropout" = 0.5,
      #                  "out_node" = 1,
      #                  "out_activation" = "rmse",
      #                  "optimizer" = "sgd",
      #                  "num.round" = 50,
      #                  "learning.rate" = 0.1,
      #                  "eval.metric" = mx.metric.rmse,
      #                  "verbose" = F,
      #                  "array.layout" = "rowmajor")

    }
  } else {
    params = params
  }

  k <- numeric(0)
  if(model_name=="elasticnet"){
    # elastic-net algorithm
    rgr <- do.call(cv.glmnet, c(list(x = X, y = Y, family = "gaussian"),params))
    ##
  } else if (model_name=="svm"){
    # support vector machine algorithm
    k <- which(apply(X, 2, var) != 0)
    srange <- sigest(x = X[,k])[2]
    if(params$kernel=="rbfdot"||params$kernel=="laplacedot"){
      rgr <- do.call(ksvm, c(list(y = Y, x = X[,k], kpar = list(sigma = srange)),params))
    } else {
      rgr <- do.call(ksvm, c(list(y = Y, x = X[,k]),params))
    }
    ##
  } else if(model_name=="ranger"){
    # random forest algorithm for large amount of features and data (e.g. ranger package)
    data.train <- data.frame(X,Y,row.names = c(1:dim(X)[1]))
    colnames(data.train) <- c(colnames(X),"prop")
    rgr <- do.call(ranger, c(list(formula = prop~.,data=data.train),params))
    ##
  } else if(model_name=="bagging"){
    # random forest algorithm for large amount of features and data (e.g. ranger package)
    data.train <- data.frame(X,Y,row.names = c(1:dim(X)[1]))
    colnames(data.train) <- c(colnames(X),"prop")
    rgr <- do.call(ranger, c(list(formula = prop~.,data=data.train, mtry = ncol(X)),params))
    ##
  } else if(model_name=="gradientboost"){
    # xgboost algorithm
    k <- which(apply(X, 2, var) != 0)
    bst.cv <- do.call(xgb.cv,c(list(data = data.matrix(X[ , k]), label = Y, verbose = FALSE),params))
    nround <- which.min(bst.cv$evaluation_log$test_rmse_mean)
    rgr <- xgboost(params = params, data = data.matrix(X[ , k]), label = Y, nrounds = nround, verbose = FALSE)
    ##
  } else if(model_name=="deeplearning"){
    # mxnet deep learning algorithm
    d0 <- cbind(X,Y)
    colnames(d0) <- c(1:(ncol(d0)-1),"target")
    varnames <- setdiff(colnames(d0),"target")
    #rgr <- do.call(mx.mlp, c(list(data = data.matrix(d0[,varnames]), label = d0[,"target"]), params))
    mx.set.seed(1)
    rgr <- do.call(mx.model.FeedForward.create, c(list(X = data.matrix(d0[,varnames]), y = d0[,"target"]), params))
    ##
  }

  return(list(rgr,params,k))
}

#' get the prediction from a Bayesian linear regression
#' @description Get the prediction from a Bayesian linear regression. For a prediction other than from a Bayesian linear
#' model, see \code{\link{use_Model}}.
#'
#' @param newx is a matrix of features (e.g. fingerprints and/or physical descriptors), issued from a test/validation set, for
#' which the properties (e.g. targeted chemical properties) are predicted.
#' @param model is the set of parameters defining a trained Bayesian linear regression model, issued from \code{\link{get_linearBayes}},
#' to be used as a predictor.
#'
#' @examples
#' \dontrun{linearBayes_pred <- use_linearBayes(newx = newx,
#' model = c(list(w0),list(V0_inv),list(a0),list(b0)))}
#'
#' @return the predicted properties with associated variance.
#'
#' @export use_linearBayes
use_linearBayes <- function(newx = NULL, model = NULL){

  Twn <- model[[1]]
  TVn <- model[[2]]
  TaN <- model[[3]]
  TbN <- model[[4]]

  predy <- t(Twn) %*% t(newx)
  temp1 <- TbN/TaN
  temp2 <- diag((1 + newx %*% TVn %*% t(newx)))
  predvar <- temp1 %*% t(temp2)

  return(list(predy,predvar))
}

#' get the prediction from a given model
#' @description Get the prediction from a desired regression model. For a Bayesian linear prediction, see \code{\link{use_linearBayes}}.
#'
#' @param newx is a matrix of features (e.g. fingerprints and/or physical descriptors), issued from a test/validation set, for
#' which the properties (e.g. targeted chemical properties) are predicted.
#' @param model_name is the model's name to be built ("elasticnet" by default).
#' @param model is the trained regression model, issued from \code{\link{get_Model}}, to be used as a predictor.
#'
#' @examples
#' \dontrun{elasticNet_pred <- use_Model(newx, model_name = "elasticnet", model = elasticnet_model)}
#'
#' @return the predicted properties.
#'
#' @import glmnet
#' @importFrom kernlab ksvm sigest
#' @import xgboost
#' @import ranger
#'
#' @importFrom stats predict
#' @export use_Model
use_Model <- function(newx = NULL, model_name = "elasticnet", model = NULL){

  params <- model[[2]]

  if(model_name=="elasticnet"){
    predy <- t(as.matrix(predict(model[[1]], s = "lambda.min", newx = as.matrix(newx))))
  } else if(model_name=="svm"){
    k <- model[[3]]
    predy <- t(as.matrix(predict(model[[1]], as.matrix(newx[,k]))))
  } else if(model_name=="ranger"){
    data.test <- data.frame(newx,row.names = c(1:dim(newx)[1]))
    colnames(data.test) <- colnames(newx)
    predy <- t(as.matrix(predict(model[[1]], dat = data.test)$predictions))
  } else if(model_name=="bagging"){
    data.test <- data.frame(newx,row.names = c(1:dim(newx)[1]))
    colnames(data.test) <- colnames(newx)
    predy <- t(as.matrix(predict(model[[1]], dat = data.test)$predictions))
  } else if(model_name=="gradientboost"){
    k <- model[[3]]
    predy <- t(as.matrix(predict(model[[1]], data.matrix(newx[,k]))))
  } else if(model_name=="deeplearning"){
    colnames(newx) <- c(1:ncol(newx))
    predy <- as.matrix(predict(model[[1]], data.matrix(newx), array.layout = "rowmajor"))
  }

  return(predy)
}
