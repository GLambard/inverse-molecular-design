## ------------------------------------------------------------------------
library(iqspr)
library(ggplot2)
library(gridExtra) # for multiple plots on a page

## ------------------------------------------------------------------------
data("qspr.data")
dim(qspr.data)
head(qspr.data)

## ------------------------------------------------------------------------
smis <- as.character(qspr.data[,1])
prop <- as.matrix(qspr.data[,c(2,3)])

## ----fig.width=4, fig.height=4, fig.align='center'-----------------------
trainidx <- sample(1:nrow(qspr.data), 1000)
testidx <- sample((1:nrow(qspr.data))[-trainidx], 500)
dt <- data.frame(prop[trainidx,])
colnames(dt) <- c("E","HOMOLUMOgap")
ggplot(data = dt, aes(x = E, y = HOMOLUMOgap)) + geom_point(size=0.4, color="black") +
  labs(x="E", y="HOMO-LUMO gap", title="Initial dataset") + ylim(c(0,max(prop[,2]))) + xlim(c(0,max(prop[,1]))) + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))

## ------------------------------------------------------------------------
qsprpred_env <- QSPRpred()
qsprpred_env$init_env(smis=smis[trainidx], prop=prop[trainidx,], v_fnames=c("standard","extended","circular"))

## ------------------------------------------------------------------------
get_descriptors()

## ------------------------------------------------------------------------
features <- qsprpred_env$get_features()
length(features)
lapply(features,dim)

## ------------------------------------------------------------------------
properties <- qsprpred_env$get_props()
length(properties)
lapply(properties,dim)

## ---- eval=FALSE, include=TRUE-------------------------------------------
#  qsprpred_env$model_training(model=c("linear_Bayes","ranger"),params=list(list(NA),list("num.trees" = 200)),n_boot=10,s_boot = 0.85,r_boot = F,parallelize=T)

## ------------------------------------------------------------------------
get_Models()

## ------------------------------------------------------------------------
get_Model_params("elasticnet")

## ------------------------------------------------------------------------
get_Model_params("ranger")

## ---- eval=FALSE, include=TRUE-------------------------------------------
#  predictions <- qsprpred_env$qspr_predict(smis[testidx])

## ------------------------------------------------------------------------
data("predictions")

## ------------------------------------------------------------------------
# Example for 4 compounds
cat("Predictions:\n")
rownames(predictions[[1]]) <- c("E","HOMO-LUMO gap")
predictions[[1]][,1:4] # predictions
cat("\nVariances:\n")
rownames(predictions[[2]]) <- c("E","HOMO-LUMO gap")
predictions[[2]][,1:4] # variances

## ---- fig.width=8, fig.height=4, fig.align='center'----------------------
d1 <- data.frame(predictions[[1]][1,],prop[testidx,"E"],sqrt(predictions[[2]][1,]))
colnames(d1) <- c("pred","prop","sd")
d2 <- data.frame(predictions[[1]][2,],prop[testidx,"HOMO-LUMO gap"],sqrt(predictions[[2]][2,]))
colnames(d2) <- c("pred","prop","sd")

prd.obs.rmse1 <- round(sqrt(sum((predictions[[1]][1,] - prop[testidx,"E"])^2) / length(testidx)),digits = 1)
prd.obs.corr1 <- round(cor(predictions[[1]][1,],prop[testidx,"E"]),digits = 2)
prd.obs.rmse2 <- round(sqrt(sum((predictions[[1]][2,] - prop[testidx,"HOMO-LUMO gap"])^2) / length(testidx)),digits = 1)
prd.obs.corr2 <- round(cor(predictions[[1]][2,],prop[testidx,"HOMO-LUMO gap"]),digits = 2)

minmaxx <- c(min(d1[,1]),max(d1[,1]))
minmaxy <- c(min(d1[,2]),max(d1[,2]))
p1 <- ggplot(data = d1, aes(x = pred, y = prop, size=sd)) + geom_point(color="cyan3") + geom_point(shape=1, alpha=0.4) +
  labs(x="predictions", y="observations", title="E") + ylim(minmaxy) + xlim(minmaxx) + 
  guides(size=guide_legend(title="s.d.")) + 
  geom_abline(intercept = 0, slope = 1, color = "red", size = 1, alpha = 0.5) +
  annotate("text", x = c(100,100), y = c(max(d1[,"prop"])-100,max(d1[,"prop"])-70), label = c(paste("RMSE:",prd.obs.rmse1), paste("COR:",prd.obs.corr1)) , color="black", size=3, angle=0, fontface="bold") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))

minmaxx <- c(min(d2[,1]),max(d2[,1]))
minmaxy <- c(min(d2[,2]),max(d2[,2]))
p2 <- ggplot(data = d2, aes(x = pred, y = prop, size = sd)) + geom_point(color="deepskyblue") + geom_point(shape=1, alpha=0.4) +
  labs(x="predictions", y="observations", title="HOMO-LUMO gap") + ylim(minmaxy) + xlim(minmaxx) + 
  guides(size=guide_legend(title="s.d.")) + 
  geom_abline(intercept = 0, slope = 1, color = "red", size = 1, alpha = 0.5) +
  annotate("text", x = c(5,5), y = c(max(d2[,"prop"])-1.5,max(d2[,"prop"])-1), label = c(paste("RMSE:",prd.obs.rmse2), paste("COR:",prd.obs.corr2)) , color="black", size=3, angle=0, fontface="bold") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))

grid.arrange(p1,p2,ncol=2)

## ---- eval=FALSE, include=TRUE-------------------------------------------
#  qsprpred_env$model_training(model=c("linear_Bayes"),params=NA) # linear_Bayes prevails here for quick execution
#  targ.min <- c(50,2)
#  targ.max <- c(250,4)
#  qsprpred_env$set_target(targ.min,targ.max)

## ---- eval=FALSE, include=TRUE-------------------------------------------
#  data("trainedSMI")
#  engram_5k <- ENgram$new(trainedSMI, order=10)

## ---- eval=FALSE, include=TRUE-------------------------------------------
#  data("engram_5k")
#  smchem <- SmcChem$new(smis = rep("c1ccccc1O", 25), v_qsprpred = qsprpred_env, v_engram = engram_5k, v_temp=c(30,3))

## ---- eval=FALSE, include=TRUE-------------------------------------------
#  for(i in 1:200){
#    smchem$smcexec(niter = 1, nsteps = 5, preorder = 0.2)
#  }

## ---- eval=FALSE, include=TRUE-------------------------------------------
#  gensmis <- smchem$get_hiscores(nsmi=2000, exsim=0.9)

## ------------------------------------------------------------------------
data("gensmis")
head(gensmis)

## ---- eval=FALSE, include=TRUE-------------------------------------------
#  pred <- qsprpred_env$qspr_predict(gensmis[,1])
#  predmat <- t(pred[[1]])
#  
#  dpred <- data.frame(predmat)
#  colnames(dpred) <- c("E","HOMOLUMOgap")
#  predinit <- (t(qsprpred_env$qspr_predict("c1ccccc1O")[[1]]))
#  colnames(predinit) <- c("E","HOMOLUMOgap")

## ---- fig.width=4, fig.height=4, fig.align='center'----------------------
data("dpred")
data("predinit")

targ.min <- c(50,2)
targ.max <- c(250,4)
ggplot(data = dt, aes(x = E, y = HOMOLUMOgap)) + geom_point(size=0.4, color="black") + 
  annotate("rect", xmin=targ.min[1], xmax=targ.max[1], ymin=targ.min[2], ymax=targ.max[2], alpha=0.2, color="blue", fill="blue") +
  geom_point(data = dpred, aes(x = E, y = HOMOLUMOgap), size=0.4, color="red") +
  geom_point(data = data.frame(predinit), aes(x = E, y = HOMOLUMOgap), size=4, color="green", shape=3) +
  labs(x="E", y="HOMO-LUMO gap", title="Initial dataset") + ylim(c(0,max(prop[,2]))) + xlim(c(0,max(prop[,1]))) + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))

## ---- fig.width=8, fig.height=16, fig.align='center'---------------------
viewstr(gensmis[1:10,1], nrow = 5, ncol = 2, legend = paste("ID:",c(1:10)))

## ---- eval=FALSE, include=TRUE-------------------------------------------
#  cut_l <- 0.9*length(smis)
#  cut <- c(1:cut_l)
#  smis1 <- smis[cut]
#  smis2 <- smis[-cut]
#  prop1 <- prop[cut,"E"]
#  prop2 <- prop[-cut,"HOMO-LUMO gap"]

## ---- eval=FALSE, include=TRUE-------------------------------------------
#  smis1_l <- length(smis1)
#  smis2_l <- length(smis2)
#  trainidx1 <- sample(1:smis1_l, 0.9*smis1_l)
#  trainidx2 <- sample(1:smis2_l, 0.9*smis2_l)
#  qsprpred_env <- QSPRpred()
#  qsprpred_env$init_env(smis=list(smis1[trainidx1],smis2[trainidx2]), prop=list(prop1[trainidx1],prop2[trainidx2]), v_fnames=c("standard","extended","circular"))

## ---- eval=FALSE, include=TRUE-------------------------------------------
#  v_fnames=list(c("standard","extended"),c("graph"))

## ---- eval=FALSE, include=TRUE-------------------------------------------
#  qsprpred_env$model_training(model=c("linear_Bayes","ranger"),params=list(list(NA),list("num.trees" = 200)),n_boot=10,s_boot = 0.85,r_boot = F,parallelize=T)

## ---- eval=FALSE, include=TRUE-------------------------------------------
#  predictions1 <- qsprpred_env$qspr_predict(smis1[-trainidx1])
#  predictions2 <- qsprpred_env$qspr_predict(smis2[-trainidx2])

## ---- eval=FALSE, include=TRUE-------------------------------------------
#  # Example for 4 compounds
#  cat("Predictions:\n")
#  rownames(predictions1[[1]]) <- c("E","HOMO-LUMO gap")
#  predictions1[[1]][,1:4] # predictions
#  cat("\nVariances:\n")
#  rownames(predictions1[[2]]) <- c("E","HOMO-LUMO gap")
#  predictions1[[2]][,1:4] # variances

## ----fig.width=8, fig.height=4, fig.align='center', eval=FALSE, include=TRUE----
#  d1 <- data.frame(predictions1[[1]][1,],prop1[-trainidx1],sqrt(predictions1[[2]][1,]))
#  colnames(d1) <- c("pred","prop","sd")
#  d2 <- data.frame(predictions2[[1]][2,],prop2[-trainidx2],sqrt(predictions2[[2]][2,]))
#  colnames(d2) <- c("pred","prop","sd")
#  
#  prd.obs.rmse1 <- round(sqrt(sum((predictions1[[1]][1,] - prop1[-trainidx1])^2) / (0.1*smis1_l)),digits = 1)
#  prd.obs.corr1 <- round(cor(predictions1[[1]][1,],prop1[-trainidx1]),digits = 2)
#  prd.obs.rmse2 <- round(sqrt(sum((predictions2[[1]][2,] - prop2[-trainidx2])^2) / (0.1*smis2_l)),digits = 1)
#  prd.obs.corr2 <- round(cor(predictions2[[1]][2,],prop2[-trainidx2]),digits = 2)
#  
#  minmaxx <- c(min(d1[,1]),max(d1[,1]))
#  minmaxy <- c(min(d1[,2]),max(d1[,2]))
#  p1 <- ggplot(data = d1, aes(x = pred, y = prop, size=sd)) + geom_point(color="cyan3") + geom_point(shape=1, alpha=0.4) +
#    labs(x="predictions", y="observations", title="E") + ylim(minmaxy) + xlim(minmaxx) +
#    guides(size=guide_legend(title="s.d.")) +
#    geom_abline(intercept = 0, slope = 1, color = "red", size = 1, alpha = 0.5) +
#    annotate("text", x = c(50,50), y = c(500,530), label = c(paste("RMSE:",prd.obs.rmse1), paste("COR:",prd.obs.corr1)) , color="black", size=3, angle=0, fontface="bold") +
#    theme(plot.title = element_text(hjust = 0.5))
#  
#  minmaxx <- c(min(d2[,1]),max(d2[,1]))
#  minmaxy <- c(min(d2[,2]),max(d2[,2]))
#  p2 <- ggplot(data = d2, aes(x = pred, y = prop, size = sd)) + geom_point(color="deepskyblue") + geom_point(shape=1, alpha=0.4) +
#    labs(x="predictions", y="observations", title="HOMO-LUMO gap") + ylim(minmaxy) + xlim(minmaxx) +
#    guides(size=guide_legend(title="s.d.")) +
#    geom_abline(intercept = 0, slope = 1, color = "red", size = 1, alpha = 0.5) +
#    annotate("text", x = c(5,5), y = c(12,13), label = c(paste("RMSE:",prd.obs.rmse2), paste("COR:",prd.obs.corr2)) , color="black", size=3, angle=0, fontface="bold") +
#    theme(plot.title = element_text(hjust = 0.5))
#  
#  grid.arrange(p1,p2,ncol=2)

## ---- eval=FALSE, include=TRUE-------------------------------------------
#  library(ranger)
#  library(rBayesianOptimization)
#  df <- data.frame(properties[[2]],features[[2]]) # 6906 SMILES, 1 property, 3073 features
#  colnames(df) <- c("property",colnames(df[,-1]))
#  id <- sample(dim(df)[1],6000)
#  df.train <- df[id,]
#  df.test <- df[-id,]
#  
#  folds <- cut(seq(1,6000), breaks = 6, labels = FALSE) # 6-folds cross-validation
#  
#  rf_cvfold <- function(mtry, ntree){ # Optimization function to be maximized for Bayesian optimization
#    mae <- c()
#    for(i in 1:6){
#      id <- which(folds == i)
#      valid <- df.train[id,]
#      train <- df.train[-id,]
#      rgr <- ranger(property~., data = train, mtry = mtry, num.trees = ntree)
#      prd <- predict(rgr, data = valid)
#      mae[i] <- mean(abs(valid$property - prd$predictions))
#    }
#    list(Score = -mean(mae), Pred = -mean(mae)) # return a negative mean absolute error (mae)
#  }
#  
#  opt_rf <- BayesianOptimization(rf_cvfold, bounds = list("mtry" = c(10L,1000L), "ntree" = c(100L,500L)), init_points = 10, n_iter = 2)

## ---- eval=FALSE, include=TRUE-------------------------------------------
#  qsprpred_env$model_training(model=c("linear_Bayes","ranger"),params=list(list(NA),list("num.trees" = opt_rf$Best_Par[2], "mtry" = opt_rf$Best_Par[1])),n_boot=100,s_boot = 0.85,r_boot = F,parallelize=T)
#  
#  predictions2 <- qsprpred_env$qspr_predict(smis2[-trainidx2])
#  
#  d2 <- data.frame(predictions2[[1]][2,],prop2[-trainidx2],sqrt(predictions2[[2]][2,]))
#  colnames(d2) <- c("pred","prop","sd")
#  
#  prd.obs.rmse2 <- round(sqrt(sum((predictions2[[1]][2,] - prop2[-trainidx2])^2) / (0.1*smis2_l)),digits = 1)
#  prd.obs.corr2 <- round(cor(predictions2[[1]][2,],prop2[-trainidx2]),digits = 2)
#  
#  minmaxx <- c(min(d2[,1]),max(d2[,1]))
#  minmaxy <- c(min(d2[,2]),max(d2[,2]))
#  ggplot(data = d2, aes(x = pred, y = prop, size = sd)) + geom_point(color="deepskyblue") + geom_point(shape=1, alpha=0.4) +
#    labs(x="predictions", y="observations", title="HOMO-LUMO gap") + ylim(minmaxy) + xlim(minmaxx) +
#    guides(size=guide_legend(title="s.d.")) +
#    geom_abline(intercept = 0, slope = 1, color = "red", size = 1, alpha = 0.5) +
#    annotate("text", x = c(5,5), y = c(12,13), label = c(paste("RMSE:",prd.obs.rmse2), paste("COR:",prd.obs.corr2)) , color="black", size=3, angle=0, fontface="bold") +
#    theme(plot.title = element_text(hjust = 0.5))

## ---- fig.align='center', fig.height=8, fig.width=4, eval= FALSE, include=TRUE----
#  library(mxnet)
#  data <- mx.symbol.Variable("data")
#  fc1 <- mx.symbol.FullyConnected(data, name = "FullyConnected 1", num.hidden = 128)
#  act1 <- mx.symbol.Activation(fc1, name = "Activation 1 tanh", act.type = 'tanh')
#  fc2 <- mx.symbol.FullyConnected(act1, name = "FullyConnected 2", num.hidden = 64)
#  act2 <- mx.symbol.Activation(fc2, name = "Activation 2 tanh", act.type = 'tanh')
#  fc3 <- mx.symbol.FullyConnected(act2, name = "FullyConnected 3", num.hidden = 1)
#  lro <- mx.symbol.LinearRegressionOutput(fc3, name = "Linear Regression Output")
#  graph.viz(lro, graph.title = "Custom graph", graph.width.px = 240, graph.height.px = 640, graph.title.font.name = "arial")

## ---- eval=FALSE, include=TRUE-------------------------------------------
#  nn_params <- list("symbol" = lro,
#                    "ctx" = mx.cpu(), # mx.gpu() recommanded
#                    "num.round" = 5,
#                    "array.batch.size" = 20,
#                    "learning.rate" = 2e-6,
#                    "momentum" = 0.9,
#                    "eval.metric" = mx.metric.rmse,
#                    "array.layout" = "rowmajor")
#  qsprpred_env$model_training(model=c("linear_Bayes","deeplearning"),params=list(list(NA),nn_params),n_boot=1,s_boot = 0.85,r_boot = F,parallelize=T)

## ---- eval=FALSE, include=TRUE-------------------------------------------
#  data <- mx.symbol.Variable("data")
#  fc1 <- mx.symbol.FullyConnected(data, name = "FullyConnected 1", num.hidden = 1)
#  lro <- mx.symbol.LinearRegressionOutput(fc1, name = "Linear Regression Output")
#  graph.viz(lro, graph.title = "Custom graph", graph.width.px = 240, graph.height.px = 640, graph.title.font.name = "arial")

## ---- eval=FALSE, include=TRUE-------------------------------------------
#  pfunc <- function(EVarlist = list()){
#    predy <- EVarlist[[1]]   # matrix with predictions for A and B
#    predvar <- Evarlist[[2]] # matrix with variances over the predictions for A and B
#    EA <- predy[1,]          # expected (predicted) values for A
#    EB <- predy[2,]          # expected (predicted) values for B
#    VarA <- predvar[1,]      # variances for A
#    VarB <- predvar[2,]      # variances for B
#    EP <- EA+EB              # calculated expected values for P
#    VarP <- VarA + VarB      # calculated variances for P, if A and B are uncorrelated
#  
#    return(list(EP,VarP))    # return a list of two vectors of expected values and variances for P
#  }

## ---- eval=FALSE, include=TRUE-------------------------------------------
#  qsprpred_env$init_env(smis=list(smis1[trainidx1],smis2[trainidx2]), prop=list(prop1[trainidx1],prop2[trainidx2]), v_fnames=c("standard"), v_func = pfunc, v_func_args = c(1,2))

## ------------------------------------------------------------------------
v_func_args = c(1,3)

## ---- eval=FALSE, include=TRUE-------------------------------------------
#  targ.min <- c(NA,NA,100)
#  targ.max <- c(NA,NA,150)
#  qsprpred_env$set_target(targ.min,targ.max)

## ---- eval=FALSE, include=TRUE-------------------------------------------
#  mw_filter <- function(smiles){
#    mols <- parse.smiles(smiles, kekulise = F)
#    hidout <- lapply(mols,do.aromaticity)
#    hidout <- lapply(mols,do.isotopes)
#    hidout <- lapply(mols,do.typing)
#    mw <- as.numeric(lapply(mols,get.exact.mass))
#    return(mw)
#  }

## ---- eval=FALSE, include=TRUE-------------------------------------------
#  qsprpred_env <- QSPRpred()
#  qsprpred_env$init_env(smis=smis[trainidx], prop=prop[trainidx,], v_fnames=c("standard"),
#                        v_filterfunc = mw_filter, v_filtermin = c(100), v_filtermax = c(200))

