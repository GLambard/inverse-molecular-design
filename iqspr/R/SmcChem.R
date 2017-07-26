#'SMILES generator
#' @description SMILES generator thanks to a sequential Monte-Carlo sampler
#'
#' @param smis is an initial vector of SMILES from which the generation of novel SMILES begins.
#' @param v_engram is an Engram object, \emph{a priori} created, encapsulating the SMILES grammar
#' (see \code{\link{ENgram}} for more details).
#' @param v_m is a positive scalar representing the order of a used ENgram model for generation.
#' @param v_qsprpred is a QSPRpred object in which a regression model, initially trained, is accessible for properties
#' predictions (i.e. physico-chemical properties here) of compounds from newly created SMILES.
#' @param v_temp is a vector of numerical values, with a length equals to the number of properties, which represents
#' the annealing parametrization in the sequential Monte-Carlo sampler.
#' @param v_decay is a positive scalar corresponding to the decay rate of temp above (\eqn{temp_{i+1}=temp_{i}^decay}).
#' @param v_ESSth is a positive scalar representing the threshold from which a re-sampling over the set of newly created SMILES is
#' done (0.5 by default). This threshold limits the degeneracy in the set of newly created SMILES. A lower (higher) value allows
#' more (less) degeneracy.
#' @param gentype is the type of the procedure used by the SMILES strings generator. For a Back-off procedure, use "ML"
#' (by default), and for a Neaser-Nay smoothing procedure, use "KN".
#' @param v_maxstock is the maximum of newly created SMILES kept in stock (2000 by default).
#' @param keeptrack is set to TRUE by default. It allows the tracking of the mean of predicted properties, and thus the plotting
#' and/or listing of the latest newly created SMILES during the generation process. It is extremely useful in order to tune the
#' annealing parameters, as to visualize the convergence speed to a targeted physico-chemical properties space.
#' @param smidatabase is a vector of known SMILES to which the generated SMILES should not match. This is useful to avoid the creation
#' of SMILES with great similarity with existing and/or un-wanted ones.
#' @examples
#' \dontrun{#sample data
#' data(qspr.data)
#' idx <- sample(nrow(qspr.data), 5000)
#' smis <- paste(qspr.data[idx,1])
#' y <- qspr.data[idx,c(2,5)]
#'
#' #learning a pattern of chemical strings
#' data(trainedSMI)
#' data(engram_5k)  #same as run => engram <- ENgram$new(trainedSMI, order=10)
#'
#' #learning QSPR model
#' data(qsprpred_EG_5k)
#' #same as run => qsprpred <- QSPRpred$new(smis=smis, y=as.matrix(y), v_fpnames="graph")
#'
#' #set target range
#' targ.min <- c(200,1.5)
#' targ.max <- c(350,2.5)
#' qsprpred_EG_5k$set_target(targ.min,targ.max)
#'
#' #getting chemical strings from the Inverse-QSPR model
#' smchem <- SmcChem$new(smis = rep("c1ccccc1O", 25), v_qsprpred=qsprpred_EG_5k,
#'                      v_engram=engram_5k,temp=3)
#'
#' smchem$smcexec(niter=5, preorder=0, nview=4)
#' #if OpenBabel (>= 2.3.1) is installed, you can use reordering for better mixing as
#' #smchem$smcexec(niter=100, preorder=0.2, nview=4)
#' #see http://openbabel.org
#'
#' #check
#' gensmis <- smchem$get_hiscores(nsmi=5, exsim=0.9)
#' pred <- qsprpred_EG_5k$qspr_predx(gensmis[,1])}
#'
#' @import methods
#' @importFrom graphics abline plot rasterImage
#' @importFrom methods new
#' @export SmcChem
#' @exportClass SmcChem

SmcChem <- setRefClass("SmcChem",
  fields = list(
    esmi = "list",
    oldesmi = "list",
    m = "numeric",
    engram = "ENgram",
    qsprpred = "QSPRpred",
    weights = "numeric",
    oldscore = "numeric",
    fkernel = "numeric",
    bkernel = "numeric",
    isvalid = "logical",
    temp = "numeric",
    ESSth = "numeric",
    decay = "numeric",
    maxstock = "numeric",
    smistock = "character",
    scorestock = "numeric",
    #fpstock = "matrix",
    pytrack = "logical",
    predystock = "matrix",
    predytrack = "list",
    smidb = "character"
  ),
  methods = list(
    initialize = function(smis=NULL, v_engram=NULL, v_m=NULL, v_qsprpred=NULL, v_temp=c(1,1),
                          v_decay=0.95, v_ESSth=0.5, gentype="ML", v_maxstock=2000,
                          keeptrack=TRUE, smidatabase=NULL){
      "Initialize the SMC chemical generator with initial SMILES strings smis, ENgram class object v_engram and QSPRpred class object v_qsprpred"
      if(is.null(v_engram)){
        cat("Need ENgram object\n")
        return(NULL)
      }else{
        engram <<- v_engram
      }

      if(is.null(v_qsprpred)){
        cat("Need QSPRpred object\n")
        return(NULL)
      }else{
        qsprpred <<- v_qsprpred
      }

      if(is.null(v_m)){
        m <<- length(engram$mat)
      }else{
        m <<- v_m
      }

      ESSth <<- v_ESSth

      if(is.null(smis)){
        esmi <<- lapply(rep("c1ccccc1O", 100), function(x) Esmi$new(x, m=m, engram, type=gentype))
      }else{
        esmi <<- lapply(smis, function(x) Esmi$new(x, m=m, engram, type=gentype))
      }

      oldesmi <<- sapply(esmi, function(x) x$copy())

      csmi <- sapply(esmi, function(x) x$get_validsmi())

      temp <<- v_temp
      num_prop <- qsprpred$propndim+length(qsprpred$func)
      # if(dim(temp)[1]==num_prop){
      #   if(dim(temp)[2]!=length(csmi)){
      #     if(dim(temp)[2]==1){
      #       temp <<- matrix(rep(temp,length(csmi)),nrow=num_prop,ncol=length(csmi))
      #       cat("The defined numerical vector of initial temperatures will be applied to the",length(csmi),"annealing processes.\n\n")
      #     } else {
      #       cat("Caution: The provided temperature matrix has a number of columns: ",dim(temp)[2]," that doesn't match the number of initial SMILES: ",length(csmi),".\n\n",sep="")
      #       return(NULL)
      #     }
      #   }
      # } else {
      #   cat("Caution: The provided temperature vector/matrix has a number of elements/lines: ",dim(temp)[1]," that doesn't match the number of properties: ",qsprpred$propndim+length(qsprpred$func),".\n",sep="")
      #   cat("Pay attention that if you have requested an empirical calculation of supplementary properties, you shoud take them into account and apply an individual temperature to each of them.\n\n")
      #   return(NULL)
      # }
      invpredtmp <- qsprpred$iqspr_predict(csmi, temp)
      oldscore <<- invpredtmp[[1]]

      weights <<- rep(1, length(esmi))
      fkernel <<- rep(1, length(esmi))
      bkernel <<- rep(1, length(esmi))
      isvalid <<- rep(TRUE, length(esmi))
      decay <<- v_decay
      maxstock <<- v_maxstock
      smistock <<- character(0)
      scorestock <<- numeric(0)
      #smistock <<- c("")
      #scorestock <<- c(0)
      #fpstock <<- matrix(0, nrow = 1, ncol=ncol(invpredtmp[[4]]))
      pytrack <<- keeptrack
      predystock <<- matrix(0, nrow = 1, ncol=nrow(invpredtmp[[2]])) # ncol = nrow(...) to take the transpose
      predytrack <<- list()
      smidb <<- as.character(smidatabase)
    },

    smcexec = function(niter, nsteps=5, preorder=0, nview=0){
      "modify chemical structures with niter SMC updates"
      for(i in 1:niter){
        localmove(nsteps)
        update_particles()
        reordering(preorder)
        if(nview>0){
          idx <- which(isvalid==T)
          viewstr(idx[1:nview])
        }
      }
    },

    localmove = function(nstep){
      oldesmi <<- sapply(esmi, function(x) x$copy())
      for(i in 1:length(esmi)){
        if(i<length(esmi)){
          cat("\rlocal move for ", i, "th molecules")
        }else{
          cat("\rlocal move for ", i, "th molecules\n")
        }
        u <- rbinom(1, nstep, 0.5)
        rep.times <- min(u,length(esmi[[i]]$vstr)-2)
        if(rep.times<0){
          rep.times <- 0 # rep.times <- u
        }
        steps <- c(rep(1, rep.times), rep(0, nstep-u))
        suppressWarnings(tryres <- try(tempk <- esmi[[i]]$chem_local(engram, steps, 1), silent=T))
        if(class(tryres)=="try-error"){
          fkernel[i] <<- 1
          bkernel[i] <<- 1
          isvalid[i] <<- F
        }else{
          fkernel[i] <<- 1
          bkernel[i] <<- 1
          isvalid[i] <<- TRUE
        }

      }
    },

    reordering = function(prob){
      idxs <- which(sapply(esmi, function(x) nchar(x$get_validsmi()))>1)
      idxs <- intersect(idxs, which(sapply(esmi, function(x) x$numn+x$bcount)==0))
      idxs <- intersect(idxs, which(isvalid==T))
      idxs <- sample(idxs, round(length(idxs)*prob, 0))
      for(idx in idxs){
        natoms <- esmi[[idx]]$get_natoms()
        prevsmi <- esmi[[idx]]$get_validsmi()
        if(.Platform$OS=="windows"){
          suppressWarnings(newsmi <- system(paste('obabel -:"', prevsmi ,'" -osmi -xf ', sample(natoms, 1) , sep=""), intern=T, show.output.on.console = F, ignore.stderr=T))
        }else{
          suppressWarnings(newsmi <- system(paste('obabel -:"', prevsmi ,'" -osmi -xf ', sample(natoms, 1) , sep=""), intern=T, ignore.stderr=T))
        }
        newsmi <- substr(newsmi[1], 1, nchar(newsmi)-1)
        if(idx!=idxs[length(idxs)]){
          cat("\rreordering prev:", prevsmi, "=> new:", newsmi, paste(rep(" ", 40), collapse=" "))
          flush.console()
        }else{
          cat("\rreordering prev:", prevsmi, "=> new:", newsmi, paste(rep(" ", 40), collapse=" "), "\n\n")
          flush.console()
        }
        suppressWarnings(reschk <- try(tempesmi<- Esmi$new(smi=newsmi, m=m, engram), silent=T))
        if(class(reschk)!="try-error"){ # exception for reordering
          temp_vstr <- tempesmi$vstr
          idxc <- grep("=\\[S", temp_vstr)
          if(length(idxc)>0) temp_vstr[idxc] <- sapply(idxc, function(x) paste(substr(temp_vstr[x], 1, 1), substr(temp_vstr[x], 3, 3), sep=""))
          idxc <- grep("=\\[C", temp_vstr)
          if(length(idxc)>0) temp_vstr[idxc] <- sapply(idxc, function(x) paste(substr(temp_vstr[x], 1, 1), substr(temp_vstr[x], 3, 3), sep=""))
          idxc <- grep("=\\[N", temp_vstr)
          if(length(idxc)>0) temp_vstr[idxc] <- sapply(idxc, function(x) paste(substr(temp_vstr[x], 1, 1), substr(temp_vstr[x], 3, 3), sep=""))
          idxc <- grep("=\\[O", temp_vstr)
          if(length(idxc)>0) temp_vstr[idxc] <- sapply(idxc, function(x) paste(substr(temp_vstr[x], 1, 1), substr(temp_vstr[x], 3, 3), sep=""))
          idxc <- grep("\\[S", temp_vstr)
          if(length(idxc)>0) temp_vstr[idxc] <- sapply(idxc, function(x) substr(temp_vstr[x], 2, 2))
          idxc <- grep("\\[C", temp_vstr)
          if(length(idxc)>0) temp_vstr[idxc] <- sapply(idxc, function(x) substr(temp_vstr[x], 2, 2))
          idxc <- grep("\\[N", temp_vstr)
          if(length(idxc)>0) temp_vstr[idxc] <- sapply(idxc, function(x) substr(temp_vstr[x], 2, 2))
          idxc <- grep("\\[O", temp_vstr)
          if(length(idxc)>0) temp_vstr[idxc] <- sapply(idxc, function(x) substr(temp_vstr[x], 2, 2))
          tempesmi$vstr <- temp_vstr
          tempesmi$vstr2smi()


          suppressWarnings(reschk <- try(tempesmi <- Esmi$new(smi=tempesmi$get_validsmi(), m=m, engram), silent=T))
          if(class(reschk)!="try-error"){ # exception for reordering
            ar <- 1
            if(runif(1, 0, 1) < ar){
              esmi[[idx]] <<- tempesmi
              weights[idx] <<- weights[idx]*ar
            }
          }
        }
      }
      weights <<- weights/sum(weights)
    },

    update_particles = function(){
      csmi <- sapply(esmi, function(x) x$get_validsmi())
      isvalid[which(sapply(csmi, function(x) class(try(parse.smiles(x, kekulise=T), silent=T))!="list"))] <<- F
      idx <- which(isvalid==T)
      mols <- parse.smiles(csmi[idx], kekulise=T)
      csmi[idx] <-  sapply(mols, function(x) get.smiles(x, aromatic=T))

      invpredtemp <- qsprpred$iqspr_predict(csmi[idx], temp)

      newscore <- numeric(length(csmi))
      newscore[idx] <- invpredtemp[[1]]
      newscore[-idx] <- 0

      smistock <<- c(smistock, csmi[idx])
      scorestock <<- c(scorestock, newscore[idx])
      #fpstock <<- rbind(fpstock, newfpstock)

      or <- order(scorestock, decreasing=T)
      smistock <<- smistock[or[1:min(length(or),maxstock)]]
      scorestock <<- scorestock[or[1:min(length(or),maxstock)]]
      #fpstock <<- fpstock[or[1:min(length(or),maxstock)],]
      if(pytrack){
        newpredystock <- data.matrix(t(invpredtemp[[2]]))
        if(dim(predystock)[1]==1 && apply(predystock,1,sum)==0){
          predystock <<- newpredystock
        } else {
          predystock <<- data.matrix(rbind(predystock, newpredystock))
        }
        predystock <<- data.matrix(predystock[or[1:min(length(or),maxstock)],])
        predytrack_size <- length(predytrack)
        #predytrack[[predytrack_size+1]] <<- predystock[1:min(length(or),100),] # keep track of 100 predy max. per step
        predytrack[[predytrack_size+1]] <<- data.frame(csmi[idx],invpredtemp[[1]],t(invpredtemp[[2]]),t(invpredtemp[[3]]))
      }

      weights <<- weights*(newscore/oldscore)
      weights[which(oldscore==0)] <<- 0
      weights <<- weights/sum(weights)

      ESS <- 1/sum(weights^2)
      if(ESS < ESSth*length(esmi)){
        cat("weight update done, ESS=", ESS, "\n")
        idx <- sample(1:length(esmi), prob=weights, replace=T)
        esmi <<- lapply(esmi[idx], function(x) x$copy())
        weights <<- rep(1/length(esmi), length(esmi))
        oldscore <<- newscore[idx]
        cat("resampling done\n\n")
      }else{
        cat("weight update done, ESS=", ESS, "\n\n")
        oldscore <<- newscore
      }
      temp <<- temp^decay
    },

    get_smiles = function(){
      "get SMILES strings from the SmcChem object (same as get_smiles function) "
      csmi <- sapply(esmi, function(x) x$get_validsmi())
      idx <- which(isvalid==T)
      mols <- parse.smiles(csmi[idx], kekulise=T)
      csmi[idx] <-  sapply(mols, function(x) get.smiles(x, aromatic=T))
      csmi
    },

    get_hiscores = function(nsmi=100, exsim=0.8){
      "get chemical structures with high QSPR score from SmcChem object (same as get_hiscores function) "
      SMILES <- character(0)
      QSPRScore <- numeric(0)
      tsmi <- smistock
      tscore <- scorestock

      if(!is.null(smidb)){
        smired <- tsmi %in% smidb # Check if newly created molecules already exist in the database
        tsmi <- tsmi[!smired]
        tscore <- tscore[!smired]
      }

      fp <- get_descriptor(smis = paste(tsmi), desctypes = c("standard","extended","circular"))[[1]]
      # To be used if
      #fpnames corresponds to fingerprints alone
      #fp <- fpstock # Not wanted for a Tanimoto similarity check
      cat("\n")
      j <- 1
      while((j<=nsmi) & (length(tsmi)>2)){
        SMILES <- c(SMILES, tsmi[1])
        QSPRScore <- c(QSPRScore, tscore[1])
        sim <- apply(fp[-1,], 1, function(x) sum(x*fp[1,])/sum((x+fp[1,])>=1)) # Tanimoto similarity
        #print(sim)
        fp <- fp[-1,]
        tsmi <- tsmi[-1]
        if(length(tsmi)>1){
          if(sum(sim>=exsim)>0){
            fp <- fp[-which(sim>=exsim),]
            tsmi <- tsmi[-which(sim>=exsim)]
            tscore <- tscore[-which(sim>=exsim)]
          }
        }
        cat("\r", j, "th molecules is chosen")
        j <- j + 1
      }
      cat("\n")
      return(cbind(SMILES, QSPRScore))
    },

    viewstr = function(idx){
      "view 2D structures from SMILES string vector with index idx (same as viewstr function) "
      tsmi <- sapply(smchem$esmi, function(x) x$get_validsmi())[idx]
      nrow = ncol = ceiling(sqrt(length(idx)))
      par(mar = c(0,0,0,0))
      par(mfrow = c(nrow,ncol))
      smisl <- length(tsmi)
      rcdk_version <- as.numeric(gsub("\\.","",packageVersion("rcdk")))
      for(i in 1:smisl){
        mol <- parse.smiles(tsmi[i], kekulise = T)[[1]] # if kekulise = F, aromatic rings are missed!
        if(rcdk_version > 338){
          dep <- get.depictor(width=500, height=500, zoom=3)
          temp <- view.image.2d(molecule=mol,depictor=dep)
        } else {
          temp <- view.image.2d(molecule=mol,width=500, height=500)
        }
        plot(NA,NA,xlim=c(1,10),ylim=c(1,10),xaxt='n',yaxt='n',xlab='',ylab='')
        rasterImage(temp,1,1,10,10)
        legend("topleft", legend[i], bty = "n", cex = 1.3)
      }
    }
  )
)
