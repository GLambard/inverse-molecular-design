#' Extend smi class 
#' @description Extend smi class
#' 
#' @import methods
#' @importFrom methods new
#' @export Esmi
#' @exportClass Esmi
#'
Esmi <- setRefClass(
  Class = "Esmi",
  
  fields = list (
    smi = "character",
    vstr = "character",
    hsubstr = "list",
    prevb = "character",
    prev = "character",
    prevbr = "list",
    prevtemp = "character",
    
    m = "numeric",
    type="character",
    btflag = "logical",
    numchk = "numeric",
    numstk = "numeric",
    bcount = "numeric",
    startbr = "logical",
    brsteps = "numeric",
    numn = "numeric",
    brpos = "numeric"
  ),
  
  methods = list (
    initialize = function(smi="CCCCCCC", m=5, mnmat=NULL, type="ML") {
      smi <<- smi
      m <<- m
      type <<- type
      
      if(!is.null(mnmat)){
        get_vstr()
        get_hsubstr(mnmat)
        
        prevb <<- hsubstr[[length(hsubstr)]]
        prev <<- prevb[length(prevb)]
        
        update_prevbr()
        prevtemp <<- character(0)
        
        btflag <<- F
        numchk <<- numeric(100)
        numstk <<- numeric(0)
        bcount <<- 0
        startbr <<- F
        brsteps <<- 0
        numn <<- 0
        brpos <<- 0
      }
    },
        
    get_vstr = function(){
      vstr <<- sapply(1:nchar(smi), function(x) substr(smi, x, x))
      vstr <<- c(vstr, "\t")
      while(1){
        p1 <- which(vstr=="[")
        p2 <- which(vstr=="]")
        if(length(p1)>0){
          vstr[p1[1]] <<- paste(vstr[(p1[1]):p2[1]], collapse="")
          vstr <<- vstr[-((p1[1]+1):p2[1])]
        }else{
          break
        }
      }
      p <- which(vstr=="=")
      if(length(p)>0){
        vstr[p] <<- paste(vstr[p], vstr[p+1], sep="")
        vstr <<- vstr[-(p+1)]
      }
      
      stknum <- c()
      idx <- grep("^[1-9]$", vstr)
      for(i in idx){
        if(sum(vstr[i]==stknum)==0){
          stknum <- c(stknum, vstr[i])
          vstr[i] <<- "#num#"
        }else{
          pos <- which(vstr[i]==stknum)
          diff <- length(stknum) - pos
          stknum <- stknum[-pos]
          vstr[i] <<- paste("#numt_", diff, sep="")
        }
      }
      idx <- which(vstr==")")
      vstr[idx] <<- "#brterm#"
      vstr <<- vstr[-length(vstr)]
    },
    
    get_hsubstr = function(mnmat){
      res <- c()
      t_bcount <- 0
      t_numn <- 0
      prevb <<- character(0)
      stkbr <- c()
      for(i in 1:length(vstr)){
        if(length(prevb)<m){
          prevb <<- c(prevb, vstr[i])
        }else{
          prevb <<- c(prevb[-1], vstr[i])
        }
        if(i>3){
          if(vstr[i-2]=="("){
            stkbr <- c(stkbr, list(c(prevb[-length(prevb)], "#brterm#")))
          }
        }
        if(vstr[i]=="#brterm#"){
          prevb <<- stkbr[[length(stkbr)]]
          stkbr <- stkbr[-length(stkbr)]
        }
        if(i < length(vstr)){
          if(vstr[i]=="("){
            t_bcount <- t_bcount + 1
          }
          if(vstr[i]=="#brterm#"){
            t_bcount <- t_bcount - 1
          }
          if(vstr[i]=="#num#"){
            t_numn <- t_numn + 1
          }
          if(length(grep("#numt_", vstr[i]))==1){
            t_numn <- t_numn - 1
          }
          matidx <- 10*(t_bcount>0)+t_numn+1
        }
        res <- c(res, list(prevb))
      }
      hsubstr <<- res
    },
    
    update_prevbr = function(){
      prevbr <<- list()
      idx <- which(vstr=="(")
      for(i in idx){
        if(i <= length(vstr)-2){
          start <- max(0, i+3-m)
          end <- i+1
          prevbr <<- c(prevbr, list(c(vstr[start:end], "#brterm")))
        }
      }
    },
    
    parse = function(){
      stknum <- c()
      idx <- grep("^[1-9]$", vstr)
      for(i in idx){
        if(sum(vstr[i]==stknum)==0){
          stknum <- c(stknum, vstr[i])
          vstr[i] <<- "#num#"
        }else{
          pos <- which(vstr[i]==stknum)
          diff <- length(stknum) - pos
          stknum <- stknum[-pos]
          vstr[i] <<- paste("#numt_", diff, sep="")
        }
      }
      idx <- which(vstr==")")
      vstr[idx] <<- "#brterm#"
      vstr <<- vstr[-length(vstr)]
    },
    
    vstr2smi = function(){
      res <- vstr
      res[which(res=="#brterm#")] <- ")"
      t_numstk <- c()
      t_numchk <- numeric(100)
      for(i in grep("#num", vstr)){
        if(vstr[i]=="#num#"){
          num <- min(which(t_numchk==0))
          t_numstk <- c(t_numstk, num)
          t_numchk[num] <- 1
          res[i] <- num
        }else{
          rn <- as.numeric(substr(res[i], 7,7))
          num <- t_numstk[length(t_numstk)-rn]
          t_numstk <- t_numstk[-(length(t_numstk)-rn)]
          t_numchk[num] <- 0
          res[i] <- num
        }
      }
      smi <<- paste(res, collapse="")
    },
    
    get_prob = function(v_prevb, mnmat, matidx){
      tm <- length(v_prevb)
      rowidx <- sapply(1:tm, function(x) find_rownames(paste(v_prevb[(tm-x+1):tm], collapse=""), mnmat$mat[[x]][[matidx]]))
      maxd <- max(which(rowidx>0))
      maxd <- max(which(sapply(1:maxd, function(x) sum(mnmat$mat[[x]][[matidx]][rowidx[x],])>0)))
      pnc <- lapply(maxd:1, function(x) mnmat$mat[[x]][[matidx]][rowidx[x],])
      ld <- sapply(maxd:1, function(x) mnmat$denom[[x]][[matidx]][rowidx[x]])
      
      ln <- sapply(maxd:1, function(x) mnmat$N1p[[x]][[matidx]][rowidx[x]])
      
      if(type=="NN"){
        list_NN <- lapply(1:length(pnc), function(x) pnc[[x]]-mnmat$mat_D[x,matidx])
        list_idx <- lapply(list_NN, function(x) which(x<0))
        for(i in 1:length(list_NN)){
          list_NN[[i]][list_idx[[i]]] <- 0
        }
        for(i in maxd:1){
          if(i==maxd){
            probNN <- list_NN[[i]]/ld[i]
            names(probNN) <- names(list_NN[[i]])
          }else{
            idx <- names(list_NN[[i]])
            probNN[idx] <- probNN[idx]*ln[i]*mnmat$mat_D[i,matidx]/ld[i] + list_NN[[i]][idx]/ld[i]
          }
        }
        probNN <- probNN/sum(probNN)
        return(probNN)
      }else{
        probNN <- pnc[[1]]
        return((probNN/sum(probNN)))
      }
    },
    
    chem_local = function(mnmat, steps, temp){
      res <- c(1,1)
      for(s in steps){
        if(s==0){
          matidx <- 10*(bcount>0)+numn+1          
          probNN <- get_prob(prevb, mnmat, matidx)
          prev <<- sample(names(probNN), 1, prob=probNN)
          res[1] <- res[1] * probNN[prev]
          if(prev!="#term#"){
            if(startbr==T){
              if(brsteps==1){
                tprevb <- c(prevb[-1], "#brterm#")
                prevbr <<- c(prevbr, list(tprevb)) 
                startbr <<- F
                brsteps <<- 0
              }else{
                brsteps <<- brsteps + 1
              }
            }
            if(prev=="("){
              startbr <<- T
              brsteps <<- 0
              bcount <<- bcount + 1
            }else if(prev=="#num#"){
              numn <<- numn + 1
            }else if(prev=="#brterm#"){
              bcount <<- bcount - 1
              btflag <<- T
            }else if(length(grep("#numt_", prev))==1){
              numn <<- numn - 1
            }
            if(btflag==T){
              prevb <<- prevbr[[length(prevbr)]]
              prevbr <<- prevbr[-length(prevbr)]
              btflag <<- F
            }else if(length(prevb)==m){
              prevb <<- c(prevb[-1], prev)
            }else{
              prevb <<- c(prevb, prev)
            }
            hsubstr <<- c(hsubstr, list(prevb))
            vstr <<- c(vstr, prev)
          }
        }else{
          if(startbr == T){
            brsteps <<- brsteps - 1
          }
          lastc <- vstr[length(vstr)]
          prev <<- lastc
          if(lastc=="#brterm#"){
            prevbr <<- c(prevbr, list(prevb))
            bcount <<- bcount + 1
          }else if(lastc=="("){
            bcount <<- bcount - 1
            if(startbr==T){
              startbr <<- F
              brsteps <<- 0
            }else{
              prevbr <<- prevbr[-length(prevbr)] 
            }
          }else if(length(grep("#numt_", lastc))==1){
            numn <<- numn + 1
          }else if(lastc=="#num#"){
            numn <<- numn - 1
          }
          hsubstr <<- hsubstr[-length(hsubstr)]
          prevb <<- hsubstr[[length(hsubstr)]]
          vstr <<- vstr[-length(vstr)]

          matidx <- 10*(bcount>0)+numn+1
          probNN <- get_prob(prevb, mnmat, matidx)
          res[2] <- res[2] * probNN[lastc]  
        }
      }
      vstr2smi()
      return(res)
    },
    
    get_validsmi = function(){
      res <- ""
      chars <- sapply(1:nchar(smi), function(x) substr(smi, x, x))
      while(1){
        p1 <- which(chars=="[")
        p2 <- which(chars=="]")
        if(length(p1)>0){
          chars[p1[1]] <- paste(chars[(p1[1]):p2[1]], collapse="")
          chars <- chars[-((p1[1]+1):p2[1])]
        }else{
          break
        }
      }  
      
      if(chars[length(chars)]=="("){
        chars <- chars[-(length(chars))]
      }
      np <- sapply(1:9, function(x) which(chars==x))
      odd <- which(sapply(np, length)%%2==1)
      if(length(odd)>0){
        rmlist <- c()
        for(i in odd){
          rmlist <- c(rmlist, np[[i]][length(np[[i]])])
        }
        chars <- chars[-(rmlist)] 
      }
      
      if(length(chars)==1){
        chars[1] <- toupper(chars[1])
      }else{
        n_num <- 0
        sq_num <- 0
        flag <- F
        numlist <- c()
        for(i in 1:(length(chars)-1)){
          if(chars[i]=="["){
            sq_num <- sq_num + 1
          }      
          if(chars[i]=="]"){
            sq_num <- sq_num - 1
          } 
          if(sq_num==0){
            if(length(grep("[1-9]", chars[i+1]))>0){
              if(sum(chars[i+1]==numlist)>0){
                idx <- which(numlist==chars[i+1])
                numlist <- numlist[-idx]
                flag <- T
              }else{
                n_num <- n_num + 1
                numlist <- c(numlist, chars[i+1])
              }
            }
            if((n_num == 0) & (sq_num == 0)){
              if(sum(chars[i]==c("c", "s" ,"n", "o"))>0){
                chars[i] <- toupper(chars[i])
              }
            }
            if(flag == T){
              n_num <- n_num -1
              flag <- F
            }
          }
        }
        if(length(grep("[1-9]", chars[length(chars)]))==0){
          if(sum(chars[length(chars)]==c("c", "s" ,"n", "o"))>0){
            chars[length(chars)] <- toupper(chars[length(chars)])
          }
        }
      }
      res <- paste(chars, collapse="")
      
      brstart <- sum(chars=="(")
      brend <- sum(chars==")")
      sbr <- brstart - brend
      while(sbr>0){
        res <- paste(res, ")", sep="")
        sbr <- sbr - 1
      }
      return(res)
    },
    
    get_natoms = function(){
      res <- 0
      res <- res + length(grep("[A-Z]" ,vstr))
      res <- res - length(grep("H" ,vstr))
      res <- res + length(grep("c" ,vstr))
      res <- res + length(grep("n" ,vstr))
      res <- res + length(grep("o" ,vstr))
      return(res-1)
    },
    
    get_nC = function(){
      res <- 0
      res <- res + sum(vstr=="C")
      res <- res + sum(vstr=="c")
      res <- res + sum(vstr=="=C")
      return(res)
    },
    
    get_pA = function(){
      res <- 0
      res <- res + sum(vstr=="N")
      res <- res + sum(vstr=="n")
      res <- res + sum(vstr=="=N")
      res <- res + sum(vstr=="F")
      res <- res + sum(vstr=="Cl")
      res <- res + sum(vstr=="Br")
      return(res)
    },
    
    find_rownames = function(input, mat, w=25){
      chars <- rownames(mat)
      leftidx <- 1
      rightidx <- length(chars)
      width <- trunc(rightidx/2, 0)
      while(width>w){
        if(input >= chars[leftidx+width]){
          leftidx <- leftidx + width
        }
        if(input <= chars[rightidx-width]){
          rightidx <- rightidx - width
        }
        width <- trunc(0.5*width, 0)
      }
      res <- (which(chars[leftidx:rightidx]==input)+leftidx-1)
      if(length(res)>0){
        return(res)
      }else{
        return(0)
      }
    }
  )
)