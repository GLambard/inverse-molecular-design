#'generate SMILES strings from extended N-gram model
#' @description Generate SMILES strings from an extended N-gram model
#' @param nsmis is the number of SMILES strings to generate
#' @param engram is an ENgram object
#' @param order n in the ENgram model
#' @param gentype is the type of the procedure used by the SMILES strings generator. For a Back-off procedure, use "ML" (by default),
#' and for a Neaser-Nay smoothing procedure, use "KN".
#' @param crange is the range of lengths, defined in a colum vector, required for the output SMILES strings (from 5 to 10 characters
#' by default)
#' @examples \dontrun{data(engram_5k)
#' smiles <- genENgram(4, engram_5k, 10)
#' viewstr(smiles)}
#'
#' @export  genENgram
genENgram <- function(nsmis, engram, order, gentype="ML", crange=c(5, 10)){
  res <- character(nsmis)
  for(i in 1:nsmis){
    cat("\r", i, "th molecules generated")
    flag <- F
    while(!flag){
      tsmi <- Esmi$new("C", m=order, engram, type=gentype)
      while(tsmi$prev!="#term#"){
        tsmi$chem_local(engram, 0, 1)
      }
      if(class(try(parse.smiles(tsmi$get_validsmi(), kekulise=T), silent=T))=="list"){
        smilength <- nchar(tsmi$get_validsmi())
        if( (smilength >= crange[1]) & (smilength <= crange[2]) ){
          res[i] <- tsmi$get_validsmi()
          flag <- T
        }
      }
    }
  }
  cat("\n")
  return(res)
}
