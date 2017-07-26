#'view 2D structures from SMILES string vector
#' @description View 2D structures from SMILES strings vector
#' @param smis is a SMILES strings vector to submit.
#'
#' @examples viewstr(c("c1ccc2ccc3c(NCCN(C)C)cc(nc3c2c1)", "c1ccc2ccc3c(NCCN(CC)CCCl)cc(nc3c2c1)",
#'  "c1ccc2ccc3c(NC(CC)CC)cc(nc3c2c1)", "c1ccc2ccc3c(c2c1)ncc(c3NCCNCC=CCCCC)"))
#'
#' @return a plot of requested 2D structures.
#'
#' @export viewstr

viewstr <- function(smis = NULL, width = 500, height = 500, nrow = 2, ncol = 2, legend = ""){
  par(mar = c(0,0,0,0))
  par(mfrow = c(nrow,ncol))
  smisl <- length(smis)
  rcdk_version <- as.numeric(gsub("\\.","",packageVersion("rcdk")))
  for(i in 1:smisl){
    mol <- parse.smiles(smis[i], kekulise = T)[[1]] # if kekulise = F, aromatic rings are missed!
    if(rcdk_version > 338){
      dep <- get.depictor(width=width, height=height, zoom=3)
      temp <- view.image.2d(molecule=mol,depictor=dep)
    } else {
      temp <- view.image.2d(molecule=mol,width=width, height=height)
    }
    plot(NA,NA,xlim=c(1,10),ylim=c(1,10),xaxt='n',yaxt='n',xlab='',ylab='')
    rasterImage(temp,1,1,10,10)
    legend("topleft", legend[i], bty = "n", cex = 1.3)
  }
}
