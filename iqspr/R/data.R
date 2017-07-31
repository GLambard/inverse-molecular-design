#' List of 5,000 SMILES used for training the n-grams model
#'
#' A dataset containing SMILES of 5,000 compounds from PubChem
#'
#' @format A character vector with 5,000 rows of 1 variable:
#' \describe{
#'   ...
#' }
#' @source \url{https://pubchem.ncbi.nlm.nih.gov/}
"trainedSMI"

#' n-grams model learned from 5,000 SMILES from PubChem
#'
#' A dataset containing up to order 10 n-grams issued from 5,000 SMILES
#'
#' @format A ENgram object only suitable for the iqspr
#' \describe{
#'   ...
#' }
#' @source \url{https://pubchem.ncbi.nlm.nih.gov/}
"engram_5k"

#' Table of 16,674 SMILES with associated internal energy E and HOMO-LUMO gap
#'
#' A dataset containing SMILES of 16,674 compounds from PubChem with
#' associated internal energy E and HOMO-LUMO gap
#'
#' @format A data.frame with 16,674 rows of 3 variables:
#' \describe{
#'   \item{SMILES}{SMILES, in string of characters}
#'   \item{E}{internal energy, in kJ/mol}
#'   \item{HOMO-LUMO gap}{HOMO-LUMO gap, in eV}
#'   ...
#' }
#' @source \url{https://pubchem.ncbi.nlm.nih.gov/}
"qspr.data"

#' List of predictions associated to the internal energy E and HOMO-LUMO gap
#' of 500 SMILES
#'
#' A list containing two matrices of predicted internal energy E and HOMO-LUMO gap
#' for 500 SMILES from linear and random forest models respectively
#'
#' @format A list with two matrices of 500 rows of 2 variables:
#' \describe{
#'   ...
#' }
#'
"predictions"

#' Table of 372 SMILES with associated QSPRScore
#'
#' A dataset containing 372 generated SMILES with their QSPRScore
#'
#' @format A matrix with 372 rows of 2 variables:
#' \describe{
#'   \item{SMILES}{SMILES, in string of characters}
#'   \item{QSPRScore}{generation score, in a.u.}
#'   ...
#' }
#'
"gensmis"

#' Table of 372 SMILES with predicted internal energy E and HOMO-LUMO gap
#'
#' A dataset containing the predicted internal energy E and HOMO-LUMO gap of
#' 372 generated compounds with iqspr
#'
#' @format A data.frame with 372 rows of 2 variables:
#' \describe{
#'   \item{E}{internal energy, in kJ/mol}
#'   \item{HOMO-LUMO gap}{HOMO-LUMO gap, in eV}
#'   ...
#' }
#'
"dpred"

#' Table of predicted internal energy E and HOMO-LUMO gap for the phenol
#'
#' A dataset containing the predicted internal energy E and HOMO-LUMO gap
#' for the phenol
#'
#' @format A matrix with 1 row of 2 variables:
#' \describe{
#'   \item{E}{internal energy, in kJ/mol}
#'   \item{HOMO-LUMO gap}{HOMO-LUMO gap, in eV}
#'   ...
#' }
#'
"predinit"
