# Header #############################################################
# 
# Author: Lisa Nicvert
# Email:  lisa.nicvert@univ-lyon1.fr
# 
# Date: 2023-03-09
#
# Script Description: compendium functions

#' Install and load a package (if not already loaded)
#'
#' @param x The package name
#'
#' @return Installs the package if it was not already installed
#' @export
require <- function(x) { 
  
  if (!base::require(x, character.only = TRUE)) {
    install.packages(x, dep = TRUE) ; 
    base::require(x, character.only = TRUE)
  }
}