# Purpose:  calculate polygenic risk scores
#
# Dependencies = BEDMatrix package
#

library(BEDMatrix)

####### helper code to delete later
base <- "./data/Height.QC.Transformed_trunc"
base_header = TRUE
base_data <- read.table(file = base, header = base_header)

target_root <- ".data/EUR.QC"
target_bed <- paste(target_root, ".bed", sep = "")



#####################################


# Primary PRS function
prs_calc <- function(base = NULL, base_header = TRUE, target_root = NULL) {
  
  # read the 'base' data file
  base_data <- read.table(file = base, header = base_header)
  
  # read the 'target' data file
  ## paste plink filenames with .bed, .bim, .fam files extensions,
  ## using 'target_root' argument
  target_bed <- paste(target_root, ".bed", sep = "")
  target_bim <- paste(target_root, ".bim", sep = "")
  target_fam <- paste(target_root, ".fam", sep = "")
  
  
  
}