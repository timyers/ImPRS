# Purpose:  calculate polygenic risk scores
#
# Dependencies = BEDMatrix package
#

library(BEDMatrix)

####### helper code to delete later
base_file <- "./data/Height.QC.Transformed_trunc"
base_header = TRUE
base_data <- read.table(file = base_file, header = base_header)

target_bed <- "./data/EUR.QC.bed"
target_data <- BEDMatrix::BEDMatrix(target_bed)



#####################################


# Primary PRS function
prs_calc <- function(base_file = NULL, base_header = TRUE, target_bed = NULL) {
  
  # read the 'base' data file
    base_data <- read.table(file = base_file, header = base_header)
  # read the 'target' data file
    target_data <- BEDMatrix::BEDMatrix(target_bed)
  
  # Rename the column names in 'table_data' to match the rsID format
  # in 'base_data'
    colnames(target_data) = gsub(pattern = "_.", 
                                 replacement = "", 
                                 x = colnames(target_data)
                                )
   
  # Transpose 'target_data; make rows rsID's and columns sampleID's; same
  # as 'base_data' ('target_data' must be a matrix)
    target_data <- t(as.matrix(target_data))
  
  # Subset 'base_data' with only rsID's and OR(weight)
    base_data_subset <- data.frame(SNP = base_data$SNP, OR = base_data$OR)

    
    
}