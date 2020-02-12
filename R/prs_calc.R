# Purpose:  calculate polygenic risk scores
#
# Dependencies = BEDMatrix package
#

# library(BEDMatrix)

####### helper code to delete later
# base_file <- "./data/Height.QC.Transformed_trunc"
# base_header = TRUE
# base_data <- read.table(file = base_file, header = base_header)

# target_bed <- "./data/EUR.QC.bed"
# target_data <- BEDMatrix::BEDMatrix(target_bed)



#####################################


########### Primary PRS function ###########

#' Calculate Polygenic Risk Score
#'
#' @param base_file path and file name to summary statistics
#' @param base_header a logical value indicating whether the first row of base_file contains column names, default is TRUE.
#' @param target_bed path and file name to plink's .bed individual-level genotype-phenotype
#' file, .bim and .bim files must also be located in the same directory with the same file name,
#' prefix, e.g. EUR.bed, EUR.bim, EUR.fam.
#'
#' @return a numeric
#' @export
#'
#' @examples 
#' prs_calc(base_file = "./data/extdat/Height.QC.Transformed", target_bed = "./inst/extdat/EUR.QC.bed")
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
    
  # The list of variants in 'base_data_subset' and 'target_data' need
  # to be the same
    base_data_subset_match <- base_data_subset[base_data_subset$SNP %in% row.names(target_data),]
    target_data_match <- target_data[rownames(target_data) %in% as.character(base_data_subset_match$SNP),]
  # Now we have a matrix, ‘target_data_subset_match’, w/ SNP-score (0, 1 or 2), and
  # a data.frame, ‘base_data_subset_match’, with the weight (transformed OR),
  # that have the same list of variants.
    
  ##### Calculating PRS #####
  
  # Multiply the 'weight' (OR) by the SNP-score (0,1,2)
    df_product <- data.frame(apply(target_data_match, 
                                   2, 
                                   function(x) x * base_data_subset_match$OR
                                  ),
                             check.names = FALSE
                            )
    
  # Calculate the PRS for each individual
    df_prs <- colSums(df_product)
    
  return(df_prs)  
}
########### End Primary Function ##########

#' @importFrom("utils", "read.table")