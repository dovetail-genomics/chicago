#' Promoter Capture data for chromosomes 18 and 19 in smESC
#'
#' Promoter Capture data to be used as a toy example to run all steps of Chicago. This data only incorporates read pairs
#' including both chromosomes 18 and 19 in order to minimize processing time and memory storage. Thus, this data set 
#' includes all cis read pairs for these two chromosomes and all trans read pairs between them.
#' 
#' 
#' @format A ChicagoData object.
#' 
#' @source Schoenfelder, S. et al. "The pluripotent regulatory circuitry connecting promoters to their long-range interacting elements." Genome research 25.4 (2015): 582-597.
#'  
#' @seealso \code{\link{chicagoData}}
#' 
#' @examples
#' data(smESC)
#' ##modifications to smESC, ensuring it uses correct design directory
#' designDir <- file.path(system.file("extdata", package="PCHiCdata"), "mm9TestDesign")
#' smESC <- modifySettings(cd = smESC, designDir = designDir)
"smESC"
