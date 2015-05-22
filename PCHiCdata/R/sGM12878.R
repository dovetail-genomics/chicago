#' Promoter Capture data for chromosomes 20 and 21 in GM12878
#'
#' Promoter Capture data to be used as a toy example to run all steps of Chicago. This data only incorporates read pairs
#' including both chromosomes 20 and 21 in order to minimize processing time and memory storage. Thus, this data set 
#' includes all cis read pairs for these two chromosomes and all trans read pairs between them.
#' 
#' @format A ChicagoData object.
#' 
#' @source Mifsud, B. et al. "Mapping long-range promoter contacts in human cells with high-resolution capture Hi-C." Nature Genetics (2015) doi:10.1038/ng.3286
#' 
#' @seealso \code{\link{chicagoData}}
#' 
#' @examples
#' data(sGM12878)
#' ##modifications to sGM12878, ensuring it uses correct design directory
#' designDir <- file.path(system.file("extdata", package="PCHiCdata"), "hg19TestDesign")
#' new.settings <- paste0(file.path(designDir, "hg19_chr20and21"), c(".rmap", ".baitmap", ".npb", ".nbpb", ".poe"))
#' new.settings <- as.list(new.settings)
#' names(new.settings) <- c("rmapfile", "baitmapfile", "nperbinfile", "nbaitsperbinfile", "proxOEfile")
#' sGM12878 <- modifySettings(cd = sGM12878, new.settings)
#' 
"sGM12878"