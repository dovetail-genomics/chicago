#!/usr/bin/env Rscript

pkgs <- c("data.table", "tidyverse", "argparser")
suppressWarnings(lapply(pkgs, library, character.only = TRUE))

args = commandArgs(trailingOnly=T)

p <- arg_parser("Text formatting function for conversion WashU file compatible with new browser", name = "WashuOld2New.R", hide.opts = TRUE)
p <- add_argument(parser = p, arg = "--washuFile", 
                  help="Full path to the washU_text.txt CHiCAGO produces as an output", default = NA)
p <- add_argument(parser = p, arg = "--outputDir", help = "Full path to the output directory", default = ".")
p <- add_argument(parser = p, arg = "--name", help = "Specify name of the output file .txt", default = "WashUNEW_text.txt")
argv <- parse_args(p, args)

WashU.file <- read.table(argv$washuFile)
outputDir <- argv$outputDir
outputName <- argv$name

## x = CHiCAGO WashU.txt output file

WashU.new.f <- function(WashU.file, outputDir, outputName) {
  
  a <- WashU.file %>% tidyr::separate(V1, into =c("chr", "start", "end"), sep = ",", remove = TRUE) %>% 
    tidyr::separate(V2, into =c("chr.2", "start.2", "end.2"), sep = ",", remove = TRUE)  %>% 
    unite(col = "V2",chr.2, start.2, sep = ":", remove = TRUE) %>%
    unite(col = "V2",V2, end.2, sep = "-", remove = TRUE) %>%
    unite(col = "V2",V2, V3, sep = ",", remove = FALSE) %>%
    mutate(strand = ".")
  
  write.table(a, paste0(outputDir, "/", outputName), quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
  
}

WashU.new.f(WashU.file, outputDir, outputName)