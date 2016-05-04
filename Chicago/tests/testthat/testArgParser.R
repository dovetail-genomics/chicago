library(argparser)

test_that("argparser works in the way that the runChicago.R script (in chicagoTools) expects", {
  
  outputFromCommandArgs <- c("--output-dir","my-out-dir","--en-feat-folder","some/dir/thing.txt",
                             "--export-format", "seqMonk,washU_text",
                              "some/directory/file.chinput", "test")
  args <- outputFromCommandArgs
    
  ##runChicago.R code as of 0.99.8
  spec = matrix(c("<input-files>", "Full path to the input file (or comma-separated list of files)", 
                  "<output-prefix>", "Output file names prefix (cannot contain folders)"),  byrow=T, ncol=2)
  
  p = arg_parser("Run Chicago from input files", name="Rscript runChicago.R")
  p = add_argument(p, arg=spec[,1], help=spec[,2])
  
  p = add_argument(p, arg="--settings-file", help = "Full path to Chicago settings file", default = NA)
  p = add_argument(p, arg="--design-dir", 
                   help = "Folder with capture design files (note the settings file has priority over these)", default = "")
  
  p = add_argument(p, arg="--print-memory", help = "Should chicagoPipeline print out memory use?", flag=T)
  
  p = add_argument(p, arg="--cutoff", help = "Score cutoff for writing out peaks and testing feature enrichment", default = 5)
  p = add_argument(p, arg="--export-format", 
                   help = "File format for writing out peaks: one or more of the following: seqMonk,interBed,washU_text,washU_track (comma-separated)", 
                   default = "washU_text")
  p = add_argument(p, arg="--export-order", help = "Should the results be ordered by \"score\" or genomic \"position\"?", 
                   default = "position")
  
  p = add_argument(p, arg="--rda", help = "Save the Chicago object as an RDa image (instead of the default RDS)", flag = T)
  p = add_argument(p, arg="--save-df-only", help = "Save only the data part of the Chicago object, as a data frame (for compatibility)", flag = T)
  
  p = add_argument(p, arg="--examples-prox-dist", help = "The distance limit for plotting \"proximal\" examples", default=1e6)
  p = add_argument(p, arg="--examples-full-range", help = "Also plot examples for the full distance range", flag = T)
  
  p = add_argument(p, arg="--output-dir", help = "The name of the output directory (can be a full path)", default="<output-prefix>")
  
  p = add_argument(p, arg="--en-feat-files", 
                   help = "A comma-separated list of files with genomic feature coordinates for computing peaks' enrichment", 
                   default = NA)
  p = add_argument(p, arg="--en-feat-list", 
                   help = "Same as above but the supplied file contains the feature names and 
                 locations of feature files (in the format <feature-name> <feature-file-location>", default = NA)
  p = add_argument(p, arg="--en-feat-folder", 
                   help = "The folder, in which all feature files are located (if provided, --en-feature-file(s) don't need to list the full path)", 
                   default=NA)
  
  p = add_argument(p, arg="--en-min-dist", help = "The lower distance limit for computing enrichment for features", default="0")
  p = add_argument(p, arg="--en-max-dist", help = "The upper distance limit for computing enrichment for features", default=1e6)
  p = add_argument(p, arg="--en-full-cis-range", help = "Assess the enrichment for features for the full distance range [same chromosome only; use --en-trans in addition to include trans-interactions] (can be very slow!)", flag = T)
  p = add_argument(p, arg="--en-sample-no", help = "The number of negative samples for computing enrichment for features", default=100)
  p = add_argument(p, arg="--en-trans", help = "Include trans-interactions into enrichment analysis", flag=T)
  
  p = add_argument(p, arg="--features-only", help = "Re-run feature enrichment analysis with Chicago output files. With this option, <input-files> must be either a single Rds file (must contain full Chicago objects) or '-', in which case the file location will be inferred automatically from <output-prefix> and files added to the corresponding folder.",  flag = T)
  
  opts = parse_args(p, args)
  
  expect_identical(opts$"<input_files>", "some/directory/file.chinput")
  expect_identical(opts$"<output_prefix>", "test")
  expect_identical(opts$"en_feat_folder", "some/dir/thing.txt")
  expect_identical(opts[["output_dir"]], "my-out-dir")
})