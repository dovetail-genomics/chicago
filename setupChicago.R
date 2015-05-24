message("\nRscript setupChicago.R [--path=<chicago-package-path>] --rlib=<r-lib-dir> --bin=<scripts-target-dir>\n\nChicago and dependencies installation script")
args = commandArgs(trailingOnly = TRUE)

pathLoc = grep("\\-\\-path\\=", args)
path = NULL
if (length(pathLoc)==1){
  path = gsub("\\-\\-path\\=", "", args[pathLoc])
}

rlibLoc = grep("\\-\\-rlib\\=", args)
rlib = NULL
if (length(rlibLoc)==1){
  rlib = gsub("\\-\\-rlib\\=", "", args[rlibLoc])
}

binLoc = grep("\\-\\-bin\\=", args)
bin = NULL
if (length(binLoc)==1){
  bin = gsub("\\-\\-bin\\=", "", args[binLoc])
}

if (length(args)>length(pathLoc)+length(binLoc)+length(rlibLoc)){
  stop(paste("Unknown options provided or no value specified:", args[!(1:length(args)) %in% c(pathLoc, binLoc, rlibLoc)  ] ,"\n"))
}

if (!length(path)){
  if (file.exists("Chicago") & file.info("Chicago")$isdir){
    message("Found uncompressed Chicago folder in the current location.")
    loc = "Chicago"
  }else{
    chicagoFiles = list.files(pattern = "Chicago")
    chicagoTarGz = grep("\\.tar\\.gz$", chicagoFiles)
    if(length(chicagoTarGz == 1)){
      message("Found compressed Chicago package in the current location.")
      loc = chicagoFiles[chicagoTarGz]
    }else{
      stop("Could not unambiguously locate Chicago package, please provide as argument.\n")
    }
  
  }
}else{
  if (file.exists(path)){
    loc = path
  }else{
    stop("Could not locate Chicago R package. Check the --path argument.\n")
  }
}

if (length(rlib)){
  if (!file.exists(rlib) | !file.info(rlib)$isdir){
    stop("The specified R library directory does not exist.\n")
  }else{
    message("Installing Chicago R package and dependencies to directory: ", rlib)
  }
}else{
  message("Installing Chicago R package and dependencies to the default directory ", .libPaths()[1])
}

if (length(bin)){
  if (!file.exists(bin) | !file.info(bin)$isdir){
    message("Creating chicagoTools target directory at ", bin)
    if (!dir.create(file.path(bin,"chicagoTools"))){
      stop("Could not create chicagoTools target directory.\n")
    }
  }else{
    message("Installing chicagoTools to directory: ", bin)
  }
}else{
  message("chicagoTools will be retained at their original location")
}

message("\nInstalling dependencies if needed...\n")

if(!"MASS" %in% rownames(installed.packages())){
  message("Installing MASS...")
  install.packages(pkgs = "MASS", lib=rlib, repos = "http://cran.rstudio.com")  
}

instDT = 0
if(! "data.table" %in% rownames(installed.packages())){ 
  instDT = 1
}else{
  if (compareVersion('1.9.4', as.character(packageVersion("data.table")))==1){
    instDT = 1
  }
}
if (instDT){
  message("Installing data.table...")
  install.packages(pkgs = "data.table", lib=rlib, repos = "http://cran.rstudio.com")  
}

instDelap = 0
if(! "Delaporte" %in% rownames(installed.packages())){
  instDelap = 1
}else{
  if (compareVersion('2.2.0', as.character(packageVersion("Delaporte")))==1){
    instDelap = 1
  }
}
if (instDelap){  
  message("Installing Delaporte...")
  install.packages(pkgs = "Delaporte", lib=rlib, repos = "http://cran.rstudio.com")  
}

if(! "Hmisc" %in% rownames(installed.packages())){
  message("Installing Hmisc...")
  install.packages(pkgs = "Hmisc", lib=rlib, repos = "http://cran.rstudio.com")  
}

if(! "argparser" %in% rownames(installed.packages())){
  message("Installing argparser...")
  install.packages(pkgs = "argparser", lib=rlib, repos = "http://cran.rstudio.com")  
}

if(! "matrixStats" %in% rownames(installed.packages())){
  message("Installing matrixStats...")
  install.packages(pkgs = "matrixStats", lib=rlib, repos = "http://cran.rstudio.com")  
}

message("\nInstalling Chicago R package...\n")
install.packages(pkgs = loc, repos=NULL, lib=rlib)

#TODO: install data package

if(!is.null(bin)){
  message("\nInstalling chicagoTools...\n")
  binct = file.path(bin, "chicagoTools")
  for (f in list.files(path="chicagoTools", pattern = ".sh$")){
    file.copy(f, binct)
    Sys.chmod(file.path(binct, f), "770")
  }
  for (f in list.files(path="chicagoTools", pattern = "(.R$)|(.py$)")){
    file.copy(f, binct)
    #Sys.chmod(file.path(binct, f), "770")
  }
}

