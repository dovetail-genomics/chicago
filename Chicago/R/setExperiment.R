setExperiment <-
function(designDir="", settings=list(), settingsFile=NULL,  
 def.settings=list(
  rmapfile= file.path(designDir, "Digest_Human_HindIII.bed"),
  baitmapfile= file.path(designDir, "Digest_Human_HindIII_baits_ID.bed"),
  nperbinfile = file.path(designDir, "Digest_Human_HindIII_NperBin.txt"),
  nbaitsperbinfile = file.path(designDir, "Digest_Human_HindIII_NbaitsPerBin.txt"),
  proxOEfile = file.path(designDir, "proxOE_out.txt"),
  Ncol = "N",
  baitmapFragIDcol=4,
  baitmapGeneIDcol=5,
  maxLBrownEst = 1.5e6,
  minFragLen = 150,
  maxFragLen = 40000,
  minNPerBait = 250,
  binsize=20000,
  removeAdjacent = TRUE,
  adjBait2bait=TRUE,
  tlb.filterTopPercent=0.01, 
  tlb.minProxOEPerBin=1000, 
  tlb.minProxB2BPerBin=100,
  techNoise.minBaitsPerBin=1000, 
  brownianNoise.subset=1000,
  brownianNoise.seed=NULL,
  baitIDcol = "baitID",
  otherEndIDcol = "otherEndID",
  otherEndLencol = "otherEndLen", 
  distcol = "distSign",
  weightAlpha = 34.1157346557331, ##from macrophage. Remove as default?
  weightBeta = -2.58688050486759,
  weightGamma = -17.1347845819659,
  weightDelta = -7.07609245521541
  )){
  
  modSettings = vector("list")
  
  if(!is.null(settingsFile)){
    
    message(paste0("Reading custom experimental settings from ", settingsFile, "..."))
    
    # http://stackoverflow.com/questions/6602881/text-file-to-list-in-r
    sf <- scan(settingsFile, what="", sep="\n")
    modSettings <- strsplit(sf, "[[:space:]]+")
    names(modSettings) <- sapply(modSettings, `[[`, 1)
    modSettings <- lapply(modSettings, `[`, -1)
    
    # convert numerically defined settings to numbers
    # do the same for logical settings
    # suppressing "NAs introduced by coercion" warnings
    suppressWarnings({
      modSettings <- lapply(modSettings, function(s){
        num_s = as.numeric(s);
        if(!is.na(num_s)) { return (num_s) };
        bool_s = as.logical(s);
        if(!is.na(bool_s)) { return(bool_s) };
        s 
      })
    })
  }

  for (s in names(settings)){
    modSettings[[s]] = settings[[s]]
  }
  
  
  for (s in names(modSettings)){
      def.settings[[s]] = modSettings[[s]]
  }
  
  cd = chicagoData(x=data.table(), params=list(), settings=def.settings)
}
