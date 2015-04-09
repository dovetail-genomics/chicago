modifySettings <-
function(cd, settings){
  
  message("Warning: settings are not checked for consistency with the previous ones.")
  
  for (s in names(settings)){
    cd@settings[[s]] = settings[[s]]
  }
  
  cd
}
