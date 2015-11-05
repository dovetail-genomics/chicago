##Accessors

setGeneric("intData", function(x) standardGeneric("intData"))
setGeneric("params", function(x) standardGeneric("params"))
setGeneric("settings", function(x) standardGeneric("settings"))

##Replacers

setGeneric("intData<-", function(x, ...) standardGeneric("intData<-")) 
setGeneric("params<-", function(x, ...) standardGeneric("params<-")) 
#setGeneric("settings<-", function(x, ...) standardGeneric("settings<-")) ##use modifySettings() instead
