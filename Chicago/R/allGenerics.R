##Accessors

setGeneric("intData", function(x) standardGeneric("intData"))
setGeneric("params", function(x) standardGeneric("params"))
setGeneric("settings", function(x) standardGeneric("settings"))

##Replacers

setGeneric("intData<-", function(x, value) standardGeneric("intData<-")) 
setGeneric("params<-", function(x, value) standardGeneric("params<-")) 
#setGeneric("settings<-", function(x, value) standardGeneric("settings<-")) ##use modifySettings() instead
