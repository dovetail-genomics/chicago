##convenience functions

describeMx <- function(x, nm)
{
  if(is.null(x)) {cat(nm, ": NULL", sep = "")} 
  else cat(nm, ": ", nrow(x), " rows, ", ncol(x), " columns.\n", sep = "")
}

describeVec <- function(x, nm)
{
  if(is.null(x)) {cat(nm, ": NULL", sep = "")} 
  else cat(nm, ": length ", length(x), ".\n", sep = "")
}

describeVecNames <- function(x, nm)
{
  if(is.null(x)) {cat(nm, ": NULL", sep = "")}
  else cat(nm, ": ", paste(names(x), collapse=", "), ".\n", sep = "")
}


describeNonDefaultSettings <- function(cd)
{
  ##FIXME If we ever introduce an argument that's a list of things, this function will need to be changed.
  #validObject(cd)
  
  x <- cd@settings
    
  ##Find which settings are not the same as the defaults
  ref <- unlist(defaultSettings())
  query <- unlist(x)

  rnames <- names(ref)
  qnames <- names(query)
  
  ##find missing/extraneous settings
  sel.missing <- rnames[!rnames %in% qnames]
  sel.extraneous <- qnames[!qnames %in% rnames]

  if(length(sel.missing) > 0)
  { 
    msg <- paste("Missing setting(s): ", paste(sel.missing, collapse= ", "))
    if ("brownianNoise.samples" %in% sel.missing){
      msg <- paste0(msg, "\n(The missing brownianNoise.samples setting likely indicates that it was equal to 1 in the original analysis).")
    }
    warning(msg)
  }

  if(length(sel.extraneous) > 0)
  {
    warning("Extraneous setting(s): ", paste(sel.extraneous, collapse= ", "))
  }
  
  ##specifically find the columns that are not uniform
  sel.shared <- rnames[rnames %in% qnames]

  sel.changed <- sapply(sel.shared, function(x) {
    !identical(ref[x], query[x])
    }
  )
  
  output <- query[sel.changed]
  
  if(any(sel.changed))
  {
    cat("Non-defaults are:\n")
    show(output)
  }else{
    cat("No non-default settings.")
  }
  
  ##Find and list missing settings
  ##FIXME
}

##validity

##check setting names are correct
##Have removed this, to avoid confusion when new settings are added in later versions of CHiCAGO.
##Now WARN in the show() method.
# setValidity("chicagoData",
#             function(object)
#             {
#               ref <- unlist(defaultSettings())
#               query <- unlist(object@settings)
#               
#               ##FIXME allow names to be *missing* from cd
#               if(length(ref) != length(query) | any(!(names(query) %in% names(ref))))
#               {
#                 return("Object's setting names do not match defaultSettings().")
#               }
#               TRUE
#             }
#           )


##show

setMethod("show", "chicagoData", 
          function(object)
          {
            cat('An object of class "chicagoData"\n')
            cat("Slot 'x':\n")
            show(object@x)
            describeVecNames(object@params, 'Slot "params"')
            describeVec(object@settings, 'Slot "settings"')
            describeNonDefaultSettings(object)
          }
)

##accessors

setMethod("intData", "chicagoData", function(x){x@x})
setMethod("settings", "chicagoData", function(x){x@settings})
setMethod("params", "chicagoData", function(x){x@params})

##replacers

setReplaceMethod("intData", c(x="chicagoData"), 
                 function(x, value)
                 {
                   x@x <- value
                   x
                 }
)

setReplaceMethod("params", c(x="chicagoData"), 
                 function(x, value)
                 {
                   x@params <- value
                   x
                 }
)
