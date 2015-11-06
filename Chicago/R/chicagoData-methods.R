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
  validObject(cd)
  
  x <- cd@settings
    
  ##Find which settings are not the same as the defaults
  ref <- unlist(defaultSettings())
  query <- unlist(x)

  output <- rbind(query, ref)
  
  ##specifically find the columns that are not uniform
  sel <- !apply(output, 2, function(x) {
    identical(unname(x[1]), unname(x[2]))
    }
  )
  
  output <- output["query", sel]
  
  cat("Non-defaults are:\n")
  show(output)
}

##validity

##check setting names are correct
setValidity("chicagoData",
            function(object)
            {
              ref <- unlist(defaultSettings())
              query <- unlist(object@settings)
              
              if(length(ref) != length(query) | any(!(names(query) %in% names(ref))))
              {
                return("Object's setting names do not match defaultSettings().")
              }
              TRUE
            }
          )


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
