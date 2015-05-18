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

##show

setMethod("show", "chicagoData", 
          function(object)
          {
            cat('An object of class "chicagoData"\n')
            cat("Slot 'x':\n")
            show(object@x)
            describeVec(object@params, "params")
            describeVec(object@settings, "settings")
          }
)