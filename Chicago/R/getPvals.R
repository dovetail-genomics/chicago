getPvals <-
function(cd){
  ## - Calls p-values
  
  # No need for this anymore
  alpha = cd@params$dispersion
  x = cd@x
  if(is.null(alpha)) {stop("getPvals: 'dispersion' parameter of x not found.")}
  
  message("Calculating p-values...") 
  
  ##p-values:
  ##(gives P(X > x-1) = P(X >= x))
  ##The "ifelse" is because pdelap cannot deal with beta=0.
  ##TODO can probably optimize this:
  #   x[,"log.p"] <- ifelse(
  #     x$Bmean < .Machine$double.eps,
  #     ppois(x[,Ncol] - 1L, lambda=x$Tmean, lower.tail=FALSE, log.p=TRUE),
  #     pdelap(x[,Ncol] - 1L, alpha, beta=x$Bmean/alpha, lambda=x$Tmean, lower.tail=FALSE, log.p=TRUE)
  #   )
  x[,log.p:= 
      ifelse(Bmean < .Machine$double.eps,
             ppois(N - 1L, lambda=Tmean, lower.tail=FALSE, log.p=TRUE),
             pdelap(N - 1L, alpha, beta=Bmean/alpha, lambda=Tmean, lower.tail=FALSE, log.p=TRUE)
      )
    ]
  
  # Large N approximation ---------------------------------------------------
  
  ##In rare cases where pdelap returns Infs, estimate the p-value magnitude
  ##using an NB approximation, through method of moments argument
  ##NaNs occur when pdelap() thinks the p-value is negative (since can have 1 - 1 != 0),
  ##thus these are also approximated.
  sel <- which(is.infinite(x$log.p) | is.nan(x$log.p))
  if(length(sel) > 0)
  {
    message("Approximating ", length(sel), " very small p-values.")
    
    gamma <- x[sel,alpha*(1+Tmean/Bmean)^2] ##gamma is the "effective" dispersion
    
    ##in the case where Bmean << Tmean, gamma becomes Inf, so cap gamma above
    ##(should make very little difference since we are basically Poisson in this case)
    ##Case where Bmean >> Tmean causes no problem.
    gamma <- pmin(gamma, 1e10)
    
    x[sel,log.p := pnbinom(N - 1L, size=gamma, mu=Bmean+Tmean, lower.tail=FALSE, log.p=TRUE)]
    
    if(any(is.infinite(x[sel,log.p]))) {warning("Some log-p-values were infinite.")}
    if(any(is.nan(x[sel,log.p]))) {warning("Some log-p-values were NaNs.")}
  }
  if(any(is.na(x$log.p))) {warning("Some log-p-values were NA.")}
  cd
}
