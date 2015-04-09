estimateDistFun <-
function (cd, method="cubic", n.obs.head=10, n.obs.tail=25, logScale=FALSE, outfile=NULL) {
  
  # Take the "refBinMean" column of the data x as f(d_b)
  # then interpolate & extrapolate to get f(d).
  # TODO output extra diagnostic information?
  # TODO optimize with data.table
  
  if (!method %in% c("lm", "cubic")){
    stop ("Unknown method.\n")
  }
  
  # Get f(d_b)
#   f.d <- unique(x[!is.na(x$refBinMean),c("distbin", "refBinMean")]) ##delete rows with NAs from baits that are too far away
  setkey(cd@x, distbin, refBinMean)
  f.d <- unique(cd@x)[is.na(refBinMean)==FALSE][, c("distbin", "refBinMean"), with=F]
  
#   f.d <- f.d[order(f.d$refBinMean, decreasing=TRUE),]
  f.d <- f.d[order(refBinMean, decreasing=TRUE)]

  setDF(f.d) # f.d is tiny, so no need to bother with it being a data.table
  f.d$midpoint <- seq(from=round(cd@settings$binsize/2), by=cd@settings$binsize, length.out=nrow(f.d))
  
  obs.min <- log(min(f.d$midpoint))
  obs.max <- log(max(f.d$midpoint))
  
  if(method == "lm") {
    ##On log-scale, do a linear interpolation.
    ##Linear models applied to first (n.obs.head) observations, and to last (n.obs.tail) observations.
    
    ##Interpolation: Estimate f(d) (NB "rule" parameter = 1 forces NAs outside of range)
    log.f.obs <- approxfun(log(f.d$midpoint), log(f.d$refBinMean), rule=c(1,1))
    
    ##Extrapolation: Fit the "head" and "tail" of f using a linear model
    head.coef <- coefficients(lm(log(refBinMean)~log(midpoint), data = head(f.d, n.obs.head))) ##Fit for small d
    tail.coef <- coefficients(lm(log(refBinMean)~log(midpoint), data = tail(f.d, n.obs.tail))) ##Fit for large d
    
    log.f.head <- function(x, head.coef) {head.coef[1] + x*head.coef[2]}
    log.f.tail <- function(x, tail.coef) {tail.coef[1] + x*tail.coef[2]} ##explicitly stated in case of later change
    
  }
  
  if(method == "cubic") {
    ##Spline - Cubic fit over observed interval, linear fit elsewhere, assume continuity of f(d) & f'(d).
    
    ##cubic fit (quadratic not immensely different TBH)
    f.d.cubic <- lm(log(refBinMean) ~ log(midpoint) + I(log(midpoint)^2) + I(log(midpoint)^3), data = f.d)
    fit <- f.d.cubic$coefficients
    
    ##Interpolation: Estimate f(d) from cubic (NB see "rule" parameter for what to do outside range)
    log.f.obs <- function(x, fit. = fit) {fit.[1] + fit.[2]*x + fit.[3]*(x^2) + fit.[4]*(x^3)}
    
    ##Extrapolation: Fit the "head" and "tail" of f using continuity
    obs.min <- log(min(f.d$midpoint))
    obs.max <- log(max(f.d$midpoint))
    
    beta <- fit[2] + 2*fit[3]*c(obs.min, obs.max) + 3*fit[4]*(c(obs.min, obs.max)^2)
    alpha <- fit[1] + (fit[2] - beta)*c(obs.min, obs.max) + fit[3]*c(obs.min, obs.max)^2 + fit[4]*c(obs.min, obs.max)^3
    
    head.coef <- c(alpha[1], beta[1])
    tail.coef <- c(alpha[2], beta[2])
    
    log.f.head <- function(x, head.coef.=head.coef) {head.coef.[1] + x*head.coef.[2]}
    log.f.tail <- function(x, tail.coef.=tail.coef) {tail.coef.[1] + x*tail.coef.[2]} ##explicitly stated in case of later change
    
  }
  
  ##Put everything together to get the final function
  ##All these appended dots (e.g. "head.conf.") are to avoid errors of form "promise already under evaluation"
  log.f <- function(x, head.coef.=head.coef, tail.coef.=tail.coef, obs.min.=obs.min, obs.max.=obs.max)
  {
    ifelse(x > obs.max.,
           log.f.tail(x, tail.coef.), ##Common case evaluated first
           ifelse(x < obs.min.,
                  log.f.head(x, head.coef.),
                  log.f.obs(x))
    )
  }
  if(logScale)
  {
    f <- log.f
  } else {
    f <- function(x) exp(log.f(log(x)))
  }

  if (!is.null(outfile)){ 
    pdf(outfile)
  }
    curve(log.f.obs, obs.min, obs.max,
          main = paste0("Distance function (points = obs, line = ", method, " fit)"),
          xlab = "log(distance)",
          ylab = "log(f(d))")
    with(f.d, points(log(midpoint), log(refBinMean)))
  if (!is.null(outfile)){ 
    dev.off()
  }
  
  cd@params$f = f
  cd
}
