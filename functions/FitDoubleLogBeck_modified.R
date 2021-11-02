#fitting the doublelog beck but allowing it to leaving the blanks in there
# modified from original phenopix implementation
# https://cran.r-project.org/web/packages/phenopix/phenopix.pdf
#for missing data
# as well as adding back in the -1 as per the original beck formulation
# as well as adding in the option to manually add weights which I have now set
# as well as adding 4 more starting value for the fitting, here trying starting
# values of 0.4 and 0.6 of doy as well

FitDoubleLogBeck.4 <-function (x, t = index(x), tout = t, weighting = TRUE, return.par = FALSE, 
          plot = FALSE, hessian = FALSE, sf = quantile(x, probs = c(0.05, 
                                                                    0.95), na.rm = TRUE), 
          normalize=TRUE, wts=wts, ...) 
{
  .normalize <- function(x, sf) (x - sf[1])/(sf[2] - sf[1])
  .backnormalize <- function(x, sf) (x + sf[1]/(sf[2] - sf[1])) * 
    (sf[2] - sf[1])
  if (normalize){
  x <- .normalize(x, sf = sf)
  }
  n <- length(x)
  avg <- mean(x, na.rm = TRUE)
  mx <- max(x, na.rm = TRUE)
  mn <- min(x, na.rm = TRUE)
  ampl <- mx - mn
  .doubleLog <- function(par, t) {
    mn <- par[1]
    mx <- par[2]
    sos <- par[3]
    rsp <- par[4]
    eos <- par[5]
    rau <- par[6]
    xpred <- mn + (mx - mn) * (1/(1 + exp(-rsp * (t - sos))) + 
                                 1/(1 + exp(rau * (t - eos)))-1)
    return(xpred)
  }
  .error <- function(par, x, weights) {
    if (any(is.infinite(par))) 
      return(99999)
    if (par[1] > par[2]) 
      return(99999)
    xpred <- .doubleLog(par, t = t)
    sse <- sum((xpred - x)^2 * weights, na.rm = TRUE)
    return(sse)
  }
  if (weighting) {
    iter <- 1:2
  }
  else {
    iter <- 1
  }
  #sce swapped these
  #weights <- rep(1, length(x))
  weights=wts
  #doy <- quantile(t, c(0.25, 0.75), na.rm = TRUE)
  #sce changed this to see if it helps the fits
  doy <- quantile(t, c(0.40, 0.60), na.rm = TRUE)
  doy_a <- quantile(t, c(0.25, 0.75), na.rm = TRUE)
  
  prior <- rbind(c(mn, mx, doy[1], 0.5, doy[2], 0.5),
                 c(mn, mx, doy[2], 0.5, doy[1], 0.5),
                 c(mn - ampl/2, mx + ampl/2,  doy[1], 0.5, doy[2], 0.5),
                 c(mn - ampl/2, mx + ampl/2,  doy[2], 0.5, doy[1], 0.5),
                 c(mn, mx, doy_a[1], 0.5, doy_a[2], 0.5),
                 c(mn, mx, doy_a[2], 0.5, doy_a[1], 0.5),
                 c(mn - ampl/2, mx + ampl/2,  doy_a[1], 0.5, doy_a[2], 0.5),
                 c(mn - ampl/2, mx + ampl/2,  doy_a[2], 0.5, doy_a[1], 0.5))
  if (plot) 
    plot(t, x)
  for (i in iter) {
    opt.l <- apply(prior, 1, optim, .error, x = x, weights = weights, 
                   method = "BFGS", control = list(maxit = 1000), 
                   hessian = hessian)
    opt.df <- cbind(cost = unlist(plyr::llply(opt.l, function(opt) opt$value)), 
                    convergence = unlist(plyr::llply(opt.l, function(opt) opt$convergence)), 
                    plyr::ldply(opt.l, function(opt) opt$par))
    #find best starting values?
    best <- which.min(opt.df$cost)
    if (opt.df$convergence[best] == 1) {
      opt <- opt.l[[best]]
      opt <- optim(opt.l[[best]]$par, .error, x = x, weights = weights, 
                   method = "BFGS", control = list(maxit = 1500), 
                   hessian = hessian)
      prior <- rbind(prior, opt$par)
      xpred <- .doubleLog(opt$par, t)
    }
    #convergence of 0 means converged
    else if (opt.df$convergence[best] == 0) {
      opt <- opt.l[[best]]
      prior <- rbind(prior, opt$par)
      xpred <- .doubleLog(opt$par, t)
    }
    if (plot) {
      plyr::llply(opt.l, function(opt) {
        xpred <- .doubleLog(opt$par, t)
        lines(t, xpred, col = "cyan")
      })
      lines(t, xpred, col = "blue", lwd = 2)
    }
    parinit <- opt$par
    mn <- opt$par[1]
    mx <- opt$par[2]
    sos <- opt$par[3]
    rsp <- opt$par[4]
    eos <- opt$par[5]
    rau <- opt$par[6]
    #sce does not understand what this is
    m <- lm(c(0, 100) ~ c(sos, eos))
    tr <- coef(m)[2] * t + coef(m)[1]
    tr[tr < 0] <- 0
    tr[tr > 100] <- 100
    res <- xpred - x
    weights <- 1/((tr * res + 1)^2)
    weights[res > 0 & res <= 0.01] <- 1
    weights[res < 0] <- 4
    #sce added this to keep downweighting the infilled ones
    weights=weights*wts
  }
  if (opt$convergence != 0) {
    opt$par[] <- NA
    xpred <- rep(NA, length(tout))
  }
  else {
    xpred <- .doubleLog(opt$par, tout)
  }
  if (normalize){
  xpred <- .backnormalize(xpred, sf = sf)
  }
  xpred.out <- zoo::zoo(xpred, order.by = t)
  names(opt$par) <- c("mn", "mx", "sos", 
                      "rsp", "eos", "rau")
  
  if (hessian) {
    opt.new <- optim(opt$par, .error, x = x, weights = weights, 
                     method = "BFGS", hessian = TRUE)
    .qr.solve <- function(a, b, tol = 1e-07, LAPACK = TRUE) {
      if (!is.qr(a)) 
        a <- qr(a, tol = tol, LAPACK = LAPACK)
      nc <- ncol(a$qr)
      nr <- nrow(a$qr)
      if (a$rank != min(nc, nr)) 
        stop("singular matrix 'a' in solve")
      if (missing(b)) {
        if (nc != nr) 
          stop("only square matrices can be inverted")
        b <- diag(1, nc)
      }
      res <- qr.coef(a, b)
      res[is.na(res)] <- 0
      res
    }
    vc <- .qr.solve(opt$hessian)
    npar <- nrow(vc)
    s2 <- opt.df$cost[best]^2/(n - npar)
    std.errors <- sqrt(diag(vc) * s2)
  }
  fit.formula <- expression(mn + (mx - mn) * (1/(1 + exp(-rsp * 
                                                           (t - sos))) + 1/(1 + exp(rau * (t - eos))))-1)
  output <- list(predicted = xpred.out, params = opt$par, formula = fit.formula, 
                 sf = sf)
  if (hessian) 
    output <- list(predicted = xpred.out, params = opt$par, 
                   formula = fit.formula, stdError = std.errors, sf = sf)
  return(output)
}
