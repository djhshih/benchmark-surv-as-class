
# generate data follwing weibull distribution 
r_weibull <- function(n, alpha, lambda) {
  # R implements weibull(x; a, b)
  # where a = alpha, b = lambda^(-1/alpha)
  b <- lambda^(-1/alpha);
  
  rweibull(n, shape=alpha, scale=b)
}

# survival function
S_weibull <- function(x, alpha, lambda) {
  lp <- -lambda * x^alpha;
  exp(lp)
}

classify_tte <- function(ftime, fstatus, t.cut) {
  ifelse(ftime > t.cut,
         0L,
         ifelse(fstatus, 1L, NA)
  )
}