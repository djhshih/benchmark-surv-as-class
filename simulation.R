library(survival)
library(io)
library(ggplot2)
library(tidyr)
library(reshape2)
source("common.R")

set.seed(1337);

N <- 1000;
alpha_v <- c(0.5, 1, 1.5, 2); # shape
lambda_v <- c(0.5, 1, 1.5, 2); # rate

# create group names
group <- levels(
  interaction(
    paste("alpha=", alpha_v, sep = ""),
    paste("lambda=", lambda_v, sep = ""),
    sep=", "
  )
);

parameters <- as.data.frame(cbind(sort(rep(lambda_v, length(lambda_v))), alpha_v))
colnames(parameters) <- c("lambda", "alpha")
rownames(parameters) <- group;
parameters <- parameters[, c(2, 1)]

theta <- c(0.5, 0.5, 0.5); #factor distributions
X <- matrix(
  unlist(lapply(theta,
                function(th) {
                  rbinom(N, 1, prob=th)
                }
  )),
  nrow=N
);
colnames(X) <- c("null", "neg", "pos");

# effects
beta <- matrix(c(0, -1, 1), nrow=3);

# log overall effect
z_v <- lapply(lambda_v, function(x) log(x) + X %*% beta);


# times to event
et<- lapply(z_v, function(x) 
  lapply(alpha_v, 
         function(y) r_weibull(N, y, exp(x))
         )
  );
et <- unlist(et, recursive=FALSE); 

# right-censor times
# uniform distribution
cr.max <- as.numeric(lapply(et, function(x) quantile(x, 0.95)));
cr <- lapply(cr.max, function(x) runif(N, 0, x));

et <- as.data.frame(et, col.names = group);
cr <- as.data.frame(cr, col.names = group);

# observed times
ot <- pmin(et, cr);

# event status (0: censored, 1: event)
s <- (et <= cr);
mode(s) <- "integer";
s <- as.data.frame(s, col.names = group);

# observed event times
tt <- seq(0, 5, by=0.1)

# event probability; truth
pe.tau <- mapply(
  function(alpha, lambda) 1 - S_weibull(tt, alpha, lambda),
  parameters$alpha,
  parameters$lambda
);
colnames(pe.tau) <- group

# event probability; estimator
pe.tau.hat <- lapply(tt, function(x) 
  mapply(
    function(ot, s) mean(classify_tte(ot, s, x), na.rm=TRUE),
    ot,
    s
    )
)
pe.tau.hat <- matrix(unlist(pe.tau.hat), nrow = length(tt), byrow = TRUE)
colnames(pe.tau.hat) <- group


pe.tau.m <- melt(pe.tau, varnames=c("idx", "group"), value.name = "truth")
pe.tau.hat.m <- melt(pe.tau.hat)

d <- cbind(pe.tau.m, estimator = pe.tau.hat.m$value, time = tt)
d.ldf <- pivot_longer(d, cols = c("truth", "estimator"))

ggplot(d.ldf, aes(x=time, y=value, linetype=name)) +
  facet_wrap(~ group, ncol = 4) +
  theme_classic() +
  geom_line() +
  labs(linetype="") +
  xlim(0, 5) +
  xlab("time") + ylab("event probability") +
  theme(legend.position = "right")
