library(survival)
library(io)
library(ggplot2)
library(tidyr)
source("common.R")

set.seed(1337)

N <- 1000;
alpha <- 0.5;      # shape
lambda <- 0.5;     # rate

out.fn <- filename("event-prob", tag=c(sprintf("shape-%.0f", alpha), sprintf("rate-%.0f", lambda)));
pdf.fn <- tag(out.fn, ext="pdf");

# Weibull distribution

# pdf
f_weibull <- function(x, alpha, lambda) {
	lp <- log(alpha) + log(lambda) + (alpha - 1)*log(x) - lambda * x^alpha;
	exp(lp)
}

# hazard rate
h_weibull <- function(x, alpha, lambda)  {
	lp <- log(alpha) + log(lambda) + (alpha - 1)*log(x);
	exp(lp)
}

# under this parameterization of Weibull,
# modifying the hazard by a multiplicative factor is
# equivalent to modifying the rate parameter lambda

f <- function(x) f_weibull(x, alpha, lambda);
S <- function(x) S_weibull(x, alpha, lambda);
h <- function(x) h_weibull(x, alpha, lambda);

# factor distributions
theta <- c(0.5, 0.5, 0.5);

X <- matrix(
	unlist(lapply(theta,
		function(th) {
			rbinom(N, 1, prob=th)
		}
	)),
	nrow=N
);
colnames(X) <- c("null", "neg", "pos");
# null: no effect
# neg: protective factor
# pos: risk factor

# effects
beta <- matrix(c(0, -1, 1), nrow=3);
#beta <- matrix(c(0, 0, 1), nrow=3);
#beta <- matrix(c(0, -1, 0), nrow=3);

# log overall effect
z <- log(lambda) + X %*% beta;

# times to event
et <- r_weibull(N, alpha, exp(z));


cr.max <- quantile(et, 0.95);

# right-censor times
cr <- runif(N, 0, cr.max);
#cr <- pmax(0, rnorm(N, cr.max * 0.5, cr.max * 0.5 * 0.25));

# observed times
ot <- pmin(et, cr);

# event status (0: censored, 1: event)
s <- as.integer(et <= cr);


# unique observed event times
tt <- sort(unique(ot[s == 1]));

tau <- quantile(et, 0.75);
cl <- classify_tte(ot, s, tau);

table(s, cl, useNA="always")

mean(s)
mean(cl, na.rm=TRUE)
1 - S(tau)

# probability of event at tau
pe.tau <- 1 - S(tt);

pe.tau.hat <- unlist(lapply(tt,
	function(tau) {
		mean(classify_tte(ot, s, tau), na.rm=TRUE)
	}
));

r <- cor(pe.tau, pe.tau.hat);
r.label <- sprintf("italic(r)^2 == %0.3f", r^2);
r.annot <- annotate("text", x=-Inf, y=Inf, label=r.label, parse=TRUE, hjust=-0.5, vjust=2);

tau.df <- data.frame(
	tau = tt,
	estimator = pe.tau.hat,
	truth = pe.tau
);
tau.ldf <- pivot_longer(tau.df, cols=c("estimator", "truth"));

qdraw(
	ggplot(tau.df, aes(x=estimator, y=truth)) +
		theme_classic() +
		geom_abline(colour="grey60", linetype=3) +
		geom_point(alpha=0.1) +
		r.annot
	,
	width = 4, height = 4,
	file = tag(pdf.fn, tag="cor")
)

qdraw(
	ggplot(tau.ldf, aes(x=tau, y=value, linetype=name)) +
		theme_classic() +
		geom_line() +
		labs(linetype="") +
		xlab("time") + ylab("event probability") +
		theme(legend.position = c(0.75, 0.25)) +
		r.annot
	,
	width = 4, height = 4,
	file = tag(pdf.fn, tag="line")
)


d <- data.frame(
	time = ot,
	status = s,
	class = cl,
	X
);

# check that the simulated data shows the correct pattern
survdiff(Surv(time, status) ~ null, data = d)
survdiff(Surv(time, status) ~ neg, data = d)
survdiff(Surv(time, status) ~ pos, data = d)
coxph(Surv(time, status) ~ null + neg + pos, data=d)

summary(glm(status ~ null + neg + pos, data=d))
summary(glm(class ~ null + neg + pos, data=d))

