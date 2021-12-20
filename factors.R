library(survival)
library(io)
library(ggplot2)
library(tidyr)

#seed <- 1334;

N <- 1000;
nu <- 3;      # shape
lambda <- 2;  # scale

out.fn <- filename("event-prob", tag=c(sprintf("shape-%.0f", nu), sprintf("scale-%.0f", lambda)));
pdf.fn <- tag(out.fn, ext="pdf");

# survival function
S <- function(x, ...) pweibull(x, shape=nu, scale=lambda, lower.tail=FALSE, ...);

# factor distributions
theta <- c(0.25, 0.5, 0.5);

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
beta <- matrix(c(0, -1, 1), nrow=3) * 0.5;

# log overall effect
z <- log(lambda) - X %*% beta;

# times to event
et <- rweibull(N, shape=nu, scale=exp(z));

cr.max <- quantile(et, 0.95);

# right-censor times
#cr <- runif(N, 0, cr.max);
cr <- pmax(0, rnorm(N, cr.max * 0.5, cr.max * 0.5 * 0.25));

# observed times
ot <- pmin(et, cr);

# event status (0: censored, 1: event)
s <- as.integer(et <= cr);


# unique observed event times
tt <- sort(unique(ot[s == 1]));

classify_tte <- function(ftime, fstatus, t.cut) {
	ifelse(ftime > t.cut,
		0L,
		ifelse(fstatus, 1L, NA)
	)
}

tau <- 0.5;
cl <- classify_tte(ot, s, tau);

table(s, cl, useNA="always")

mean(s)
mean(cl, na.rm=TRUE)
1 - S(tau)

# probability of event at tau
pe.tau <- 1 - S(tt);

pe.tau.hat <- unlist(lapply(tt,
	function(tau) {
		mean(surv_to_class(ot, s, tau), na.rm=TRUE)
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
coxph(Surv(time, status) ~ ., data=d)

summary(glm(status ~ . - time - class, data=d))
summary(glm(class ~ . - time - status, data=d))

