library(survival)
library(io)
library(ggplot2)
library(tidyr)
library(dplyr)
library(reshape2)
source("common.R")


N <- 1000;
alpha_v <- c(0.5, 1, 1.5, 2); # shape
lambda_v <- c(0.5, 1, 1.5, 2); # rate
quantile.tau = 0.75; # quantitle for tau
theta <- c(0.5, 0.5, 0.5); #factor distributions to generate X

simu_cox_glm <- function(N, alpha_v, lambda_v, quantile.tau, theta, niter){
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
  
  # effects
  beta <- matrix(c(0, -1, 1), nrow=3); # null, negative, positive
  
  cox.est.list <- list()
  glm.est.list <- list()
  
  for(i in seq(niter)){
    X <- matrix(
      unlist(lapply(theta,
                    function(th) {
                      rbinom(N, 1, prob=th)
                    }
      )),
      nrow=N
    );
    colnames(X) <- c("null", "neg", "pos");
    
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
    colnames(ot) <- group;
    
    # event status (0: censored, 1: event)
    s <- (et <= cr);
    mode(s) <- "integer";
    s <- as.data.frame(s, col.names = group);
    colnames(s) <- group
    
    tau <- lapply(et, function(x) quantile(x, quantile.tau));
    
    cl <- mapply(function(time, status, cut) classify_tte(time, status, cut), 
                 ot, s, tau);
    
    ot.list <- as.list(ot);
    s.list <- as.list(s);
    cl.list <- as.list(as.data.frame(cl));
    
    d <- mapply(function(time, status, class) data.frame(time, status, class, X), 
                ot.list, s.list, cl.list, SIMPLIFY = FALSE);
    
    
    # fit the cox PH model
    cox.fit <- lapply(d, function(x) coxph(Surv(time, status) ~ null + neg + pos, data = x));
    
    cox.est.m <- lapply(cox.fit, function(x) 
      melt(cbind(
          est = x$coefficients, 
          p_value = coef(summary(x))[, 5])
      )
    );
    
    cox.est <- lapply(cox.est.m, function(x) 
      t(data.frame(
        value = x$value,
        row.names = paste(x$Var1, x$Var2, sep = ".")
      ))
    );
    
    cox.est.list[[i]]  <- do.call(rbind, Map(data.frame, cox.est));
    cox.est.list[[i]]$group <- group
   
    # fit the glm model
    glm.fit <- lapply(d, function(x) glm(class ~ null + neg + pos, data = x, family=binomial(logit)));
    
    glm.est.m <- lapply(glm.fit, function(x) 
      melt(
        cbind(
        est = x$coefficients, 
        p_value = coef(summary(x))[, 4]
        )[-1, ]
        )
      );
    
    glm.est <- lapply(glm.est.m, function(x) 
      t(data.frame(
        value = x$value,
        row.names = paste(x$Var1, x$Var2, sep = ".")
        ))
      );
    
    
    glm.est.list[[i]]  <- do.call(rbind, Map(data.frame, glm.est));
    
    glm.est.list[[i]]$group <- group
  }
  
  cox.est.df <- do.call("rbind", cox.est.list)
  cox.est.df$model <- "survival";
  cox.est.df$null_p_idx <- ifelse(cox.est.df$null.p_value > 0.05, 1, 0)
  cox.est.df$neg_est_idx <- ifelse(cox.est.df$neg.est < 0, 1, 0)
  cox.est.df$neg_p_idx <- ifelse(cox.est.df$neg.p_value < 0.05, 1, 0)
  cox.est.df$pos_est_idx <- ifelse(cox.est.df$pos.est > 0, 1, 0)
  cox.est.df$pos_p_idx <- ifelse(cox.est.df$pos.p_value < 0.05, 1, 0)
  
  glm.est.df <- do.call("rbind", glm.est.list)
  glm.est.df$model <- "classification";
  glm.est.df$null_p_idx <- ifelse(glm.est.df$null.p_value > 0.05, 1, 0)
  glm.est.df$neg_est_idx <- ifelse(glm.est.df$neg.est < 0, 1, 0)
  glm.est.df$neg_p_idx <- ifelse(glm.est.df$neg.p_value < 0.05, 1, 0)
  glm.est.df$pos_est_idx <- ifelse(glm.est.df$pos.est > 0, 1, 0)
  glm.est.df$pos_p_idx <- ifelse(glm.est.df$pos.p_value < 0.05, 1, 0)
  
  plot.d <- rbind(cox.est.df, glm.est.df);
  
  #binom_se <- function(p, n, z=1.96) {
  #  z * sqrt(p * (1- p) / n)
  #}
  
  accu.est.d <- plot.d %>% group_by(group, model) %>% 
    summarise(
      `negative` = mean(neg_est_idx),
      `positive` = mean(pos_est_idx),
      );
  
 # negative_se <- binom_se(mean(neg_est_idx), niter)
  
  accu.p.d <- plot.d %>% group_by(group, model) %>% 
    summarise(
      `null` = mean(null_p_idx),
      `negative` = mean(neg_p_idx), 
      `positive` = mean(pos_p_idx)
    );

  accu.est.d.m <- melt(accu.est.d)
  
  accu.p.d.m <- melt(accu.p.d)
  
  
  qdraw(
    ggplot(accu.est.d.m, aes(x = variable, y = value, fill = model)) +
      facet_wrap(~ group, ncol = 4) +
      theme_classic() +
      geom_bar(stat="identity", position = "dodge", width=0.5) + 
      labs(title = "Risk factor identification based on coefficient sign",
           x = "risk factor",
           y = "accuracy")
    , width = 8, height = 8, "risk-factor-id_sign.pdf"
  )
  
  qdraw(
    ggplot(accu.p.d.m, aes(x = variable, y = value, fill = model)) +
      facet_wrap(~ group, ncol = 4) +
      theme_classic() +
      geom_bar(stat="identity", position = "dodge", width=0.5) +
      labs(title = "Risk factor identification based on statistical significance",
           x = "risk factor",
           y = "accuracy")
    , width = 8, height = 8, "risk-factor-id_p-value.pdf"
  )
}

simu_cox_glm(N, alpha_v = alpha_v, lambda_v = lambda_v, 
             quantile.tau = quantile.tau, theta = theta, niter = 100);