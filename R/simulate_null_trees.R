library(tidyverse)
library(diversitree)

## data simulations
d <- read.csv("output/summary_tree_results.csv")

## calculate a mean of lambda and mu from old trees
tree_lambda <- d %>% filter(tree.max.age > 150) %>% 
  summarise(mean(mean.clade.lambda)) %>% as.numeric()

tree_mu <- d %>% filter(tree.max.age > 150) %>% 
  summarise(mean(mean.clade.mu)) %>% as.numeric()

pars <- c(tree_lambda, tree_mu)

## do n independent simulations and estimate slope for each of these
tree.bd.mintax <- function(pars, max.t, min.taxa){
  while(1){
    tree <- tree.bd(pars, max.t=max.t)
    if (is.null(tree)){
      ntip <- 0
    } else {
      ntip <- Ntip(tree)
    }
    if (ntip > min.taxa) break()
  }
  tree
}

tree_bd_slope <- function(pars, ages, min.taxa){
  trees <- lapply(ages, function(x) tree.bd.mintax(pars, x, min.taxa))
  lambdas <- vector(length=length(ages))
  for (i in 1:length(ages)){
    lik <- make.bd(trees[[i]])
    fit <- find.mle(lik, x.init=pars)
    lambdas[i] <- fit$par["lambda"]
  }
  fit_lm <- lm(log(lambdas)~log(ages))
  as.numeric(fit_lm$coefficients["log(ages)"])
}

## filter so only look at ages less than 200 mya
d_200 <- filter(d, tree.max.age < 200)
ages <- d_200$tree.max.age
slopes <- sapply(c(1:1000), function(x) tree_bd_slope(pars, ages, 6))
saveRDS(slopes,file="output/tree_bd_slope.rds")

## estimate empirical slope
emp <- lm(log(d_200$mean.clade.lambda)~log(d_200$tree.max.age))
emp <- as.numeric(emp$coefficients[2])

## plot distribution of recovered slopes (due to sampling alone) and include empircal slope
df <- data.frame(slopes=slopes)
ggplot(df, aes(x=slopes)) + geom_histogram() + theme_bw() + xlab("slope estimates") + geom_vline(xintercept=emp)

## repeat with mu = 0
pars_pb <- c(pars[1], 0)
slopes_pb <- sapply(c(1:1000), function(x) tree_bd_slope(pars_pb, ages, 6))
df_pb <- data.frame(slopes=slopes_pb)
ggplot(df_pb, aes(x=slopes)) + geom_histogram() + theme_bw() + xlab("slope estimates") + geom_vline(xintercept=emp)
saveRDS(slopes,file="output/tree_bd_slope_mu0.rds")






## Now fixing N and looking at distribution of ages
tree.bd.time <- function(pars, max.taxa, n){
  out <- vector(length=n)
  for (i in 1:n){
    tree <- tree.bd(pars, max.taxa=taxa)
    if (is.null(tree)){
      out[i] <- 0
    } else {
      out[i] <- max(branching.times(tree))
    }
  }
}

n_t <- lapply(d$n.clade, function(x) 
  tree.bd.time(pars=c(tree_lambda, tree_mu), max.taxa=x, n=1000))
saveRDS(n_t,file="output/n_t.rds")

## calculate where in the distribution of possible branching times does our tree land
perc_t <- vector(length=nrow(d))
for (j in 1:nrow(d)){
  perc_t[j] <- length(which(n_t[[j]] < d$tree.max.age[j]))
}

## for every age simulate 1000 trees to get a range of possible N
ages <- d$tree.max.age

tree.bd.n <- function(pars, t, n){
  out <- vector(length=n)
  for (i in 1:n){
    tree <- tree.bd(pars, max.t=t)
    if (is.null(tree)){ ## NULL trees have 0 taxa
      out[i] <- 0
    } else {
      out[i] <- Ntip(tree)
    }
  }
}

n_age <- lapply(ages, function(x) 
  tree.bd.n(pars=c(tree_lambda, tree_mu), t=x, n=1000))
saveRDS(n_age,file="output/n_age.rds")

## calculate where in the distribution of possible Ns does our tree land
perc_n <- vector(length=nrow(d))
for (j in 1:nrow(d)){
  perc_n[j] <- length(which(n_age[[j]] < d$n.clade[j]))
}
