library(tidyverse)
library(diversitree)

## data simulations
d <- read.csv("output/summary_tree_results.csv")

## calculate a mean of lambda and mu from old trees
tree_lambda <- d %>% filter(tree.max.age > 200) %>% 
  summarise(mean(mean.clade.lambda)) %>% as.numeric()

tree_mu <- d %>% filter(tree.max.age > 200) %>% 
  summarise(mean(mean.clade.mu)) %>% as.numeric()

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

## calculate where in the distribution of possible Ns does our tree land
perc_n <- vector(length=nrow(d))
for (j in 1:nrow(d)){
  perc_n[j] <- length(which(n_age[[j]] < d$n.clade[j]))
}



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

## calculate where in the distribution of possible branching times does our tree land
perc_t <- vector(length=nrow(d))
for (j in 1:nrow(d)){
  perc_t[j] <- length(which(n_t[[j]] < d$tree.max.age[j]))
}

