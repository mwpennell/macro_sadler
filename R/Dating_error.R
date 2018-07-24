library(tidyverse)
library(diversitree)
library(castor)

## data simulations
d <- read.csv("output/summary_tree_results.csv")

## calculate a mean of lambda and mu from old trees
tree_lambda <- d %>% filter(tree.max.age > 150) %>% 
  summarise(mean(mean.clade.lambda)) %>% as.numeric()

tree_mu <- d %>% filter(tree.max.age > 150) %>% 
  summarise(mean(mean.clade.mu)) %>% as.numeric()

#pars <- c(tree_lambda, tree_mu)
pars<-list(birth_rate_factor = tree_lambda,
           death_rate_factor = tree_mu)


## do n independent simulations and estimate slope for each of these
tree.bd.taxa <- function(pars, taxa, max.tries = 1000){
  while(max.tries > 0){
    tree <- generate_random_tree(parameters = pars, max_tips = taxa)$tree
    if (!is.null(tree)){
      if (Ntip(tree) == taxa) return(tree)
    }
    max.tries = max.tries - 1
  }
  return(NULL)
}

tree.bd.time <- function(pars, time, min.taxa = 10, max.tries = 1000){
  while(max.tries > 0){
    tree <- generate_random_tree(parameters = pars, max_time = time)$tree
    if (!is.null(tree)){
      if (Ntip(tree) == min.taxa) return(tree)
    }
    max.tries = max.tries - 1
  }
  return(NULL)
}

tree_bd_slope <- function(pars, taxa = NULL, time = NULL, error, min.taxa = 10, max.tries = 1000){
  if(!is.null(taxa)){
    trees <- lapply(taxa, function(x)
      tree.bd.taxa(pars, taxa = x, max.tries = max.tries))
  }else if(!is.null(time)){
    trees <- lapply(time, function(x)
      tree.bd.time(pars, time = x, min.taxa = min.taxa, max.tries = max.tries))
  } else {stop("Must provide stop condidtion.")}
  
  ind <- sapply(trees, class) == "phylo"
  trees <- trees[ind]
  if(is.null(time)){
    ages <- sapply(trees, function(x) max(branching.times(x)) )
  } else{
    ages <- time[ind]
  }
  
  lambdas <- vector(length = length(ages), mode = "numeric")
  for (i in 1:length(ages)){
    lik <- make.bd(trees[[i]])
    fit <- find.mle(lik, x.init=as.numeric(pars))
    if(fit$code %in% c(1,2)){
      lambdas[i] <- fit$par["lambda"]
    } else{
      lambdas[i] <- NA
    }
  }
  
  out = c()
  for(i in 1:length(error)){
    these_ages <- ages*runif(length(ages), min=1-error[i], max=1+error[i])
    fit_lm <- lm(log(lambdas)~log(these_ages), na.action = na.exclude)
    out[i] <- as.numeric(fit_lm$coefficients["log(these_ages)"])
  }
  names(out) <- as.character(error)
  

  return(out)
}

## filter so only look at trees with more than 10 tips
error <- seq(0.1, 0.9, length.out = 9)
Nreps <- 10
# Ntips
d_10 <- filter(d, ntips > 10)
ntips <- as.numeric(d_10$ntips)
slopes_taxa <- sapply(c(1:Nreps), function(x) tree_bd_slope(pars, taxa = ntips, error = error)) %>% t()



## estimate empirical slope
emp_10 <- lm(log(d_10$mean.clade.lambda)~log(d_10$tree.max.age))
emp_taxa <- as.numeric(emp_10$coefficients[2])


## calculate where in the distribution of possible slopes does our estimate land
perc_taxa <-  vector(length = length(error), mode = "numeric")
for (j in 1:length(error)){
  perc_taxa[j] <- sum(slopes_taxa[ , j] < emp_taxa)
}

## plot distribution of recovered slopes (due to dating bias) and include empircal slope
ggplot(perc, aes(x = value)) + geom_histogram(bins = 100) + theme_bw() + 
  xlab("slope estimates") + geom_vline(xintercept=0.05, color="red") +
  facet_grid(cond ~ .)
# plot the most extreme bias (i.e., slope = 10)
df <- data.frame(time = slopes_ages[ , "10"], ntip = slopes_taxa[ , "10"]) %>% gather(cond, value) %>% mutate(cond = as.factor(cond))
vline.df <- data.frame(cond = c("time", "ntip"), value = c(emp_ages, emp_taxa))
ggplot(df, aes(x = value)) + geom_histogram(bins = 100) + theme_bw() + 
  xlab("slope estimates") + 
  geom_vline(data = vline.df, aes(xintercept = value), color="red") +
  facet_grid(. ~ cond)






