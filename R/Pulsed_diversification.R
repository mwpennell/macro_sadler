library(tidyverse)
library(diversitree)
summary_tree_results <- read.csv("output/summary_tree_results.csv")
tree_lambda <- summary_tree_results %>% filter(tree.max.age > 50) %>% summarise(mean(mean.clade.lambda)) %>% as.numeric()
tree_mu <- summary_tree_results %>% filter(tree.max.age > 50) %>% summarise(mean(mean.clade.mu)) %>% as.numeric()




times <- rpois(n=1, lambda = 26/5)
total_depth <- round(max(summary_tree_results$tree.max.age)) / 5
all_times <- c(1:total_depth)
new_birth_rates <- rep(tree_lambda, length(all_times))
new_death_rates <- rep(tree_mu, length(all_times))
## interval distance of two
while(1){
  temp <- rpois(n=1, lambda=26/5) + times[length(times)]
  if(temp+1 > total_depth){
    break()
  }
  aux <- rexp(1, rate=1/tree_lambda)
  new_birth_rates[c(temp, temp+1)] <- aux + tree_lambda
  new_death_rates[c(temp, temp+1)] <- new_birth_rates[c(temp, temp+1)] * 0.7
  times[length(times)+1] <- temp
}

birth_fun = stepfun(x = all_times, y = c(tree_lambda, new_birth_rates))
death_fun = stepfun(x = all_times, y = c(tree_mu, new_death_rates))
curve(birth_fun(x), 0, total_depth, col = "blue", ylim = c(0, 1))
curve(death_fun(x), 0, total_depth, col = "red", add = T)

ages <- summary_tree_results$tree.max.age / 5

trees <- list()
for (i in 1:length(ages)){
  trees[[i]] <- try(rbdtree(birth = birth_fun, death = death_fun, Tmax = ages[i]))
}
phylos = trees[sapply(trees, class) == "phylo"]
sapply(phylos, Ntip)

aux = lapply(phylos, drop.fossil)
sapply(aux, Ntip)
