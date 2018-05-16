library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(ggpubr)

summary_tree_results<-read.delim("output/summary_tree_results.csv", sep = ",", dec=".")
express<-expression(paste("Ln mean ", lambda, " (species ", Myr^-1,")"),sep=" ")
express1<-expression(paste("Ln mean ", mu, " (species ", Myr^-1,")"),sep=" ")

### Descriptive co-relations
a1<-ggplot(summary_tree_results, aes(x=log(tree.max.age),y=log(n.clade))) + geom_point() + theme_bw() + 
  labs(x="Log Clade age (Myr)", y= "Log Clade Richness") + theme(axis.title = element_text(size=15)) + 
  geom_smooth(method=lm, se=T,alpha=.2) 

a2<-ggplot(summary_tree_results, aes(x=log(tree.max.age),y=(gamma.stat))) + geom_point() + theme_bw() + 
  labs(x="Log Clade age (Myr)", y= "Gamma statistic") + theme(axis.title = element_text(size=15)) + 
  geom_smooth(method=lm, se=T,alpha=.2) 

ggarrange(a1, a2,  
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

g1<-ggplot(summary_tree_results, aes(x=(gamma.stat),y=log(mean.clade.lambda))) + 
  geom_point() + theme_bw() + 
  labs(x= "Gamma statistic", y=express) + theme(axis.title = element_text(size=15)) + 
  geom_smooth(method=lm, se=T,alpha=.2)

g2<-ggplot(summary_tree_results, aes(x=(gamma.stat),y=log(mean.clade.mu))) + 
  geom_point() + theme_bw() + 
  labs(x= "Gamma statistic", y=express1) + theme(axis.title = element_text(size=15)) + 
  geom_smooth(method=lm, se=T,alpha=.2)

ggarrange(g1, g2, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

g3<-ggplot(summary_tree_results, aes(x=(gamma.stat),y=log(net.diver))) + 
  geom_point() + theme_bw() + 
  labs(x= "Gamma statistic", y=expression(paste("Ln Net diversification (species ", Myr^-1,")"),sep=" ")) + theme(axis.title = element_text(size=15)) + 
  geom_smooth(method=lm, se=T,alpha=.2)

ggarrange(g1, g2, g3,  
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)

### Time-dependent diversification rates
p1<-ggplot(summary_tree_results, aes(x=log(tree.max.age),y=log(mean.clade.lambda))) + 
  geom_point(aes(size=n.clade,colour=gamma.stat)) + theme_bw() + 
  labs(x="Ln clade age (Myr)", y=express) + theme(axis.title = element_text(size=15)) + 
  geom_smooth(method=lm, se=F,alpha=.1) +
  scale_size_continuous(name = "Clade\nrichness") + 
 scale_colour_gradientn(name="Gamma\nstatistic",colours = rev(brewer.pal(4,"Spectral"))) + guides(size = guide_legend(order=1))

pp1<-ggplot(summary_tree_results, aes(x=(tree.max.age),y=(mean.clade.lambda))) + 
  geom_point(aes(size=n.clade,colour=gamma.stat)) + theme_bw() + 
  labs(x="Ln clade age (Myr)", y=express) + theme(axis.title = element_text(size=15)) + 
  geom_smooth(method='glm',method.args=list(family="Gamma"), se=F) + #formula=y~splines::bs(x,4)
  scale_size_continuous(name = "Clade\nrichness") + 
  scale_colour_gradientn(name="Gamma\nstatistic",colours = rev(brewer.pal(4,"Spectral"))) + 
  guides(size = guide_legend(order=1))

summary_paleo_results<-read.delim("output/summary_paleo_results.csv", sep = ",", dec=".")
long_paleo_results<-(gather(summary_paleo_results, key="Rate", value = "Estimate",mean.clade.mu:mean.clade.origination))
long_paleo_results$Rate<-recode(long_paleo_results$Rate, `mean.clade.mu`="Extinction", `mean.clade.origination`="Origination")

express2<-expression(paste("Ln estimate (genera ", Myr^-1,")"),sep=" ")
p2<-ggplot(long_paleo_results,aes(x=log(Duration),y=log(Estimate),shape=Rate, group=Rate)) + 
  geom_point(aes(size=Ngen,colour=freqRat)) + theme_bw() +labs(x="Ln duration (Myr)", y=express2) + 
  theme(axis.title = element_text(size=15)) + geom_smooth(aes(linetype=Rate),method=lm, se=F,alpha=.1) +
  scale_size_continuous(name = "Clade\ngenera") + 
  scale_colour_gradientn(name="Record\ncompletedness",colours = rev(brewer.pal(4,"Spectral"))) + 
  guides(size = guide_legend(order=1)) 

### Null trees (push of the past)

empirical_bd_slope<-readRDS("null_trees/empirical/output/tree_bd_slope.rds")
m50_bd_slope<-readRDS("null_trees/mu50/output/tree_bd_slope.rds")
m75_bd_slope<-readRDS("null_trees/mu75/output/tree_bd_slope.rds")
null_bd_slope<-data.frame(empirical_bd_slope,m50_bd_slope,m75_bd_slope)

h1<-ggplot(null_bd_slope, aes(x=empirical_bd_slope)) + geom_histogram(color="darkgreen", fill="white") + 
  theme_bw() + geom_vline(aes(xintercept=-0.436),color="orange", linetype="dashed", size=1) + 
  labs(title="Empirical parameters",x="Slope", y = "Count")

h2<-ggplot(null_bd_slope, aes(x=m50_bd_slope)) + geom_histogram(color="darkgreen", fill="white") + 
  theme_bw() + geom_vline(aes(xintercept=-0.436),color="orange", linetype="dashed", size=1) + 
  labs(title="mu = 0.5 x Lambda",x="Slope", y = "Count")

h3<-ggplot(null_bd_slope, aes(x=m75_bd_slope)) + geom_histogram(color="darkgreen", fill="white") + 
  theme_bw() + geom_vline(aes(xintercept=-0.436),color="orange", linetype="dashed", size=1) + 
  labs(title="mu = 0.75 x Lambda",x="Slope", y = "Count")

ggarrange(h1, h2, h3, 
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)
