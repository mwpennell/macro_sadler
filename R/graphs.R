library(ggplot2)
library(RColorBrewer)
library(tidyverse)

summary_tree_results<-read.delim("output/summary_tree_results.csv", sep = ",", dec=".")

express<-expression(paste("Ln mean ", lambda, " (species ", Myr^-1,")"),sep=" ")
p1<-ggplot(summary_tree_results, aes(x=log(tree.max.age),y=log(mean.clade.lambda))) + 
  geom_point(aes(size=n.clade,colour=gamma.stat)) + theme_bw() + 
  labs(x="Ln clade age (Myr)", y=express) + theme(axis.title = element_text(size=15)) + 
  geom_smooth(method=lm, se=F,alpha=.1) +
  scale_size_continuous(name = "Clade\nrichness") + 
 scale_colour_gradientn(name="Gamma\nstatistic",colours = rev(brewer.pal(4,"Spectral"))) + guides(size = guide_legend(order=1))

pp1<-ggplot(summary_tree_results, aes(x=(tree.max.age),y=(mean.clade.lambda))) + 
  geom_point(aes(size=n.clade,colour=gamma.stat)) + theme_bw() + 
  labs(x="Ln clade age (Myr)", y=express) + theme(axis.title = element_text(size=15)) + 
  geom_smooth(method='glm',formula=y~splines::bs(x,4), se=F) + #method.args=list(family="Gamma")
  scale_size_continuous(name = "Clade\nrichness") + 
  scale_colour_gradientn(name="Gamma\nstatistic",colours = rev(brewer.pal(4,"Spectral"))) + guides(size = guide_legend(order=1))

summary_paleo_results<-read.delim("output/summary_paleo_results.csv", sep = ",", dec=".")
long_paleo_results<-(gather(summary_paleo_results, key="Rate", value = "Estimate",mean.clade.mu:mean.clade.origination))
long_paleo_results$Rate<-recode(long_paleo_results$Rate, `mean.clade.mu`="Extinction", `mean.clade.origination`="Origination")

express2<-expression(paste("Ln estimate (genera ", Myr^-1,")"),sep=" ")
p2<-ggplot(long_paleo_results,aes(x=log(Duration),y=log(Estimate),shape=Rate, group=Rate)) + 
  geom_point(aes(size=Ngen,colour=freqRat)) + theme_bw()

p2+labs(x="Ln duration (Myr)", y=express2) + theme(axis.title = element_text(size=15)) + 
  geom_smooth(aes(linetype=Rate),method=lm, se=F,alpha=.1)+
  scale_size_continuous(name = "Clade\ngenera") + 
  scale_colour_gradientn(name="Record\ncompletedness",colours = rev(brewer.pal(4,"Spectral"))) + 
  guides(size = guide_legend(order=1)) 
