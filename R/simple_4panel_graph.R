library(ggplot2)
library(cowplot)
library(MASS)
library(scales)
library(tidyverse)
library(ggpubr)

summary_tree_results<-read.table("output/summary_tree_results.csv", header = T, dec = ".", sep=",")
summary_tree_results<-filter(summary_tree_results,prop.samp >= 0.01)
summary_tree_results<-filter(str,ntips >6)
summary(summary_tree_results)

summary(lm(log(mean.clade.lambda)~log(tree.max.age),data=summary_tree_results))
summary(lm(log(mean.clade.mu)~log(tree.max.age),data=summary_tree_results))
summary(lm(log(net.diver+0.03)~log(tree.max.age),data=summary_tree_results))

summary(lm(log(var.clade.lambda)~log(tree.max.age),data=summary_tree_results))
summary(lm(log(var.clade.mu)~log(tree.max.age),data=summary_tree_results))

ph1.1<-ggplot(summary_tree_results, aes(x=(tree.max.age),y=(mean.clade.lambda))) + 
  geom_point(show.legend = F, colour= "#0E233E", size=2) + theme_cowplot() + 
  labs(x="Clade age (Myr)", y=expression(paste("Speciation rate (species ", Myr^-1,")"),sep=" ")) + 
  theme(axis.title = element_text(size=15)) + geom_smooth(method = 'nls', formula = y~a*x^b, method.args = list(start = c(a=1,b=2)), se=F, na.rm = T, colour="#EA3770", size=1.5)
  # + geom_smooth(method='glm',method.args=list(family="Gamma"), se=T,alpha=.15) 

ph3<-ggplot(summary_tree_results, aes(x=(tree.max.age),y=(mean.clade.mu))) + 
  geom_point() + theme_cowplot() + 
  labs(x="Clade age (Myr)", y=expression(paste("Extinction rate (species ", Myr^-1,")"),sep=" ")) + 
  theme(axis.title = element_text(size=15)) + geom_smooth(method = 'nls', formula = y~a*x^b, method.args = list(start = c(a=1,b=1)), se=F, na.rm = T, colour="#EA3770", size=1.5)
  # + geom_smooth(method='glm',method.args=list(family="Gamma"), se=T,alpha=.15) 

l.inset<-ggplot(summary_tree_results, aes(x=log(tree.max.age),y=log(mean.clade.lambda))) + 
  geom_point(colour= "#0E233E") + theme_cowplot() + labs(x="", y="") + 
  scale_x_continuous(breaks= c(1,3.5,6), labels=as.character(round(exp(c(1,3.5,6)),0))) +
  scale_y_continuous(breaks= c(-3.5,-1.5,0), labels=as.character(round(exp(c(-3.5,-1.5,0)),2))) +
  # labs(x="Ln clade age (Myr)", y=express) +
  theme(axis.title = element_text(size=10), axis.text.x = element_text(size=10), axis.text.y = element_text(size=10)) + 
  geom_smooth(method=lm, se=T,alpha=.1,colour="#EA3770")

m.inset<-ggplot(summary_tree_results, aes(x=log(tree.max.age),y=log(mean.clade.mu))) + 
  geom_point(colour= "#0E233E") + theme_cowplot() + labs(x="", y="") +
  scale_x_continuous(breaks= c(1,3.5,6), labels=as.character(round(exp(c(1,3.5,6)),0))) +
  scale_y_continuous(breaks= c(-6,-3,0), labels=as.character(round(exp(c(-6,-3,0)),2))) +
  #labs(x="Ln clade age (Myr)", y=e.express) + 
  theme(axis.title = element_text(size=10), axis.text.x = element_text(size=10), axis.text.y = element_text(size=10)) + 
  geom_smooth(method=lm, se=T,alpha=.1,colour="#EA3770")

gg1<-ggdraw()+ draw_plot(ph1.1, x = 0, y = 0) + draw_plot(l.inset, x = 0.5, y = 0.5, width = 0.45, height = 0.45)
gg2<-ggdraw()+ draw_plot(ph3, x = 0, y = 0) + draw_plot(m.inset, x = 0.5, y = 0.5, width = 0.45, height = 0.45)

ggarrange(gg1, gg2, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

o.express<-expression(paste("Origination rate (genera ", Myr^-1,")"),sep=" ")
p1.1<-ggplot(summary_paleo_results,aes(x=(Duration),y=(mean.clade.origination))) + 
  geom_point(colour= "#0E233E", size=2, na.rm = T, show.legend = F) + theme_cowplot() + 
  labs(x="Duration (Myr)", y=expression(paste("Origination rate (genera ", Myr^-1,")"),sep=" ")) + geom_smooth(method = 'nls', formula = y~a*x^b, method.args = list(start = c(a=1,b=1)), se=F, na.rm = T, colour="#EA3770", size=1.5)
#+ geom_smooth(method='glm',method.args=list(family="Gamma"), se=T,alpha=.15)
  
p2<-ggplot(summary_paleo_results,aes(x=(Duration),y=(mean.clade.mu))) + 
  geom_point(colour= "#0E233E", size=2, na.rm = T) + theme_cowplot() + 
  labs(x="Duration (Myr)", y=expression(paste("Extinction rate (genera ", Myr^-1,")"),sep=" ")) + geom_smooth(method = 'nls', formula = y~a*x^b, method.args = list(start = c(a=1,b=1)), se=F, na.rm = T, colour="#EA3770", size=1.5)
#+ geom_smooth(method='glm',method.args=list(family="Gamma"), se=T,alpha=.15)

o.inset<-ggplot(summary_paleo_results, aes(x=log(Duration),y=log(mean.clade.origination))) + 
  geom_point(colour= "#0E233E") + theme_cowplot() + labs(x="", y="") +
  scale_x_continuous(breaks= c(1,3.5,6), labels=as.character(round(exp(c(1,3.5,6)),0))) +
  scale_y_continuous(breaks= c(-3.5:-1), labels=as.character(round(exp(c(-3.5:-1)),2))) +
  theme(axis.title = element_text(size=10), axis.text.x = element_text(size=10), axis.text.y = element_text(size=10)) + 
  #labs(x="Ln clade duration (Myr)", y=o.express) + 
  geom_smooth(method=lm, se=T,alpha=.1, colour="#EA3770", size=1.5)

mf.inset<-ggplot(summary_paleo_results, aes(x=log(Duration),y=log(mean.clade.mu))) + 
  geom_point(colour= "#0E233E") + theme_cowplot() + labs(x="", y="") +
  scale_x_continuous(breaks= c(1,3.5,6), labels=as.character(round(exp(c(1,3.5,6)),0))) +
  scale_y_continuous(breaks= c(-3:-1), labels=as.character(round(exp(c(-3:-1)),2))) +
  theme(axis.title = element_text(size=10), axis.text.x = element_text(size=10), axis.text.y = element_text(size=10)) +
  #labs(x="Ln clade duration (Myr)", y=e.express) +
  theme(axis.title = element_text(size=10)) + 
  geom_smooth(method=lm, se=T,alpha=.1, colour="#EA3770", size=1.5)

gg3<-ggdraw()+ draw_plot(p1.1, x = 0, y = 0) + draw_plot(o.inset, x = 0.5, y = 0.5, width = 0.45, height = 0.45)
gg4<-ggdraw()+ draw_plot(p2, x = 0, y = 0) + draw_plot(mf.inset, x = 0.5, y = 0.5, width = 0.45, height = 0.45)

ggarrange(gg3, gg4, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

# THE GRAPH
ggarrange(gg1, gg2, gg3, gg4, 
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)
