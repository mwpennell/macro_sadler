ph1.1<-ggplot(summary_tree_results, aes(x=(tree.max.age),y=(mean.clade.lambda))) + 
  geom_point(show.legend = F) + theme_cowplot() + 
  labs(x="Clade age (Myr)", y=expression(paste("Mean ", lambda, " (species ", Myr^-1,")"),sep=" ")) + 
  theme(axis.title = element_text(size=15)) + 
  geom_smooth(method='glm',method.args=list(family="Gamma"), se=T,alpha=.15) 

ph3<-ggplot(summary_tree_results, aes(x=(tree.max.age),y=(mean.clade.mu))) + 
  geom_point() + theme_cowplot() + 
  labs(x="Clade age (Myr)", y=expression(paste("Mean ", mu, " (species ", Myr^-1,")"),sep=" ")) + theme(axis.title = element_text(size=15)) + 
  geom_smooth(method='glm',method.args=list(family="Gamma"), se=T,alpha=.15) 

l.inset<-ggplot(summary_tree_results, aes(x=log(tree.max.age),y=log(mean.clade.lambda))) + 
  geom_point() + theme_cowplot() + labs(x="", y="") +
  #labs(x="Ln clade age (Myr)", y=express) +
  theme(axis.title = element_text(size=10)) + 
  geom_smooth(method=lm, se=F,alpha=.1)

m.inset<-ggplot(summary_tree_results, aes(x=log(tree.max.age),y=log(mean.clade.mu))) + 
  geom_point() + theme_cowplot() + labs(x="", y="") +
  #labs(x="Ln clade age (Myr)", y=e.express) + 
  theme(axis.title = element_text(size=10)) + 
  geom_smooth(method=lm, se=F,alpha=.1)

gg1<-ggdraw()+ draw_plot(ph1.1, x = 0, y = 0) + draw_plot(l.inset, x = 0.6, y = 0.6, width = 0.4, height = 0.4)
gg2<-ggdraw()+ draw_plot(ph3, x = 0, y = 0) + draw_plot(m.inset, x = 0.6, y = 0.6, width = 0.4, height = 0.4)

ggarrange(gg1, gg2, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)


o.express<-expression(paste("Ln origination rate (genera ", Myr^-1,")"),sep=" ")
p1.1<-ggplot(summary_paleo_results,aes(x=(Duration),y=(mean.clade.origination))) + 
  geom_point( na.rm = T, show.legend = F) + theme_cowplot() + 
  labs(x="Duration (Myr)", y=expression(paste("Mean origination rate (genera ", Myr^-1,")"),sep=" ")) + 
  geom_smooth(method='glm',method.args=list(family="Gamma"), se=T,alpha=.15)
  
p2<-ggplot(summary_paleo_results,aes(x=(Duration),y=(mean.clade.mu))) + 
  geom_point(na.rm = T) + theme_cowplot() + 
  labs(x="Duration (Myr)", y=expression(paste("Mean ", mu, " (species ", Myr^-1,")"),sep=" ")) + 
  geom_smooth(method='glm',method.args=list(family="Gamma"), se=T,alpha=.15)

o.inset<-ggplot(summary_paleo_results, aes(x=log(Duration),y=log(mean.clade.origination))) + 
  geom_point() + theme_cowplot() + labs(x="", y="") +
  #labs(x="Ln clade duration (Myr)", y=o.express) + 
  theme(axis.title = element_text(size=10)) + 
  geom_smooth(method=lm, se=F,alpha=.1)

mf.inset<-ggplot(summary_paleo_results, aes(x=log(Duration),y=log(mean.clade.mu))) + 
  geom_point() + theme_cowplot() + labs(x="", y="") +
  #labs(x="Ln clade duration (Myr)", y=e.express) +
  theme(axis.title = element_text(size=10)) + 
  geom_smooth(method=lm, se=F,alpha=.1)

gg3<-ggdraw()+ draw_plot(p1.1, x = 0, y = 0) + draw_plot(o.inset, x = 0.6, y = 0.6, width = 0.4, height = 0.4)
gg4<-ggdraw()+ draw_plot(p2, x = 0, y = 0) + draw_plot(mf.inset, x = 0.6, y = 0.6, width = 0.4, height = 0.4)

ggarrange(gg3, gg4, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

# THE GRAPH
ggarrange(gg1, gg2, gg3, gg4, 
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)
