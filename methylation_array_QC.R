#after processing beta values through sesame pipe
s$na_count<-apply(betas,2, function(x) sum(is.na(x)))
#density plot
test<-t(betas)
test<-cbind(p[,c(8,10)],test)
p2<-melt(test,id.vars = c('PassageNumber','condition'))
g<-ggplot(data=p2,aes(x=value,col=condition))
pdf('pilotEE_density.pdf')
g+geom_density()+
  facet_wrap(~PassageNumber)+
  theme_dark()+
  scale_color_brewer(palette="Pastel1")
dev.off()
              
