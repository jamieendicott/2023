library(ggplot2)
library(GEOquery)
library(purrr)
library(reshape2)
library(transcripTools)
library(pheatmap)
library(viridisLite)
library(psych)

#looking at clusters of gene expression in fibroblast TERT experiment (bulk RNAseq)
g<-getGEO('GSE197471')
p<-(pData(g[[1]]))
#beautify pdata
p<-p[,c(1,2,27,55:60,62:64)]
cols<-(as.character(map(strsplit(colnames(p), split = ":"), 1)))
colnames(p)<-cols
p<-subset(p,p$coriell_id=="AG06561") #will contain TERT immortalization experiment and normal serial passaging

download.file('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE197nnn/GSE197471/suppl/GSE197471_norm_cts_GEO.tsv.gz','./GSE197471_norm_cts_GEO.tsv.gz')
cts<-read.table('GSE197471_norm_cts_GEO.tsv.gz',row.names=2)
dim(cts)
#[1] 19384    80
ncts<-cts[,c(match(p$description,colnames(cts)))]  
ncts<-cbind(cts[,1:4],ncts)


#trying to cluster meaningful groups of genes based on different expression patterns
#revisiting 020623
dat<-p[,c(9,10)]
rownames(dat)<-p$description
dat<-subset(dat,dat$expression_vector=='TERT')    # & dat$population_doublings>=60) #only after immortalization timepoint, not many timepoints so not super useful...
dat$population_doublings<-as.numeric(dat$population_doublings)

v.ncts<-ncts[,5:24]
rownames(v.ncts)<-ncts$Uniq_syms
v.ncts<-v.ncts[,c(match(rownames(dat),colnames(v.ncts)))]
v.ncts<-mostVar(v.ncts,1000) #larger sample clusters possibly driven by the single senescence timepoint

pheatmap(v.ncts[,c(match(rownames(dat),colnames(v.ncts)))],
         show_rownames = F,annotation_col = dat,cluster_cols = F,
         color=viridis(100))

#pull clusters based on corr distance
mat<-as.matrix(v.ncts)
cor<-cor(t(mat),method="pearson")
cor <- cor2dist(cor)

my_hclust <- hclust(dist(cor),method = 'ward.D2') 
my_clust <- as.data.frame(cutree(my_hclust, k = 4))
colnames(my_clust)<-"C"
my_clust$C<-as.factor(my_clust$C)

png("hm.png", res = 250, width = 1400, height = 1800)
pheatmap(v.ncts[,c(match(rownames(dat),colnames(v.ncts)))],
         annotation_row = my_clust,clustering_method = 'ward.D2',
         show_rownames = F,annotation_col = dat,cluster_cols = F,
         color=viridis(100),border_color = NA
        )
dev.off()

#here cluster 2 is of interest/oscillates (check for each run, n)

dat2<-cbind(dat,t(v.ncts))
mdat<-melt(dat2,id.vars = c('expression_vector','population_doublings'))
mdat$C<-my_clust[c(match(mdat$variable,rownames(my_clust))),1]

g<-ggplot(data=mdat,aes(x=population_doublings,y=value,col=C))
g+geom_point()+geom_smooth(method='loess')+
  theme_bw()+
  facet_wrap(~C)
#not impressive... can kind of see oscillation? heatmap more pronounced
#fuzzy c means might be more suitable than corr distance 

#running GO
#from basic template https://github.com/mousepixels/sanbomics_scripts/blob/main/GO_in_R.Rmd
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

#instead looking at specific clusters, here looks like cluster 1
genes_to_test<- subset(my_clust,my_clust$C==2) #check cluster...
GO_results <- enrichGO(gene = rownames(genes_to_test), OrgDb = "org.Hs.eg.db", keyType="SYMBOL", ont = "BP")
head(as.data.frame(GO_results))

fit <- plot(barplot(GO_results, showCategory = 20))
png("GOres.png", res = 250, width = 1400, height = 1800)
print(fit)
dev.off()
