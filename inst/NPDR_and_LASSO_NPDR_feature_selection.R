####### run npdr ####################
#install.packages("devtools")
#library(devtools)
#install_github("insilico/npdr")  

library(npdr)
SxG.dat <- t(GxS.covfilter)
min.group.size <- min( table ( pheno.factor ) )
k <- npdr :: knnSURF( 2 * min.group.size - 1 , 0.5 )
npdr.MDD.results <- npdr(phenos.df$age, SxG.dat, 
                         regression.type="lm",
                         attr.diff.type="numeric-abs",
                         #nbd.method="multisurf",
                         nbd.metric = "manhattan", 
                         #msurf.sd.frac=.5, 
                         k=k,
                         dopar.nn = F, dopar.reg=F,
                         padj.method="bonferroni", verbose = T)

colnames(npdr.MDD.results)
dim(npdr.MDD.results)
library(dplyr)
# attributes with FDR-adjusted p-value<.05
top.p05.npdr <- npdr.MDD.results %>% filter(pval.adj<.05) %>% pull(att) 
top.p05.npdr <- npdr.MDD.results %>% filter(pval.adj<.05)  %>% 
  #mutate_at("pval.adj", as.character) %>%
  mutate_at("pval.adj", as.numeric) %>%
  arrange(pval.adj) 

write.table(top.p05.npdr%>% pull(att) ,row.names=F,col.names=F,quote=F)

# grab top 200, remove NA, remove "", get att col
top.npdr <- npdr.MDD.results %>% dplyr::slice(1:400) %>% 
  filter(att!="") %>% pull(att)
write.table(top.npdr,row.names=F,col.names=F,quote=F)

npdr.MDD.results %>% filter(pval.adj<.05) %>% dplyr::select(att,pval.adj)
library(dplyr)
NPDR_plot <- top.p05.npdr %>%
  filter( pval.adj<.05 ) %>%
  ggplot( aes(x=pval.adj)) +
  geom_histogram( binwidth=0.01, fill="#69b3a2", color="#e9ecef", alpha=0.9,boundary=0) +
  ggtitle("Age associated genes selected by NPDR") +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=15)
  )+  
  stat_bin(binwidth=0.01, geom="text", aes(label=..count..), vjust=1,boundary=0)+
  scale_x_continuous(name = "Adjusted p threshold",
                     breaks = seq(0, 0.05, 0.01)) +
  scale_y_continuous(name = "Count")+
  theme_bw() +
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 13, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma"),
        axis.text.x=element_text(colour="black", size = 13),
        axis.text.y=element_text(colour="black", size = 13))
NPDR_plot

# get the plot for the gene has lowest adj p value 
plot_row_npdr <- which(rownames(GxS.covfilter)=="SESTD1")
plot(phenos.df$age,GxS.covfilter[plot_row_npdr,],xlab="age",ylab="SESTD1")


######### lasso-NDPR ##############################
library(npdr)
SxG.dat <- t(GxS.covfilter)
min.group.size <- min( table ( pheno.factor ) )
k <- npdr :: knnSURF( 2 * min.group.size - 1 , 0.5 )

npdr_lasso_results <- npdr::npdr(phenos.df$age, SxG.dat, 
                                 regression.type="lm",
                                 attr.diff.type="numeric-abs",
                                 nbd.method="multisurf",
                                 nbd.metric = "manhattan",
                                 msurf.sd.frac=.5,
                                 knn=k,
                                 use.glmnet = T,
                                 glmnet.alpha = 1,
                                 neighbor.sampling="none", dopar.nn = F,
                                 padj.method="bonferroni", verbose=T)

# attributes with score/coefficient > 0
top.score.npdr <- rownames(npdr_lasso_results)[npdr_lasso_results$scores>0]
write.table(top.score.npdr ,row.names=F,col.names=F,quote=F)
top.score.npdr <- npdr_lasso_results %>% filter(scores>0)  %>% 
  #mutate_at("pval.adj", as.character) %>%
  mutate_at("scores", as.numeric) %>% data.frame()


#intersect(intersect(intersect(glmnet.age.nonzero.coeffs, lm.results.sorted.0.05$gene), top.p05.npdr),top.score.npdr)
write.table(top.score.npdr%>% pull(scores) ,row.names=F,col.names=F,quote=F)
LassoNPDR_plot <- npdr_lasso_results %>%
  filter( scores>0 ) %>%
  ggplot( aes(x=scores)) +
  geom_histogram( binwidth=1, fill="#69b3a2", color="#e9ecef", alpha=0.9,boundary=0) +
  ggtitle("Age associated genes selected by LASSO NPDR") +
  #theme_ipsum() +
  theme(
    plot.title = element_text(size=15)
  )+  
  stat_bin(binwidth=1, geom="text", aes(label=..count..), vjust=1,boundary=0)+
  scale_x_continuous(name = "Coefficient") +
  scale_y_continuous(name = "Count")+
  theme_bw() +
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 13, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma"),
        axis.text.x=element_text(colour="black", size = 13),
        axis.text.y=element_text(colour="black", size = 13))
LassoNPDR_plot

