####### linear regression of age ~ gene expression ##########

# could do upper and lower quartiles of age and then
# could do logistic regression of gene on hi/lo age

lm.age.fn <- function(i){
  gene=rownames(GxS.covfilter)[i]
  gene.expr <- GxS.covfilter[i,] 
  gene.fit <- lm(phenos.df$age~gene.expr) # gene.expr uses the GxS.covfilter for gene i
  coeff.mat <- coef(summary(gene.fit)) #function of age~gene.expr
  b1 <- coeff.mat[2,1]
  b1.pval <- coeff.mat[2,4]
  coefvec <- gene.fit$estimate # intercept, gene
  pvec <- gene.fit$p.value     # intercept, gene
  c(gene, b1, b1.pval)    
} 

lm.age.fn(1)

# for loop the function to all genes
num.genes<-nrow(GxS.covfilter)
lm.results.mat <- matrix(0,nrow=nrow(GxS.covfilter),
                         ncol=3)
for (i in 1:num.genes){
  lm.results.mat[i,] <- lm.age.fn(i) 
}
lm.results.df <- data.frame(lm.results.mat)
colnames(lm.results.df) <- c("gene", "b1", "p.val")

# sort results by slope coefficient p-value
# sort by p-value
library(dplyr)
lm.results.sorted <- lm.results.df %>% 
  mutate_at("p.val", as.character) %>%
  mutate_at("p.val", as.numeric) %>%
  arrange(p.val) 
lm.results.sorted$p.adj <- p.adjust(lm.results.sorted$p.val,method="fdr")

# adjusted p-value feature selection
#lm.results.sorted.0.05 <- subset(lm.results.sorted, lm.results.sorted$p.val<= 0.05)
library(ggplot2)
library(tidyverse)
library(hrbrthemes)
LR_plot <- lm.results.sorted %>%
  filter( p.val<= 0.05 ) %>%
  ggplot( aes(x=p.adj)) +
  geom_histogram( binwidth=0.05, fill="#69b3a2", color="#e9ecef", alpha=0.9,boundary=0) +
  ggtitle("464 age associated genes selected by linear regression with p<0.05") +
  theme_ipsum() +
  theme(plot.title = element_text(size=15))+  
  stat_bin(boundary=0,binwidth=0.05, geom="text", aes(label=..count..), vjust=2)+
  geom_vline(xintercept = 0.05, size = 1, colour = "#FF3721",
             linetype = "dashed")+
  scale_x_continuous(name = "Adjusted p threshold",
                     breaks = seq(0, 0.7, 0.05)) +
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


# make a list of genes with p<0.05
lm.results.sorted.0.05 <- lm.results.sorted %>% filter(p.val<=0.05)%>% dplyr::select(`gene`)%>% c()
write.table(lm.results.sorted.0.05$gene[1:200],row.names=F,col.names=F,quote=F)
#lm.results.sorted.0.05$gene[1:200]
# import the reactome result
setwd("/Users/liyijie/Documents/Gene age/")
lm.results.sorted.age.overlap <- read.table(file = 'lm.gene.age.overlap.tsv', skip=24, sep = '\t', header = TRUE)
colnames(lm.results.sorted.age.overlap)

#Genes in the infectious disease pathway
lm.results.sorted.age.overlap %>% filter(lm.results.sorted.age.overlap$REACTOME_INFECTIOUS_DISEASE == "REACTOME_INFECTIOUS_DISEASE")%>% dplyr::select(`Gene.Symbol`)%>% c()

#Genes in the SARS-CoV Infections pathway
lm.results.sorted.age.overlap %>% filter(lm.results.sorted.age.overlap$REACTOME_SARS_COV_INFECTIONS == "REACTOME_SARS_COV_INFECTIONS")%>% dplyr::select(`Gene.Symbol`)%>% c()

#intersection
lm.results.sorted.age.overlap %>% filter(lm.results.sorted.age.overlap$REACTOME_INFECTIOUS_DISEASE == "REACTOME_INFECTIOUS_DISEASE",
                                         lm.results.sorted.age.overlap$REACTOME_SARS_COV_INFECTIONS == "REACTOME_SARS_COV_INFECTIONS")%>% dplyr::select(`Gene.Symbol`)%>% c()

