#Gene age paper code
#2023
#Jamie Li

#=========================================================
# Data processing #######
#=========================================================

# setwd to data/ directory

# load gene expression data
load("sense.filtered.cpm.Rdata")

dim(sense.filtered.cpm) # 157 subjects 
# load phenotype (mdd/hc) data
subject.attrs <- read.csv("Demographic_symptom.csv", 
                          stringsAsFactors = FALSE)
#install.packages("dplyr")
library(dplyr)
# grab intersecting X (subject ids) and Diag (Diagnosis) from columns
phenos.df <- subject.attrs %>% 
  filter(X %in% colnames(sense.filtered.cpm)) %>%
  dplyr::select(X, Diag, age, sex)  
mddPheno <- as.factor(phenos.df$Diag)  

# Normalized and transform
#install.packages("preprocessCore")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("preprocessCore")
library(preprocessCore)
mddExprData_quantile <- normalize.quantiles(sense.filtered.cpm)
mddExprData_quantileLog2 <- log2(mddExprData_quantile)
colnames(mddExprData_quantileLog2) <- mddPheno  # add phenotype names to matrix
rownames(mddExprData_quantileLog2) <- rownames(sense.filtered.cpm) 

# B. COV filter
# coefficient of variation filter sd(x)/abs(mean(x))
# there are a lot of genes that have very low signal to noise that we can get rid of.
# 1 apply to rows
CoV_values <- apply(mddExprData_quantileLog2,1,
                    function(x) {sd(x)/abs(mean(x))})
# there is one gene that has 0 variation -- remove
sd_values <- apply(mddExprData_quantileLog2,1,
                   function(x) {sd(x)})
rownames(mddExprData_quantileLog2)[sd_values==0]

# the smaller the threshold, the higher the experimental effect relative to the measurement precision
sum(CoV_values<.045)  # 5,588 genes
# filter the data matrix and transpose
GxS.covfilter <- mddExprData_quantileLog2[CoV_values<.045 &
                                            sd_values>0,]
dim(GxS.covfilter)
colnames(GxS.covfilter)

############# Should we remove the outlier?
mddCorr<-cor(GxS.covfilter)  # distance based on correlation
d <- 1-mddCorr
dim(d)
rownames(d)
mddTree = hclust(as.dist(d))
#mddTree$labels <- mddPheno
mddTree$labels <- phenos.df$X
plot(mddTree)
table(phenos.df$sex)
table(phenos.df$Diag)
table(phenos.df$age)

# convert phenotype to factor
pheno.factor <- as.factor(colnames(GxS.covfilter))
pheno.factor
str(pheno.factor)
levels(pheno.factor)

####
t.test(phenos.df$age~as.factor(mddPheno)) #p-value = 0.1627
t.test(phenos.df$age~phenos.df$Diag) #p-value = 0.1627
