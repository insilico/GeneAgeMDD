# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0
################################################################
#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(hrbrthemes)

#library(umap)

# load series and platform data from GEO
options(timeout = 300)
gset <- getGEO("GSE98793", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- paste0("11111111111111111111111111111111000000000000000000",
               "00000000000000XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX10X0",
               "X1X100X1X01XX000000011100XXXX0X10X1XX11X1111X1XX01",
               "001X111X000001XX0XXXX010001XX1X0111X001X11")
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("MDD","HC"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

gset <- gset[complete.cases(exprs(gset)), ] # skip missing values
#gset.expression <- data.frame(t(gset@assayData$exprs))
gset.expression <- data.frame(gset@assayData$exprs)
is.numeric(gset.expression)
gset.expression <- as.matrix(gset.expression)
dim(gset.expression)
is.numeric(t(gset.expression))


# Convert all columns to numeric
#gset.expression <- gset.expression %>%
#mutate_all(as.numeric)

gset.gene <- gset@featureData@data
gset.pheno <- gset@phenoData@data
gset.pheno$age <- as.numeric(gset.pheno$`age:ch1`)
dim(gset.expression)


####### map our top LGAGE genes to GSE Gene Symbols. #####
#Gene age paper code
#2023
#Jamie Li

#=========================================================
# Data processing #######
#=========================================================

setwd("/Users/liyijie/Documents/bioinform class")
# gene_age.R do with or after lab5

# load gene expression data
load("sense.filtered.cpm.Rdata")

dim(sense.filtered.cpm) # 8923 157 subjects 
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
# 行名：基因 列名：症状

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

####### Lasso of age ~ gene expression ####################

#install.packages("glmnet")
library(glmnet)
#remotes::install_github("const-ae/ggupset")
library(ggupset)
library(ggplot2)
library(dplyr)

# alpha=0 means ridge, alpha=1 means lasso
set.seed(2000)
glmnet.age<-cv.glmnet(t(GxS.covfilter),phenos.df$age,alpha=1,family="gaussian")
glmnet.age$lambda.1se #use 1.636048 to rerun the model
glmnet.age$lambda.min


glmnet.age.coeffs<-predict(glmnet.age,type="coefficients")
plot(glmnet.age) # MSE be lowest when log lambda around 0
glmnet.age.nonzero.coeffs <- glmnet.age.coeffs@Dimnames[[1]][which(glmnet.age.coeffs!=0)]

write.table(glmnet.age.nonzero.coeffs,row.names=F,col.names=F,quote=F)
length(glmnet.age.nonzero.coeffs)

glmnet.age.nonzero.coeffs.df <- data.frame(gene=glmnet.age.nonzero.coeffs[-1],
                                           coefficient=glmnet.age.coeffs[which(glmnet.age.coeffs!=0)][-1])

#intersect(glmnet.age.nonzero.coeffs.df$gene, glmnet.gene$Gene.symbol[-1])


#construct the GSE L-GAGE only using the genes we found
# Modify matching_genes to select only exact matches when multiple results are found
matching_genes_old_top <- sapply(glmnet.age.nonzero.coeffs.df$gene, function(gene) {
  # Get all matches for the current gene
  matches <- gset.gene$Gene.symbol[grep(gene, gset.gene$Gene.symbol)]
  
  # If there's more than one match, select the exact match
  if (length(matches) > 1) {
    exact_match <- matches[matches == gene]
    if (length(exact_match) > 0) return(exact_match)
  }
  
  # Return all matches if there's only one or no exact match
  return(matches)
})

# Filter to keep only the non-NA results
matching_genes_old_top <- matching_genes_old_top[!is.na(matching_genes_old_top)]
matching_genes_df <- unlist(matching_genes_old_top)
matching_genes_df <- unique(unlist(matching_genes_old_top))
# Display matching results
#matching_genes

# Display matching results
# Load necessary package
library(tidyr)
#matching_genes$AACS
# Convert matching_genes list to a data frame
matching_genes_df <- data.frame(gene = names(matching_genes_old_top), matched_symbols = I(matching_genes_old_top))

# Unnest the data frame to create a long format where each match is a separate row
matching_genes_df <- unnest(matching_genes_df, matched_symbols)

# View the resulting data frame
head(matching_genes_df)

intersect(glmnet.age.nonzero.coeffs.df$gene, matching_genes_df$gene)
gset.expression.old.gene <- subset(gset.gene, gset.gene$Gene.symbol %in%(matching_genes_df$matched_symbols))
dim(gset.expression.old.gene)

# rerun the model by only using the matching_genes_df
gset.expression.matching <- gset.expression[gset.expression.old.gene$ID,]

gset.expression.matching <- data.frame(Gene.symbol=gset.expression.old.gene$Gene.symbol,
                                       gset.expression[gset.expression.old.gene$ID,])
gset.expression.matching$Gene.symbol

library(dplyr)

gset.expression.matching.avg <- gset.expression.matching %>%
  group_by(Gene.symbol) %>%
  summarize(across(starts_with("GSM"), mean, na.rm = TRUE))

gset.expression.matching.avg <- gset.expression.matching.avg %>%
  mutate(Gene.symbol = ifelse(Gene.symbol == "AK6///TAF9", "TAF9", Gene.symbol))%>%
  mutate(Gene.symbol = ifelse(Gene.symbol == "TGIF2-C20orf24///C20orf24", "TGIF2-C20orf24", Gene.symbol))%>%
  filter(Gene.symbol != "TAF9B")

intersect(glmnet.age.nonzero.coeffs.df$gene,gset.expression.matching.avg$Gene.symbol)

# alpha=0 means ridge, alpha=1 means lasso
set.seed(2000)
dim(t(gset.expression.matching.avg[,-1]))
colnames(t(gset.expression.matching.avg[,-1]))
#glmnet.age.samegene<-cv.glmnet(t(gset.expression.matching.avg[,-1]),gset.pheno$age,alpha=0,family="gaussian")
glmnet.age.samegene <- cv.glmnet(t(gset.expression.matching.avg[,-1]), gset.pheno$age, alpha = 0, family = "gaussian", lambda = seq(0.001, 2, by = 0.01))

glmnet.age.samegene$lambda.1se #use 1.991 to rerun the model
glmnet.age.samegene$lambda.min #1.991

#gset.expression.matching.avg[,-1]
gset.expression.matching.avg$Gene.symbol
glmnet.age.samegene.coeffs<-predict(glmnet.age.samegene,type="coefficients")
# Example of reducing lambda
glmnet.age.samegene.coeffs.df <- data.frame(c(gset.expression.matching.avg$Gene.symbol),
                                            glmnet.age.samegene.coeffs@x[-1])


# I run the 19 genes only but didn't find any coefficient

glmnet_predict <- predict(glmnet.age.samegene, newx=t(gset.expression.matching.avg[,-1]), type="response")

library(tidyverse)
library(caret)

R2(glmnet_predict,gset.pheno$age)
RMSE(glmnet_predict,gset.pheno$age) # 7.679305
#RMSE(glmnet_predict,phenos.df$age)/mean(phenos.df$age) #0.2412267
plot(gset.pheno$age,glmnet_predict)
ggplot(data=data.frame(age=gset.pheno$age,glmnet_age=c(glmnet_predict)),
       aes(x=age,y=glmnet_age)) + geom_point() + geom_smooth(method='lm', formula= y~x) + ggtitle("glmnet model of age")

#############
gset.data1 <- data.frame(age=gset.pheno$age, gene=c(glmnet_predict), sex = ifelse(gset.pheno$`gender:ch1`=="F", "Female", "Male"),
                         predicted = predict(lm(glmnet_predict[,1]~gset.pheno$age)),
                         res = residuals(lm(glmnet_predict[,1]~gset.pheno$age)),
                         res.sign = ifelse(residuals(lm(glmnet_predict[,1]~gset.pheno$age))>=0, 1, 0),
                         phenotype = if_else(gset.pheno$`subject group:ch1`=="CNTL; healthy control", "HC", "MDD"))
gset.data1$abres = abs(gset.data1$res)


#heteroskedasticity by plotting the residuals on the y axis and the predicted values on the x axis. 
plot(gset.data1$gene,gset.data1$res)
#install.packages("lmtest")
library(lmtest)
lmtest::bptest(fit.g)
car::ncvTest(fit.g)

ggplot(data = gset.data1, aes(y = res, x = gene)) + geom_point(col = 'blue') + geom_abline(slope = 0) + ggtitle("Test heteroskedasticity")
#par(mfrow=c(2,2))
plot(fit.g)
dev.off()
#and then test whether the residual is associated with MDD, like this
#t.test(gset.data1$phenotype~gset.data1$res)
#and add also see what happens when you adjust for age again
#glm(as.numeric(gset.data1$phenotype)~gset.data1$res+gset.data1$age) 

table(gset.data1$res.sign,gset.data1$phenotype)
prop.table(table(gset.data1$res.sign,gset.data1$phenotype), margin = 2)
table(gset.data1$res.sign,gset.data1$sex)
# histogram of the ages with a bin size of 1 or 2
h <- ggplot(gset.data1, aes(x=age)) + 
  geom_histogram(color="black", fill="pink",binwidth=1)+ 
  ggtitle("Histogram of the ages with a bin size of 1")  +
  theme(plot.title = element_text(size=20,face="bold"),
        axis.text=element_text(size=18),
        axis.title=element_text(size=18,face="bold"),
        legend.text = element_text(size=18),
        legend.title = element_text(size=18,face="bold"))+  
  stat_bin(binwidth=1, geom="text", aes(label=..count..), vjust=-1) 
h
table(gset.data1$age)

# Separate plot 

ggplot(gset.data1,
       aes(x=age,y=gene)) + 
  xlab("chronological age")+
  ylab("Prediction age")+
  geom_smooth(method='lm', formula= y~x) + 
  geom_abline(intercept = 0, slope = 1,color="darkred") + 
  ggtitle("glmnet model of age---by sex") +
  geom_segment(aes(color = sex,xend = as.numeric(age), yend = as.numeric(predicted)), alpha = .5) +
  geom_point(aes(col=sex,size = 0.2)) +
  theme(plot.title = element_text(size=20,face="bold"),
        axis.text=element_text(size=18),
        axis.title=element_text(size=18,face="bold"),
        legend.text = element_text(size=18),
        legend.title = element_text(size=18,face="bold")) +
  guides(colour = guide_legend(override.aes = list(size=5))) + labs(fill = "Sex and phenotype")


ggplot(gset.data1,
       aes(x=age,y=gene)) + 
  xlab("Chronological age (years)")+
  ylab("Biological age (years)")+
  geom_smooth(method='lm', formula= y~x) + 
  #geom_abline(intercept = 0, slope = 1,color="darkred") + 
  #ggtitle("glmnet model of age---by phenotype") +
  #geom_segment(aes(color = phenotype,xend = as.numeric(age), yend = as.numeric(predicted)), alpha = .5) +
  geom_point(aes(col=phenotype,shape=sex),size = 4) + #size = abs(res)
  guides(colour = guide_legend(override.aes = list(size=5))) + labs(fill = "Sex and phenotype")+
  #geom_vline(aes(xintercept=20),colour="#BB0000",linetype="dashed")+
  #geom_vline(aes(xintercept=30),colour="#BB0000",linetype="dashed")+
  #geom_vline(aes(xintercept=40),colour="#BB0000",linetype="dashed")+
  theme_bw()+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 15, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma"),
        axis.text.x=element_text(colour="black", size = 11),
        axis.text.y=element_text(colour="black", size = 11),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text=element_text(size=18),
        axis.title=element_text(size=18,face="bold"),
        legend.text = element_text(size=18),
        legend.title = element_text(size=18,face="bold"))


table(gset.data1$phenotype, gset.data1$res.sign)
table(gset.data1$sex, gset.data1$res.sign)
# Ensure phenotype is a factor with exactly two levels
#data1$phenotype <- as.factor(data1$phenotype)
# Perform the t-test
t.test(res ~ phenotype, data = gset.data1)
# Test the significance of the correlation between x and y
spearman_test <- cor.test(gset.data1$age,gset.data1$gene, method = "spearman")
#print(spearman_test)
summary(lm(res ~ phenotype, data = gset.data1))

# Density plots 
g.phenotype <- ggplot(data=gset.data1,
                      aes(x=res, colour=phenotype, fill=phenotype)) +
  geom_density(alpha = .3) + 
  #ggtitle("A") +
  xlab(label = "L-GAGE")+
  theme_bw()+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 15, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma"),
        axis.text.x=element_text(colour="black", size = 11),
        axis.text.y=element_text(colour="black", size = 11),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"),
        legend.text = element_text(size=13),
        legend.title = element_text(size=13,face="bold"))

g.sex <- ggplot(data=gset.data1,
                aes(x=res, colour=sex, fill=sex)) +
  geom_density(alpha = .3) + 
  ggtitle("B") +
  xlab(label = "L-GAGE")+
  theme_bw()+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 15, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma"),
        axis.text.x=element_text(colour="black", size = 11),
        axis.text.y=element_text(colour="black", size = 11),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"),
        legend.text = element_text(size=13),
        legend.title = element_text(size=13,face="bold"))
library(ggpubr)
ggarrange(g.phenotype, g.sex,
          #labels = c("(a)", "(b)", "(c)","(d)","(e)"),
          ncol = 2, nrow = 1,widths = 1,vjust = 1.5)

# Test the significance of the correlation between x and y

t.test(gset.data1$age~as.factor(gset.data1$phenotype))
t.test(gset.data1$res~as.factor(gset.data1$phenotype))
summary(gset.data1$res)

nameUI <- function(id) {
  ns <- NS(id)
  tagList(
    
  )
}

nameServer <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      
    }
  )
}
###### for all genes ########
corr.GxS.covfilter <- cor(t(GxS.covfilter))
corr.GxS.covfilter.old.gene <- corr.GxS.covfilter[,glmnet.age.nonzero.coeffs.df$gene]

# Create an empty list to store the top 20 row names for each column
top20_row_names <- list()

# Loop over each column in the matrix
for (col_name in colnames(corr.GxS.covfilter.old.gene)) {
  # Sort the column in descending order and get the top 20 row names
  top20 <- head(order(corr.GxS.covfilter.old.gene[, col_name], decreasing = TRUE), 20)
  
  # Store the row names in the list
  top20_row_names[[col_name]] <- rownames(corr.GxS.covfilter.old.gene)[top20]
}

# Display the top 20 row names for each column
top20_row_names
# Combine all the lists into a single vector
all_top20_row_names <- unlist(top20_row_names)

# Find the unique values
unique_top20_row_names <- unique(all_top20_row_names)

# Display the unique row names
unique_top20_row_names



#construct the GSE L-GAGE only using the genes we found
# Modify matching_genes to select only exact matches when multiple results are found
matching_genes_old_topcorr <- sapply(unique_top20_row_names, function(gene) {
  # Get all matches for the current gene
  matches <- gset.gene$Gene.symbol[grep(gene, gset.gene$Gene.symbol)]
  
  # If there's more than one match, select the exact match
  if (length(matches) > 1) {
    exact_match <- matches[matches == gene]
    if (length(exact_match) > 0) return(exact_match)
  }
  
  # Return all matches if there's only one or no exact match
  return(matches)
})

# Filter to keep only the non-NA results
matching_genes_old_topcorr <- matching_genes_old_topcorr[!is.na(matching_genes_old_topcorr)]
matching_genes_df_old_topcorr <- unlist(matching_genes_old_topcorr)
matching_genes_df_old_topcorr <- unique(unlist(matching_genes_old_topcorr))
# Display matching results
#matching_genes

# Display matching results
# Load necessary package
library(tidyr)
#matching_genes$AACS
# Convert matching_genes list to a data frame
matching_genes_df_old_topcorr <- data.frame(gene = names(matching_genes_old_topcorr), matched_symbols = I(matching_genes_old_topcorr))

# Unnest the data frame to create a long format where each match is a separate row
matching_genes_df_old_topcorr <- unnest(matching_genes_df_old_topcorr, matched_symbols)

# View the resulting data frame
head(matching_genes_df_old_topcorr)

intersect(unique_top20_row_names, matching_genes_df_old_topcorr$gene)
gset.expression.old.gene_topcorr <- subset(gset.gene, gset.gene$Gene.symbol %in%(matching_genes_df_old_topcorr$matched_symbols))
dim(gset.expression.old.gene_topcorr)

# rerun the model by only using the matching_genes_df_old_topcorr
gset.expression.matching_topcorr <- gset.expression[gset.expression.old.gene_topcorr$ID,]

gset.expression.matching_topcorr <- data.frame(Gene.symbol=gset.expression.old.gene_topcorr$Gene.symbol,
                                               gset.expression[gset.expression.old.gene_topcorr$ID,])
gset.expression.matching_topcorr$Gene.symbol

library(dplyr)

gset.expression.matching_topcorr.avg <- gset.expression.matching_topcorr %>%
  group_by(Gene.symbol) %>%
  summarize(across(starts_with("GSM"), mean, na.rm = TRUE))

gset.expression.matching_topcorr.avg <- gset.expression.matching_topcorr.avg %>%
  mutate(Gene.symbol = ifelse(Gene.symbol == "AK6///TAF9", "TAF9", Gene.symbol))%>%
  mutate(Gene.symbol = ifelse(Gene.symbol == "TGIF2-C20orf24///C20orf24", "TGIF2-C20orf24", Gene.symbol))%>%
  filter(Gene.symbol != "TAF9B")

intersect(unique_top20_row_names,gset.expression.matching_topcorr.avg$Gene.symbol)

# alpha=0 means ridge, alpha=1 means lasso
set.seed(2000)
rownames(gset.expression.matching_topcorr.avg) <- gset.expression.matching_topcorr.avg$Gene.symbol
dim(t(gset.expression.matching_topcorr.avg[,-1]))
colnames(t(gset.expression.matching_topcorr.avg[,-1])) #<- gset.expression.matching_topcorr.avg$Gene.symbol

# Set row names to the 'Gene.symbol' column
#rownames(gset.expression.matching_topcorr.avg) <- gset.expression.matching_topcorr.avg$Gene.symbol

# Remove the 'Gene.symbol' column to keep only expression values
#gset.expression.matching_topcorr.avg <- gset.expression.matching_topcorr.avg[,-1]

# Transpose the matrix
transposed_matrix <- t(gset.expression.matching_topcorr.avg[,-1])

# Check the dimensions and column names of the transposed matrix
dim(transposed_matrix)
colnames(transposed_matrix) <- gset.expression.matching_topcorr.avg$Gene.symbol


glmnet.age.samegene<-cv.glmnet(transposed_matrix,gset.pheno$age,alpha=1,family="gaussian")

#glmnet.age.samegene <- cv.glmnet(t(gset.expression.matching_topcorr.avg[,-1]), gset.pheno$age, alpha = 0, family = "gaussian", lambda = seq(0.001, 2, by = 0.01))

glmnet.age.samegene$lambda.1se #use 1.296205 to rerun the model
glmnet.age.samegene$lambda.min #0.7417361

glmnet.age.samegene.coeffs<-predict(glmnet.age.samegene,type="coefficients")
#gset.expression.matching_topcorr.avg$Gene.symbol
glmnet.age.samegene.coeffs@Dimnames[[1]][which(glmnet.age.samegene.coeffs!=0)]

#glmnet.age.samegene.coeffs
# Example of reducing lambda

# I run the 19 genes only but didn't find any coefficient

glmnet_predict <- predict(glmnet.age.samegene, newx=t(gset.expression.matching_topcorr.avg[,-1]), type="response")

library(tidyverse)
library(caret)

R2(glmnet_predict,gset.pheno$age) # 0.4711623
RMSE(glmnet_predict,gset.pheno$age) # 9.146458
#RMSE(glmnet_predict,phenos.df$age)/mean(phenos.df$age) #0.2412267
plot(gset.pheno$age,glmnet_predict)
ggplot(data=data.frame(age=gset.pheno$age,glmnet_age=c(glmnet_predict)),
       aes(x=age,y=glmnet_age)) + geom_point() + geom_smooth(method='lm', formula= y~x) + ggtitle("glmnet model of age")

#############
gset.data1 <- data.frame(age=gset.pheno$age, gene=c(glmnet_predict), sex = ifelse(gset.pheno$`gender:ch1`=="F", "Female", "Male"),
                         predicted = predict(lm(glmnet_predict[,1]~gset.pheno$age)),
                         res = residuals(lm(glmnet_predict[,1]~gset.pheno$age)),
                         res.sign = ifelse(residuals(lm(glmnet_predict[,1]~gset.pheno$age))>=0, 1, 0),
                         phenotype = if_else(gset.pheno$`subject group:ch1`=="CNTL; healthy control", "HC", "MDD"))
gset.data1$abres = abs(gset.data1$res)


#heteroskedasticity by plotting the residuals on the y axis and the predicted values on the x axis. 
plot(gset.data1$gene,gset.data1$res)
#install.packages("lmtest")
library(lmtest)
lmtest::bptest(fit.g)
car::ncvTest(fit.g)

ggplot(data = gset.data1, aes(y = res, x = gene)) + geom_point(col = 'blue') + geom_abline(slope = 0) + ggtitle("Test heteroskedasticity")
#par(mfrow=c(2,2))
plot(fit.g)
dev.off()
#and then test whether the residual is associated with MDD, like this
#t.test(gset.data1$phenotype~gset.data1$res)
#and add also see what happens when you adjust for age again
#glm(as.numeric(gset.data1$phenotype)~gset.data1$res+gset.data1$age) 

table(gset.data1$res.sign,gset.data1$phenotype)
prop.table(table(gset.data1$res.sign,gset.data1$phenotype), margin = 2)
table(gset.data1$res.sign,gset.data1$sex)
# histogram of the ages with a bin size of 1 or 2
h <- ggplot(gset.data1, aes(x=age)) + 
  geom_histogram(color="black", fill="pink",binwidth=1)+ 
  ggtitle("Histogram of the ages with a bin size of 1")  +
  theme(plot.title = element_text(size=20,face="bold"),
        axis.text=element_text(size=18),
        axis.title=element_text(size=18,face="bold"),
        legend.text = element_text(size=18),
        legend.title = element_text(size=18,face="bold"))+  
  stat_bin(binwidth=1, geom="text", aes(label=..count..), vjust=-1) 
h
table(gset.data1$age)

# Separate plot 

ggplot(gset.data1,
       aes(x=age,y=gene)) + 
  xlab("Chronological age (years)")+
  ylab("Biological age (years)")+
  geom_smooth(method='lm', formula= y~x) + 
  #geom_abline(intercept = 0, slope = 1,color="darkred") + 
  #ggtitle("glmnet model of age---by phenotype") +
  #geom_segment(aes(color = phenotype,xend = as.numeric(age), yend = as.numeric(predicted)), alpha = .5) +
  geom_point(aes(col=phenotype,shape=sex),size = 4) + #size = abs(res)
  guides(colour = guide_legend(override.aes = list(size=5))) + labs(fill = "Sex and phenotype")+
  geom_vline(aes(xintercept=20),colour="#BB0000",linetype="dashed")+
  geom_vline(aes(xintercept=30),colour="#BB0000",linetype="dashed")+
  geom_vline(aes(xintercept=40),colour="#BB0000",linetype="dashed")+
  theme_bw()+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 15, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma"),
        axis.text.x=element_text(colour="black", size = 11),
        axis.text.y=element_text(colour="black", size = 11),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text=element_text(size=18),
        axis.title=element_text(size=18,face="bold"),
        legend.text = element_text(size=18),
        legend.title = element_text(size=18,face="bold"))


table(gset.data1$phenotype, gset.data1$res.sign)
table(gset.data1$sex, gset.data1$res.sign)
# Ensure phenotype is a factor with exactly two levels
#data1$phenotype <- as.factor(data1$phenotype)
# Perform the t-test
t.test(res ~ phenotype, data = gset.data1)
# Test the significance of the correlation between x and y
cor.test(gset.data1$age,gset.data1$gene, method = "spearman")
print(spearman_test)
summary(lm(res ~ phenotype, data = data1))

# Density plots 
g.phenotype <- ggplot(data=gset.data1,
                      aes(x=res, colour=phenotype, fill=phenotype)) +
  geom_density(alpha = .3) + 
  ggtitle("A") +
  xlab(label = "L-GAGE")+
  theme_bw()+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 15, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma"),
        axis.text.x=element_text(colour="black", size = 11),
        axis.text.y=element_text(colour="black", size = 11),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"),
        legend.text = element_text(size=13),
        legend.title = element_text(size=13,face="bold"))

g.sex <- ggplot(data=gset.data1,
                aes(x=res, colour=sex, fill=sex)) +
  geom_density(alpha = .3) + 
  ggtitle("B") +
  xlab(label = "L-GAGE")+
  theme_bw()+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 15, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma"),
        axis.text.x=element_text(colour="black", size = 11),
        axis.text.y=element_text(colour="black", size = 11),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"),
        legend.text = element_text(size=13),
        legend.title = element_text(size=13,face="bold"))

ggarrange(g.phenotype, g.sex,
          #labels = c("(a)", "(b)", "(c)","(d)","(e)"),
          ncol = 2, nrow = 1,widths = 1,vjust = 1.5)

# Test the significance of the correlation between x and y
spearman_test <- cor.test(x, y, method = "spearman")
print(spearman_test)












###########################################################
summary(glmnet.age)
gset.expression.matching.avg.1<- data.frame(t(gset.expression.matching.avg[,-1]))
colnames(gset.expression.matching.avg.1) = gset.expression.matching.avg$Gene.symbol
ncol(gset.expression.matching.avg.1)


missing_cols <- setdiff(colnames(t(GxS.covfilter)), colnames(gset.expression.matching.avg.1))
gset.expression.matching.avg.1[missing_cols] <- 0  # or NA, depending on your needs
#gset.expression.matching.avg.1 <- scale(gset.expression.matching.avg.1)
#apply(gset.expression.matching.avg.1, 2, var)

# Get column names that differ
which(colnames(as.matrix(gset.expression.matching.avg.1)) != colnames(t(GxS.covfilter)))
# Assume matrix1 and matrix2 are your matrices
all(colnames(gset.expression.matching.avg.1) == colnames(t(GxS.covfilter)))

#setdiff(colnames(gset.expression.matching.avg.1), colnames(t(GxS.covfilter)))
#setdiff(colnames(t(GxS.covfilter)), colnames(gset.expression.matching.avg.1))

#glmnet.age$glmnet.fit
#colnames(as.matrix(gset.expression.matching.avg.1))
all(sort(colnames(gset.expression.matching.avg.1)) == sort(colnames(t(GxS.covfilter))))
# Reorder the columns of gset.expression.matching.avg.1
gset.expression.matching.avg.1 <- gset.expression.matching.avg.1[, colnames(t(GxS.covfilter))]
all(colnames(gset.expression.matching.avg.1) == colnames(t(GxS.covfilter)))

#predict(glmnet.age, newx = as.matrix(gset.expression.matching.avg.1), s = "lambda.1se")

#glmnet.age <- cv.glmnet(as.matrix(gset.expression.matching.avg.1), response, alpha = 0.5)
predict(glmnet.age, newx = as.matrix(gset.expression.matching.avg.1), s = "lambda.1se")
#coef(glmnet.age)


glmnet_predict <- predict(glmnet.age, newx = as.matrix(gset.expression.matching.avg.1), type = "response")
coef(glmnet.age, s = "lambda.1se")
glmnet.age.nonzero.coeffs.df <- data.frame(gene=glmnet.age.nonzero.coeffs[-1],
                                           coefficient=glmnet.age.coeffs[which(glmnet.age.coeffs!=0)][-1])


GxS.covfilter.samegene <- t(GxS.covfilter[glmnet.age.nonzero.coeffs.df$gene,])
gset.expression.matching.avg.1.same



