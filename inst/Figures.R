library(ggpubr)
# histogram of the ages with a bin size of 1 or 2
ggplot(data1, aes(x=age)) + 
  geom_histogram(aes(fill=phenotype), color="black", binwidth=1)+ 
  #ggtitle("Histogram of the ages with a bin size of 1")  +
  theme_bw() +
  theme(#plot.title = element_text(size=20,face="bold"),
        axis.text=element_text(size=18),
        axis.title=element_text(size=18,face="bold"),
        legend.text = element_text(size=18),
        legend.title = element_text(size=18,face="bold"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        #plot.title = element_text(size = 15, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma"),
        axis.text.x=element_text(colour="black", size = 11),
        axis.text.y=element_text(colour="black", size = 11),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15))+
  scale_x_continuous(name = "Age (years)") +
  scale_y_continuous(name = "Count")
  #stat_bin(binwidth=1, geom="text",  aes(label=..count..),vjust=-1) #aes(label=..count..),

library(ggpubr)

LR_plot <- lm.results.sorted %>%
  filter( p.val<= 0.05 ) %>%
  ggplot( aes(x=p.adj)) +
  geom_histogram( binwidth=0.05, fill="#69b3a2", color="#e9ecef", alpha=0.9,boundary=0) +
  #ggtitle("Age associated genes selected by linear regression") +
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
        plot.title = element_text(size = 15, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma"),
        axis.text.x=element_text(colour="black", size = 11),
        axis.text.y=element_text(colour="black", size = 11),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15))+ ggtitle("A")

LR_plot

glmnet_plot <- glmnet.age.nonzero.coeffs.df %>%
  ggplot( aes(x=coefficient)) +
  geom_histogram( binwidth=5, fill="#69b3a2", color="#e9ecef", alpha=0.9,boundary=-20) +
  #ggtitle("Age associated genes selected by Lasso") +
  theme_ipsum() +
  theme(plot.title = element_text(size=15))+  
  stat_bin(boundary=-20, binwidth=5, geom="text", aes(label=..count..), vjust=2)+
  geom_vline(xintercept = 0, size = 1, colour = "#FF3721",
             linetype = "dashed")+
  scale_x_continuous(name = "Coefficient",
                     breaks = seq(-10, 15, 5)) +
  scale_y_continuous(name = "Count")+
  theme_bw() +
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 15, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma"),
        axis.text.x=element_text(colour="black", size = 11),
        axis.text.y=element_text(colour="black", size = 11),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15))+ ggtitle("B")

glmnet_plot


NPDR_plot <- top.p05.npdr %>%
  filter( pval.adj<.05 ) %>%
  ggplot( aes(x=pval.adj)) +
  geom_histogram( binwidth=0.01, fill="#69b3a2", color="#e9ecef", alpha=0.9,boundary=0) +
  #ggtitle("Age associated genes selected by NPDR") +
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
        plot.title = element_text(size = 15, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma"),
        axis.text.x=element_text(colour="black", size = 11),
        axis.text.y=element_text(colour="black", size = 11),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15))+ ggtitle("C")
NPDR_plot

LassoNPDR_plot <- npdr_lasso_results %>%
  filter( scores>0 ) %>%
  ggplot( aes(x=scores)) +
  geom_histogram( binwidth=1, fill="#69b3a2", color="#e9ecef", alpha=0.9,boundary=-20) +
  #ggtitle("Age associated genes selected by LASSO NPDR") +
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
        plot.title = element_text(size = 15, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma"),
        axis.text.x=element_text(colour="black", size = 11),
        axis.text.y=element_text(colour="black", size = 11),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15))+ ggtitle("D")
LassoNPDR_plot

## Create a plot of variable importance
RF_plot <- var_imp_top20 %>%
  
  ## Sort the data by importance
  arrange(importance) %>%
  
  ## Create a ggplot object for aesthetic
  ggplot(aes(x=reorder(variables, importance), y=importance)) + 
  
  ## Plot the bar graph
  geom_bar(stat='identity') + 
  
  ## Flip the graph to make a horizontal bar plot
  coord_flip() + 
  
  ## Add x-axis label
  xlab('Variables') +
  ylab('Importance') +
  
  ## Add a title
  #labs(title='Top 20 genes in Random forest') + 
  
  ## Some layout for the plot
  theme_minimal() + 
  theme(axis.text = element_text(size = 11), 
        #axis.title = element_text(size = 12), 
        plot.title = element_text(size = 15, family = "Tahoma", face = "bold"),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)
  )+ ggtitle("E")
RF_plot

ggarrange(LR_plot, glmnet_plot, NPDR_plot, LassoNPDR_plot,#RF_plot, 
          #labels = c("(a)", "(b)", "(c)","(d)","(e)"),
          ncol = 2, nrow = 2,widths = 1,vjust = 1.5)


### DID ########################
#〖𝑅𝑒𝑠𝑖𝑑𝑢𝑎𝑙𝑠〗_𝑖=𝛽_0+𝛽_1 〖𝑀𝐷𝐷〗_𝑖+𝛽_2 〖𝑆𝑒𝑥〗_𝑖+𝛽_3 〖𝑀𝐷𝐷〗_𝑖∗〖𝑆𝑒𝑥〗_𝑖+𝜇

#### DID age ##################################################
data1$old <- ifelse(data1$age<40,"<40",">=40")
data1$old <- factor(data1$old , levels=c("<40",">=40"))
did.age <- summary(lm(data1$res~data1$phenotype+data1$old+data1$phenotype*data1$old))

A <- sum(did.age$coefficients[,1])
B <- sum(did.age$coefficients[1,1], did.age$coefficients[2,1])
C <- sum(did.age$coefficients[1,1], did.age$coefficients[3,1])
D <- did.age$coefficients[1,1]
A
B
C
D

B-D
A-C
(B-D) - (A-C)

d1 <- data.frame(`Age>=40`=c(A,C), `Age<40`=c(B,D),row.names = c("MDD","HC"))
d1 <- data.frame(value=c(A,C,B,D),Age = c(">=40",">=40","<40","<40"),
                 MDD = c("MDD","HC","MDD","HC"),
                 label1 = c("Beta0+beta1+beta2+beta3","Beta0+beta2","Beta0+beta1","Beta0"))

library(ggplot2)
library(ggrepel)
DID_age <- ggplot(d1, aes(x=Age, y=value, colour=MDD, group=MDD)) +
  geom_line() +
  geom_point() + 
  theme_minimal() + 
  ylab(label = 'L-GAGE') +
  xlab(label = 'Age (years)')+
  ggtitle(label = "A")+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 15, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma"),
        axis.text.x=element_text(colour="black", size = 13),
        axis.text.y=element_text(colour="black", size = 13),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.title=element_text(size=15,face="bold"),
        legend.text = element_text(size=13),
        legend.title = element_text(size=13,face="bold"))+
  geom_label_repel(aes(label = d1$label1),
                   box.padding   = 1, 
                   point.padding = 1,
                   segment.color = 'grey50')

DID_age


box_age <- ggboxplot(data = data1, x = "old",y = "res",
                     combine = TRUE,
                     #color = "phenotype", #palette = c("#00AFBB", "#E7B800"),
                     color = "phenotype",
                     ylab = "L-GAGE", 
                     add = "jitter",                              # Add jittered points
                     add.params = list(size = 1, jitter = 0.1))+  # Point size and the amount of jittering
  ggtitle("B")+
  xlab("Age (years)")+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 15, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma"),
        axis.text.x=element_text(colour="black", size = 13),
        axis.text.y=element_text(colour="black", size = 13),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.title=element_text(size=15,face="bold"),
        legend.position="right",
        legend.text = element_text(size=13),
        legend.title = element_text(size=13,face="bold"))


ggarrange(DID_age, box_age,
          #labels = c("(a)", "(b)", "(c)","(d)","(e)"),
          ncol = 2, nrow = 1,widths = 1,vjust = 1.5)

######## sex #########################
#### DID sex ##################################################
### multiple linear model, DID gene sex
DID.sex <- glm(res~phenotype+sex+phenotype*sex, data = data1)
DID.sex
library(modelsummary)
msummary(DID.sex, stars = c('*' = .1, '**' = .05, '***' = .01),
         statistic = "p.value", notes = "Coefficient p values are added in parentheses")

### test of H_0: all regression coefficients are zero 
### (ignore intercept)

### define coefficients of linear function directly
K <- diag(length(coef(DID.sex)))[-1,]
rownames(K) <- names(coef(DID.sex))[-1]
K

### set up general linear hypothesis
#install.packsexs("multcomp")
library(multcomp)
glht(DID.sex, linfct = K)

### define coefficients of linear function directly
#K <- diag(length(coef(DID.sex)))[-1,]
#rownames(K) <- names(coef(DID.sex))[-1]
K <- rbind("phenotype=HC   sex=Female"=c(1,0,0,0),
           "phenotype=MDD  sex=Female"=c(1,1,0,0),
           "phenotype=HC   sex=Male"=c(1,0,1,0),
           "phenotype=MDD  sex=Male"=c(1,1,1,1))
did.result <- glht(DID.sex, linfct = K)
summary(did.result)

#data1$sex.num <- ifelse(data1$sex=="Female", 1, 0)
data1$sex.num <- ifelse(data1$sex=="Male", 1, 0)
data1$mdd.num <- ifelse(data1$phenotype=="MDD", 1, 0)

did.sex <- summary(glm(data1$res~data1$phenotype+data1$sex+data1$phenotype*data1$sex))
#did.sex <- summary(lm(data1$res~data1$mdd.num+data1$sex.num+data1$mdd.num*data1$sex.num))

did.sex

datamale = subset(data1,sex=="Male")
boxplot(res~phenotype,data = subset(data1,sex=="Male"))
boxplot(res~phenotype,data = subset(data1,sex=="Female"))

summary(lm(data1$res~data1$phenotype+data1$sex))
summary(lm(res~phenotype, data = subset(data1,sex=="Male")))
summary(lm(res~phenotype, data = subset(data1,sex=="Female")))

boxplot(res~sex,data = subset(data1,phenotype=="HC"))
boxplot(res~sex,data = subset(data1,phenotype=="MDD"))

summary(lm(res~sex, data = subset(data1,phenotype=="HC")))
summary(lm(res~sex, data = subset(data1,phenotype=="MDD")))
cor(data1$res,data1$sex.num)
summary(lm(res~sex.num, data = data1))

#plot(data1$phenotype,data1$res)
#corrplot(cor(data.frame(data1$age,data1$res,data1$sex.num)))
ggboxplot(data = subset(data1,sex=="Female"), x = "phenotype",y = "res",
          combine = TRUE,
          #color = "phenotype", #palette = c("#00AFBB", "#E7B800"),
          color = "phenotype",
          ylab = "L-GAGE", 
          add = "jitter",                              # Add jittered points
          add.params = list(size = 1, jitter = 0.1))

ggboxplot(data = subset(data1,phenotype=="MDD"), x = "sex",y = "res",
          combine = TRUE,
          #color = "phenotype", #palette = c("#00AFBB", "#E7B800"),
          color = "sex",
          ylab = "L-GAGE", 
          add = "jitter",                              # Add jittered points
          add.params = list(size = 1, jitter = 0.1))

#Female <- ifelse(data1$sex=="Female", 1, 0)
#did.sex <- summary(glm(data1$res~data1$phenotype+Female+data1$phenotype*Female))


A <- sum(did.sex$coefficients[,1]) # Beta0+beta1+beta2+beta3
B <- sum(did.sex$coefficients[1,1], did.sex$coefficients[2,1])
C <- sum(did.sex$coefficients[1,1], did.sex$coefficients[3,1])
D <- did.sex$coefficients[1,1]
A
B
C
D

B-D
A-C
(B-D) - (A-C)

d1 <- data.frame(`Male`=c(A,C), `Female`=c(B,D),row.names = c("MDD","HC"))
d1 <- data.frame(sex = c("Male","Male","Female","Female"), value=c(A,C,B,D),
                 MDD = c("MDD","HC","MDD","HC"),
                 label1 = c("Beta0+beta1+beta2+beta3","Beta0+beta2","Beta0+beta1","Beta0"))
d1

levels(data1$sex)
as.factor(data1$phenotype)
library(ggplot2)
library(ggrepel)
DID_sex <- ggplot(d1, aes(x=sex, y=value, colour=MDD, group=MDD)) +
  geom_line() +
  geom_point() + 
  theme_minimal() + 
  ylab(label = 'L-GAGE') +
  ggtitle(label = "A")+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 15, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma"),
        axis.text.x=element_text(colour="black", size = 13),
        axis.text.y=element_text(colour="black", size = 13),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.text = element_text(size=13),
        legend.title = element_text(size=13,face="bold"),
        axis.title=element_text(size=15,face="bold"))+
  geom_label_repel(aes(label = d1$label1),
                   box.padding   = 1, 
                   point.padding = 1,
                   segment.color = 'grey50')
DID_sex


box_sex <- ggboxplot(data = data1, x = "sex",y = "res",
                     combine = TRUE,
                     #color = "phenotype", #palette = c("#00AFBB", "#E7B800"),
                     color = "phenotype",
                     ylab = "L-GAGE", 
                     add = "jitter",                              # Add jittered points
                     add.params = list(size = 1, jitter = 0.1))+  # Point size and the amount of jittering
  ggtitle("B")+
  xlab("sex")+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 15, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma"),
        axis.text.x=element_text(colour="black", size = 13),
        axis.text.y=element_text(colour="black", size = 13),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.position="right",
        legend.text = element_text(size=13),
        legend.title = element_text(size=13,face="bold"),
        axis.title=element_text(size=15,face="bold"))


ggarrange(DID_sex, box_sex,
          #labels = c("(a)", "(b)", "(c)","(d)","(e)"),
          ncol = 2, nrow = 1,widths = 1,vjust = 1.5)


lm.sex <- summary(glm(data1$res~data1$phenotype+data1$sex))
lm.sex

A <- sum(lm.sex$coefficients[,1]) # Beta0+beta1+beta2+beta3
B <- sum(lm.sex$coefficients[1,1], lm.sex$coefficients[2,1])
C <- sum(lm.sex$coefficients[1,1], lm.sex$coefficients[3,1])
D <- lm.sex$coefficients[1,1]
A
B
C
D
d1 <- data.frame(sex = c("Male","Male","Female","Female"), value=c(A,C,B,D),
                 MDD = c("MDD","HC","MDD","HC"),
                 label1 = c("Beta0+beta1+beta2","Beta0+beta2","Beta0+beta1","Beta0"))
d1

# =============================================================
# use lambda.1se(1.636048) to rerun the glmnet model
# =============================================================
#install.packages("glmnet")
library(glmnet)
# alpha=0 means ridge, alpha=1 means lasso
set.seed(2000)
glmnet.age<-cv.glmnet(t(GxS.covfilter),phenos.df$age,alpha=1,family="gaussian", lambda = 1.636048)
glmnet.age<-cv.glmnet(t(GxS.covfilter),phenos.df$age,alpha=1,family="gaussian")

glmnet.age$lambda.1se #use 1.636048 to rerun the model
glmnet.age$lambda.min #1.127664

glmnet.age.coeffs<-predict(glmnet.age,type="coefficients", s = "lambda.1se")
plot(glmnet.age) # MSE be lowest when log lambda around 0
glmnet.age.nonzero.coeffs <- glmnet.age.coeffs@Dimnames[[1]][which(glmnet.age.coeffs!=0)]
write.table(glmnet.age.nonzero.coeffs,row.names=F,col.names=F,quote=F)
length(glmnet.age.nonzero.coeffs)



glmnet.age.coeffs<-predict(glmnet.age,type="coefficients",s = "lambda.min")
plot(glmnet.age) # MSE be lowest when log lambda around 0
glmnet.age.nonzero.coeffs <- glmnet.age.coeffs@Dimnames[[1]][which(glmnet.age.coeffs!=0)]
write.table(glmnet.age.nonzero.coeffs,row.names=F,col.names=F,quote=F)
length(glmnet.age.nonzero.coeffs)

#intersect(glmnet.age.nonzero.coeffs, lm.results.sorted.0.05$gene)
#Add coefficient column. 
glmnet.age.nonzero.coeffs.df <- data.frame(gene=glmnet.age.nonzero.coeffs,
                                           coefficient=glmnet.age.coeffs[which(glmnet.age.coeffs!=0)])

glmnet.age$lambda.min
summary(glmnet.age$cvm)
summary(glmnet.age$cvup)
summary(glmnet.age$lambda)
# genes in common with glmnet and lm p.adjust
#intersect(glmnet.age.nonzero.coeffs, lm.results.sorted.0.05$gene)
glmnet_predict <- predict(glmnet.age, newx=t(GxS.covfilter), type="response",s = "lambda.min")
#install.packages("tidyverse")
library(tidyverse)
library(caret)
R2(glmnet_predict,phenos.df$age) # 0.6014565
RMSE(glmnet_predict,phenos.df$age) # 7.679305
RMSE(glmnet_predict,phenos.df$age)/mean(phenos.df$age) #0.2412267
#plot(phenos.df$age,glmnet_predict)
ggplot(data=data.frame(age=phenos.df$age,glmnet_age=glmnet_predict[,1]),
       aes(x=age,y=glmnet_age)) + geom_point() + geom_smooth(method='lm', formula= y~x) + ggtitle("glmnet model of age")

fit.g <- lm(glmnet_predict[,1]~phenos.df$age)
summary(fit.g)
res = residuals(fit.g)
plot(res, type = "l")
plot(sort(res), type = "l")

data1 <- data.frame(age=phenos.df$age, gene=glmnet_predict[,1], sex = phenos.df$sex,
                    predicted = predict(lm(glmnet_predict[,1]~phenos.df$age)),
                    res = residuals(lm(glmnet_predict[,1]~phenos.df$age)),
                    res.sign = ifelse(res>=0, 1, 0),
                    phenotype = colnames(GxS.covfilter))
data1$abres = abs(data1$res)
#heteroskedasticity by plotting the residuals on the y axis and the predicted values on the x axis. 
plot(data1$gene,data1$res)
#install.packages("lmtest")
library(lmtest)
lmtest::bptest(fit.g)
car::ncvTest(fit.g)

ggplot(data = data1, aes(y = res, x = gene)) + geom_point(col = 'blue') + geom_abline(slope = 0) + ggtitle("Test heteroskedasticity")
#par(mfrow=c(2,2))
plot(fit.g)
dev.off()
#and then test whether the residual is associated with MDD, like this
#t.test(data1$phenotype~data1$res)
#and add also see what happens when you adjust for age again
#glm(as.numeric(data1$phenotype)~data1$res+data1$age) 

table(data1$res.sign,data1$phenotype)
prop.table(table(data1$res.sign,data1$phenotype), margin = 2)
table(data1$res.sign,data1$sex)
# histogram of the ages with a bin size of 1 or 2
h <- ggplot(data1, aes(x=age)) + 
  geom_histogram(color="black", fill="pink",binwidth=1)+ 
  ggtitle("Histogram of the ages with a bin size of 1")  +
  theme(plot.title = element_text(size=20,face="bold"),
        axis.text=element_text(size=18),
        axis.title=element_text(size=18,face="bold"),
        legend.text = element_text(size=18),
        legend.title = element_text(size=18,face="bold"))+  
  stat_bin(binwidth=1, geom="text", aes(label=..count..), vjust=-1) 
h
table(data1$age)

# Separate plot 

ggplot(data=data.frame(age=phenos.df$age, gene=glmnet_predict[,1], sex = phenos.df$sex,
                       predicted = predict(lm(glmnet_predict[,1]~phenos.df$age)),
                       res = glmnet_predict[,1]-predict(lm(glmnet_predict[,1]~phenos.df$age)),
                       phenotype = colnames(GxS.covfilter)),
       aes(x=age,y=gene)) + 
  xlab("chronological age")+
  ylab("Prediction age")+
  geom_smooth(method='lm', formula= y~x) + 
  geom_abline(intercept = 0, slope = 1,color="darkred") + 
  ggtitle("glmnet model of age---by sex") +
  geom_segment(aes(color = sex,xend = as.numeric(age), yend = as.numeric(predicted)), alpha = .5) +
  geom_point(aes(col=sex,size = abs(res))) +
  theme(plot.title = element_text(size=20,face="bold"),
        axis.text=element_text(size=18),
        axis.title=element_text(size=18,face="bold"),
        legend.text = element_text(size=18),
        legend.title = element_text(size=18,face="bold")) +
  guides(colour = guide_legend(override.aes = list(size=5))) + labs(fill = "Sex and phenotype")



ggplot(data=data.frame(age=phenos.df$age, gene=glmnet_predict[,1], sex = phenos.df$sex,
                       predicted = predict(lm(glmnet_predict[,1]~phenos.df$age)),
                       res = glmnet_predict[,1]-predict(lm(glmnet_predict[,1]~phenos.df$age)),
                       phenotype = colnames(GxS.covfilter)),
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

# Density plots 
g.phenotype <- ggplot(data=data.frame(age=phenos.df$age, gene=glmnet_predict[,1], sex = phenos.df$sex,
                                      predicted = predict(lm(glmnet_predict[,1]~phenos.df$age)),
                                      res = glmnet_predict[,1]-predict(lm(glmnet_predict[,1]~phenos.df$age)),
                                      phenotype = colnames(GxS.covfilter)),
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

g.sex <- ggplot(data=data.frame(age=phenos.df$age, gene=glmnet_predict[,1], sex = phenos.df$sex,
                                predicted = predict(lm(glmnet_predict[,1]~phenos.df$age)),
                                res = glmnet_predict[,1]-predict(lm(glmnet_predict[,1]~phenos.df$age)),
                                phenotype = colnames(GxS.covfilter)),
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
