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

hist(glmnet.age.nonzero.coeffs.df$coefficient)

glmnet_plot <- glmnet.age.nonzero.coeffs.df %>%
  ggplot( aes(x=coefficient)) +
  geom_histogram( binwidth=5, fill="#69b3a2", color="#e9ecef", alpha=0.9,boundary=-20) +
  ggtitle("21 age associated genes selected by Lasso with non-zero coefficient") +
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
        plot.title = element_text(size = 13, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma"),
        axis.text.x=element_text(colour="black", size = 13),
        axis.text.y=element_text(colour="black", size = 13))

glmnet_plot

glmnet_predict <- predict(glmnet.age, newx=t(GxS.covfilter), type="response")

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
t.test(glmnet_predict[,1],phenos.df$age)

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
  geom_point(aes(col=sex,size = 0.2)) +
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
  xlab("Chronological age")+
  ylab("Biological age")+
  geom_smooth(method='lm', formula= y~x) + 
  #geom_abline(intercept = 0, slope = 1,color="darkred") + 
  ggtitle("glmnet model of age---by phenotype") +
  #geom_segment(aes(color = phenotype,xend = as.numeric(age), yend = as.numeric(predicted)), alpha = .5) +
  geom_point(aes(col=phenotype,shape=sex),size = 4) + #size = abs(res)
  theme(plot.title = element_text(size=20,face="bold"),
        axis.text=element_text(size=18),
        axis.title=element_text(size=18,face="bold"),
        legend.text = element_text(size=18),
        legend.title = element_text(size=18,face="bold")) +
  guides(colour = guide_legend(override.aes = list(size=5))) + labs(fill = "Sex and phenotype")+
  geom_vline(aes(xintercept=20),colour="#BB0000",linetype="dashed")+
  geom_vline(aes(xintercept=30),colour="#BB0000",linetype="dashed")+
  geom_vline(aes(xintercept=40),colour="#BB0000",linetype="dashed")


# Density plots 
g.phenotype <- ggplot(data=data.frame(age=phenos.df$age, gene=glmnet_predict[,1], sex = phenos.df$sex,
                                      predicted = predict(lm(glmnet_predict[,1]~phenos.df$age)),
                                      res = glmnet_predict[,1]-predict(lm(glmnet_predict[,1]~phenos.df$age)),
                                      phenotype = colnames(GxS.covfilter)),
                      aes(x=res, colour=phenotype, fill=phenotype)) +
  geom_density(alpha = .3) + 
  ggtitle("glmnet model of age") +
  xlab(label = "glmnet Gene age ~ Chronological age (Residuals)")+
  theme(plot.title = element_text(size=20,face="bold"),
        axis.text=element_text(size=18),
        axis.title=element_text(size=18,face="bold"),
        legend.text = element_text(size=18),
        legend.title = element_text(size=18,face="bold"))

g.sex <- ggplot(data=data.frame(age=phenos.df$age, gene=glmnet_predict[,1], sex = phenos.df$sex,
                                predicted = predict(lm(glmnet_predict[,1]~phenos.df$age)),
                                res = glmnet_predict[,1]-predict(lm(glmnet_predict[,1]~phenos.df$age)),
                                phenotype = colnames(GxS.covfilter)),
                aes(x=res, colour=sex, fill=sex)) +
  geom_density(alpha = .3) + 
  ggtitle("glmnet model of age") +
  xlab(label = "glmnet Gene age ~ Chronological age (Residuals)")+
  theme(plot.title = element_text(size=20,face="bold"),
        axis.text=element_text(size=18),
        axis.title=element_text(size=18,face="bold"),
        legend.text = element_text(size=18),
        legend.title = element_text(size=18,face="bold"))

ggplot(data1,
       aes(x=abres, colour=phenotype, fill=phenotype)) +
  geom_density(alpha = .3) + 
  ggtitle("glmnet model of age") +
  xlab(label = "glmnet Gene age ~ Chronological age (Residuals)")+
  theme(plot.title = element_text(size=20,face="bold"),
        axis.text=element_text(size=18),
        axis.title=element_text(size=18,face="bold"),
        legend.text = element_text(size=18),
        legend.title = element_text(size=18,face="bold"))

library(ggpubr)
ggarrange(g.phenotype,g.sex)

# Spearman rank correlation coefficient
cor.test(glmnet_predict[,1],phenos.df$age,  method = "spearman",exact=FALSE)

#t.test(data1$phenotype~data1$predicted)

#fit.g <- lm(glmnet_predict[,1]~phenos.df$age)

summary(fit.g)
#res = residuals(fit.g)
#abres = abs(data1$res)
cor.test(res, as.numeric(mddPheno))

#res = glmnet_predict[,1]-predict(lm(glmnet_predict[,1]~phenos.df$age))
t.test(phenos.df$age~as.factor(mddPheno)) #p-value = 0.1627
t.test(data1$predicted~as.factor(mddPheno)) #p-value = 0.1627
t.test(data1$res~as.factor(mddPheno)) #p-value = 0.1627
t.test(data1$res~data1$phenotype) #p-value = 0.1627

#data1 <- data.frame(age=phenos.df$age, gene=glmnet_predict[,1], sex = phenos.df$sex,
# predicted = predict(lm(glmnet_predict[,1]~phenos.df$age)),
# res = glmnet_predict[,1]-predict(lm(glmnet_predict[,1]~phenos.df$age)),
#  res.sign = ifelse(res>=0, 1, 0),
# phenotype = colnames(GxS.covfilter), abres = abs(res))
#pairs(data.frame(chr_age = data1_2$age, bio_age = data1_2$gene, residuals = data1_2$res))
#pairs(data.frame(data1$res,data1$age, as.numeric(data1$sex), as.numeric(data1$phenotype)))

############## test residuals ##############################
ggboxplot(data = data1, x = "phenotype", y = "res", 
          color = "phenotype", palette = c("#00AFBB", "#E7B800"),
          xlab = "phenotype", ylab = "Residuals") +
  geom_jitter(shape=16, position=position_jitter(0.2))+ 
  ggtitle("Boxplots of residuals in HC and MDD")

ggboxplot(data = data1, x = "phenotype", y = "abres", 
          color = "phenotype", palette = c("#00AFBB", "#E7B800"),
          xlab = "phenotype", ylab = "absolute Residuals") +
  geom_jitter(shape=16, position=position_jitter(0.2))+ 
  ggtitle("Boxplots of residuals in HC and MDD")

data1$age40 <- ifelse(data1$age<40,0,1)
ggboxplot(data = data1, x = "age40", y = "res", 
          color = "phenotype", palette = c("#00AFBB", "#E7B800"),
          xlab = "age40", ylab = "Residuals") +
  geom_jitter(shape=16, position=position_jitter(0.2))+ 
  ggtitle("Boxplots of residuals in HC and MDD --- by Age 40")

ggboxplot(data = data1, x = "sex", y = "res", 
          color = "phenotype", palette = c("#00AFBB", "#E7B800"),
          xlab = "sex", ylab = "Residuals") +
  geom_jitter(shape=16, position=position_jitter(0.2))+ 
  ggtitle("Boxplots of residuals in HC and MDD --- by Sex")

library(ggpubr)
data1$age40 <- ifelse(data1$age >=40, "Age>=40", "Age<40")
data1$age40 <- ifelse(data1$age<40,0,1)
ggboxplot(data = data1, x = "age40",y = "res",
          combine = TRUE,
          color = "phenotype", palette = c("#00AFBB", "#E7B800"),
          ylab = "Residuals(GAGE)", 
          add = "jitter",                              # Add jittered points
          add.params = list(size = 1, jitter = 0.1))+  # Point size and the amount of jittering
  ggtitle("Boxplots of GAGE in HC and MDD --- by Age (Lasso)")+
  xlab("Age")

ggboxplot(data = data1, x = "sex",
          y = c("res"),
          combine = TRUE,
          color = "phenotype", palette = c("#00AFBB", "#E7B800"),
          ylab = "Residuals", 
          add = "jitter",                              # Add jittered points
          add.params = list(size = 1, jitter = 0.1))+  # Point size and the amount of jittering
  ggtitle("Boxplots of GAGE in HC and MDD --- by Sex (Lasso)")

t.test(abres~phenotype, data = data1)
t.test(res~phenotype, data = data1)
cor.test(data1$res, as.numeric(mddPheno))


#glm(phenoMDD~residual)
summary(glm(as.factor(mddPheno) ~ res, data = data1, family=binomial))

#glm(phenoMDD~residual+phenos.df$age) 
summary(glm(as.factor(mddPheno) ~ age+res,data = data1, family=binomial))

#glm(as.factor(mddPheno) ~ age + res + sex
summary(glm(as.factor(mddPheno) ~ age + res + sex, data = data1, family=binomial))

#pairs(data1$res, as.factor(mddPheno))#, ,as.factor(data1$sex),data1$age)

#we tested the residuals for MDD on subjects whose chronological age was greater than 40 years
data1_1 <- subset(data1, age>=40)
t.test(res~phenotype, data = data1_1)
#t.test(abres~phenotype, data = data1_1) very not significant  

b1 <- ggboxplot(data = subset(data1, age<=20), x = "phenotype", y = "res", 
                color = "phenotype", palette = c("#00AFBB", "#E7B800"),
                xlab = "phenotype", ylab = "Residuals") +
  geom_jitter(shape=16, position=position_jitter(0.2))+ 
  ggtitle("Residuals in HC and MDD (Age <= 20)")

b2 <- ggboxplot(data = subset(data1, age>20&age<=30), x = "phenotype", y = "res", 
                color = "phenotype", palette = c("#00AFBB", "#E7B800"),
                xlab = "phenotype", ylab = "Residuals") +
  geom_jitter(shape=16, position=position_jitter(0.2))+ 
  ggtitle("Residuals in HC and MDD (Age>20 & Age<=30)")

b3 <- ggboxplot(data = subset(data1, age>30&age<=40), x = "phenotype", y = "res", 
                color = "phenotype", palette = c("#00AFBB", "#E7B800"),
                xlab = "phenotype", ylab = "Residuals") +
  geom_jitter(shape=16, position=position_jitter(0.2))+ 
  ggtitle("Residuals in HC and MDD (Age>30 & Age<40)")

b4 <- ggboxplot(data = subset(data1, age>=40), x = "phenotype", y = "res", 
                color = "phenotype", palette = c("#00AFBB", "#E7B800"),
                xlab = "phenotype", ylab = "Residuals") +
  geom_jitter(shape=16, position=position_jitter(0.2))+ 
  ggtitle("Rresiduals in HC and MDD (Age >= 40)")

ggarrange(b1,b2,b3, b4)

# t test for subset dataset 
data1_1 <- subset(data1, age>=40)
data1_1$old <- ifelse(data1_1$age>=40, 1, 0)
t.test(res~phenotype, data = subset(data1, age<=20))
#t.test(res~phenotype, data = subset(data1, age==24))
t.test(res~phenotype, data = subset(data1, age>20&age<=30))
t.test(res~phenotype, data = subset(data1, age>30&age<40))
t.test(res~phenotype, data = subset(data1, age>=40))

#####<=20#########
summary(glm(as.factor(phenotype) ~ res, data = subset(data1, age<=20), family=binomial))

#glm(phenoMDD~residual+phenos.df$age) 
summary(glm(as.factor(phenotype) ~ age+res,data = subset(data1, age<=20), family=binomial))

#glm(as.factor(mddPheno) ~ age + res + sex
summary(glm(as.factor(phenotype) ~  age + res + sex, data = subset(data1, age<=20), family=binomial))


#####<=30 ########
summary(glm(as.factor(phenotype) ~ res, data = subset(data1, age<=30), family=binomial))

#glm(phenoMDD~residual+phenos.df$age) 
summary(glm(as.factor(phenotype) ~ age+res,data = subset(data1, age<=30), family=binomial))

#glm(as.factor(mddPheno) ~ age + res + sex
summary(glm(as.factor(phenotype) ~  age + res + sex, data = subset(data1, age<=30), family=binomial))


#####<=40 #####
summary(glm(as.factor(phenotype) ~ res, data = subset(data1, age<=40), family=binomial))

#glm(phenoMDD~residual+phenos.df$age) 
summary(glm(as.factor(phenotype) ~ age+res,data = subset(data1, age<=40), family=binomial))

#glm(as.factor(mddPheno) ~ age + res + sex
summary(glm(as.factor(phenotype) ~  age + res + sex, data = subset(data1, age<=40), family=binomial))


#####>=30 ########
summary(glm(as.factor(phenotype) ~ res, data = subset(data1, age>=30), family=binomial))

#glm(phenoMDD~residual+phenos.df$age) 
summary(glm(as.factor(phenotype) ~ age+res,data = subset(data1, age>=30), family=binomial))

#glm(as.factor(mddPheno) ~ age + res + sex
summary(glm(as.factor(phenotype) ~  age + res + sex, data = subset(data1, age>=30), family=binomial))


#####>=40
#glm(phenoMDD~residual)
summary(glm(as.factor(phenotype) ~ res, data = subset(data1, age>=40), family=binomial))

#glm(phenoMDD~residual+phenos.df$age) 
summary(glm(as.factor(phenotype) ~ age+res,data = subset(data1, age>=40), family=binomial))

#glm(as.factor(mddPheno) ~ age + res + sex
summary(glm(as.factor(phenotype) ~  age + res + sex, data = subset(data1, age>=40), family=binomial))

summary(glm(as.factor(phenotype) ~  age40 + res + sex, data = data1, family=binomial))

summary(lm(data1$res~data1$phenotype+data1$age+
             data1$phenotype*data1$age))

summary(lm(data1$res~data1$phenotype))

summary(lm(data1$res~data1$phenotype+data1$age40+
             data1$phenotype*data1$age40+
             data1$sex))

table(data1$sex,data1$phenotype)
pairs(data.frame(data1_1$age,data1_1$predicted,data1_1$res))
# replace the residual to the residual sign (dummy variable)
lm4 <- glm(as.factor(mddPheno) ~ data1$res.sign, family=binomial)
summary(lm4)
t.test(data1$res.sign~mddPheno)

chisq.test(data1$res.sign,as.factor(mddPheno))
chisq.test(data1$res,as.factor(mddPheno))

# test if residuals still have correlation with age
t.test(res, data1$age, paired = TRUE, alternative = "two.sided")
summary(lm(res~data1$age))

t.test(age~phenotype, data = data1)
#t.test(as.factor(sex)~phenotype, data = data1)
#t.test(data1$sex, data1$phenotype)

chisq.test(data1$sex,as.factor(mddPheno))
chisq.test(data1$age,as.factor(mddPheno))
chisq.test(data1$age40,as.factor(mddPheno))

##################################################
library(modelsummary)
data1$GAGE <- data1$res
Age40.GAGE1 <- glm(as.factor(phenotype) ~ GAGE,data = subset(data1, age<40), family=binomial)
Age40.GAGE2 <- glm(as.factor(phenotype) ~ GAGE,data = subset(data1, age>=40), family=binomial)
Age40.GAGE3 <- glm(as.factor(phenotype) ~ GAGE+age,data = subset(data1, age<40), family=binomial)
Age40.GAGE4 <- glm(as.factor(phenotype) ~ GAGE+age,data = subset(data1, age>=40), family=binomial)
Age40.GAGE5 <- glm(as.factor(phenotype) ~ GAGE+age+sex,data = subset(data1, age<40), family=binomial)
Age40.GAGE6 <- glm(as.factor(phenotype) ~ GAGE+age+sex,data = subset(data1, age>=40), family=binomial)

adj.table <- modelsummary(list("Age<40"=Age40.GAGE1,"Age 40+"=Age40.GAGE2,
                               "Age<40"=Age40.GAGE3,"Age 40+"=Age40.GAGE4,
                               "Age<40"=Age40.GAGE5,"Age 40+"=Age40.GAGE6), 
                          statistic = 'p.value',
                          stars = c('*' = .1, '**' = .05, '***' = .01),output = "data.frame")

