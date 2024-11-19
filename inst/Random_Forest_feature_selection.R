setwd("/Users/liyijie/Documents/Gene age")

library(randomForest)
library(ggplot2)

### Import libraries

#pheno.factor.relevel~gene.expr sample
#data.age.rf.45 <- data.frame(pheno.factor.relevel, t(GxS.covfilter[gene.45,]))
data.age.rf <- data.frame(data.age = phenos.df$age, t(GxS.covfilter))
#Select mtry value with minimum out of bag(OOB) error.
set.seed(1234)
#mtry <- tuneRF(data.age.rf[-1],data.age.rf$data.age, ntreeTry=500,stepFactor=1.5,improve=0.01, trace=TRUE, plot=TRUE)
#best.m <- mtry[mtry[, 2] == min(mtry[, 2]), 1]
#print(mtry)
#print(best.m)

set.seed(1234)
rf <-randomForest(data.age.rf$data.age~.,data=data.age.rf, #mtry=best.m,
                  importance=TRUE,ntree=500)
print(rf)

## Get variable importance, and turn into a data frame
#var_imp <- varImp(rf, scale=FALSE)$importance
var_imp <- varImp(rf)
var_imp <- data.frame(variables=row.names(var_imp), importance=var_imp$Overall)
var_imp_nonzero <- subset(var_imp, var_imp$importance >0)
nrow(var_imp_nonzero) #2812 genes with importance > 0 
hist(var_imp$importance)

#reactome overlap of the top 200 genes
write.table(var_imp_nonzero[order(var_imp_nonzero$importance, decreasing = T),][1:400,]$variables,row.names=F,col.names=F,quote=F)

#generate 200 random genes to feed into reactor, repeated several times to see if some pathways show up by chance. 
write.table(sample(var_imp_nonzero$variables, 200),row.names=F,col.names=F,quote=F)

var_imp_top20 <- var_imp[order(var_imp$importance, decreasing = T),][1:20,]
## Create a plot of variable importance
var_imp_top20 %>%
  
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
  
  ## Add a title
  labs(title='Top 20 genes in Random forest') + 
  
  ## Some layout for the plot
  theme_minimal() + 
  theme(axis.text = element_text(size = 10), 
        axis.title = element_text(size = 15), 
        plot.title = element_text(size = 20), 
  )

var_imp_top20$variables

# predict the biological age by RF

length(var_imp_nonzero$variables)
gene.expression.rf <- data.age.rf %>% select(var_imp_nonzero$variables)


predict.rf <- predict(rf, data.age.rf)

plot(data.age.rf$data.age, predict.rf)
ggplot(data.frame(real.age = data.age.rf$data.age, pred.age = predict.rf,
                  phenotype = data1$phenotype, sex = data1$sex),
       aes(x=real.age, 
           y=pred.age)) +
  geom_point() + 
  #geom_point(aes(col=phenotype,shape=sex),size = 2) + #size = abs(res)
  geom_smooth(method=lm,color = 'blue', fullrange=TRUE) +
  #geom_abline(intercept = 0, slope = 1,color="darkred") + 
  ggtitle(label = ' Random Forest model of age---by phenotype') + 
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab(label = 'Real age') +
  ylab(label = 'Pred age') +
  #geom_point(aes(col=phenotype, shape = sex)) +
  geom_point(aes(col=phenotype,shape=sex),size = 4) + #size = abs(res)
  theme(plot.title = element_text(size=20,face="bold"),
        axis.text=element_text(size=18),
        axis.title=element_text(size=18,face="bold"),
        legend.text = element_text(size=18),
        legend.title = element_text(size=18,face="bold")) +
  guides(colour = guide_legend(override.aes = list(size=5))) + labs(fill = "phenotype")+
  geom_vline(aes(xintercept=20),colour="#BB0000",linetype="dashed")+
  geom_vline(aes(xintercept=30),colour="#BB0000",linetype="dashed")+
  geom_vline(aes(xintercept=40),colour="#BB0000",linetype="dashed")

write.table(sample(var_imp_nonzero$variables, 200),row.names=F,col.names=F,quote=F)
var_imp_top200 <- var_imp[order(var_imp$importance, decreasing = T),][1:200,]

rf_top200 <- data.age.rf %>% select(var_imp_top200$variables)
length(rf_top200)
rf_top200_age <- data.frame(data.age = data1$age, rf_top200)
write.csv(rf_top200_age, "data_rf_top200_age.csv")
