
############
# Regression Tree
library(rpart)
library(partykit)

set.seed(1234)

data1$"L-GAGE" <- data1$res
mod_TRG_age <- rpart(phenotype ~ age, data = data1)
mod_TRG_ageGEGE <- rpart(phenotype ~ age+`L-GAGE`, data = data1)
mod_TRG_ageGAGEsex <- rpart(phenotype ~ age+`L-GAGE`+sex, data = data1)

summary(mod_TRG_age)

df <- data.frame(imp = mod_TRG_ageGAGEsex$variable.importance)
df2 <- df %>% 
  tibble::rownames_to_column() %>% 
  dplyr::rename("variable" = rowname) %>% 
  dplyr::arrange(imp) %>%
  dplyr::mutate(variable = forcats::fct_inorder(variable))
ggplot2::ggplot(df2) +
  geom_col(aes(x = variable, y = imp),
           col = "black", show.legend = F) +
  coord_flip() +
  scale_fill_grey() +
  theme_bw()+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 13, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma"),
        axis.text.x=element_text(colour="black", size = 13),
        axis.text.y=element_text(colour="black", size = 13))


#plot(as.party(mod_TRG))
#install.packages("ggparty")
library(ggparty)
ggparty(as.party(mod_TRG_ageGAGEsex)) +
  geom_edge() +
  geom_edge_label(mapping = aes(label = paste(substr(breaks_label, start = 1, stop = 15)))) +
  geom_node_splitvar() +
  # pass list to gglist containing all ggplot components we want to plot for each
  # (default: terminal) node
  geom_node_plot(gglist = list(geom_bar(aes(x = "", fill = phenotype),
                                        position = position_fill()),
                               xlab("phenotype")))





autoplot(as.party(mod_TRG_ageGAGEsex))

plot(as.party(mod_TRG_age),main="Regression tree for predicting MDD based on age")
plot(as.party(mod_TRG_ageGEGE),main="Regression tree for predicting MDD based on age and GAGE")
plot(as.party(mod_TRG_ageGAGEsex),main="Regression tree for predicting MDD based on age, GAGE and sex")
plot(as.party(mod_TRG_ageGAGEsex),cex=2)

rpart.plot::rpart.plot(mod_TRG_age,main="Regression tree for predicting MDD based on age")
rpart.plot::rpart.plot(mod_TRG_ageGEGE,main="Regression tree for predicting MDD based on age and GAGE")
rpart.plot::rpart.plot(mod_TRG_ageGAGEsex,main="Regression tree for predicting MDD based on age, GAGE and sex")

library(ggplot2)
library(tree)
#qplot(Petal.Width, Sepal.Width, data=iris, colour=Species, size=I(4))
#graph + geom_hline(aes(yintercept=2.65)) + geom_vline(aes(xintercept=0.8)) + geom_vline(aes(xintercept=1.75)) + geom_vline(aes(xintercept=1.35))


#######  Continuous threshold regression models ######
#install.packages("chngpt")
library(chngpt)

#types=c("hinge", "segmented", "step", "stegmented") 
data1$MDD <- ifelse(data1$phenotype=="MDD", 1, 0)
fit.mdd=chngptm(formula.1=MDD~1, # antibody.positive+res+sex, 
                formula.2=~age, family="binomial", 
                data=data1,
                #data=subset(data1.antibody, data1.antibody$CMV=="Pos"),
                type="step", 
                #var.type="bootstrap", 
                weights=NULL)
summary(fit.mdd)

fit.mdd$threshold.type
par(mfrow=c(2,1))
plot(fit.mdd)




plot(fit.mdd[["logliks"]])
summary(data1$age)


