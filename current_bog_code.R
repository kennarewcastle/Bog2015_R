### Currrent code for bog paper analyses and figures
### Kenna Rewcastle
### March 6th, 2018
### Revised June 14, 2019

library(ggplot2)
library(dplyr)
library(MASS)
library(stats)
library(gridExtra)
library(knitr)
library(kableExtra)


# Master Data Frame ******SKIP THIS SECTION ***** ------------------------------------------
#bog1<-read.csv("bog_2015.csv")
#bog<-filter(bog1,GWC!="NA")

# Compress enzymes into a total nutrient enzyme activity parameter and a total carbon enzyme activity parameter.

# nutrient enzymes
#NAG<-bog$NAG
#PHOS<-bog$PHOS
#LAP<-bog$LAP 
#enzyNut<-NAG+PHOS+LAP

# carbon enzymes
#CBH<-bog$CBH
#AG<-bog$AG
#BG<-bog$BG
#enzyC<-CBH+AG+BG

# add total enzyme columns to main dataframe
#bog<-cbind(bog,enzyNut,enzyC)

#bog[51,31]<-NA # Excludes one very clear outlier (screwed up analysis?) MBC point

# data frame without NAs for step AIC function
#realbog<-filter(bog,MBC.mg.C.g.1.dry.soil!="NA")

#write.csv(bog,file="Master_BOG2015_forR.csv",row.names=FALSE)
#write.csv(realbog,file="BOG2015_for_stepAIC.csv",row.names=FALSE)


# ANALYSES FOR MANUSCRIPT BODY --------------------------------------------
bog<-read.csv("Master_BOG2015_forR.csv")
realbog<-read.csv("BOG2015_for_stepAIC.csv")
bog_starch<-filter(realbog,carbon=="label")

# Scale C and N values  ----------------

# Reviewer 2 didn't like the idea of adding enzyme activities together to come up with a composite enzyme activity. To address this, I've scaled the enzyme activities so that the magnitude each of the three individual activities in the composite variable are comparable. Resp x each enzyme and d13 by each enzyme will be included in the supplemental material.


# RESPIRATION --------------------------------------------
resp<-bog$D5.CO2

# Carbon enzymes
bg_resp<-lm(resp~bog$BG)
summary(bg_resp) # p = 0.003, R2 = 0.1327
qplot(x=BG,y=D5.CO2,data=bog)

ag_resp<-lm(resp~bog$AG)
summary(ag_resp) # p = 0.990, R2 = -0.01754
qplot(x=BG,y=D5.CO2,data=bog)

cbh_resp<-lm(resp~bog$CBH)
summary(cbh_resp) # p = 0.007 R2 = 0.1067
qplot(x=BG,y=D5.CO2,data=bog)


# Nutrient enzymes
nag_resp<-lm(resp~bog$NAG) 
summary(nag_resp) # p = 0.063 R2 = 0.042
qplot(bog$NAG,resp)

phos_resp<-lm(resp~bog$PHOS)
summary(phos_resp) # p = 0.0001 R2 = 0.2157
qplot(bog$PHOS,resp)

lap_resp<-lm(resp~bog$LAP)
summary(lap_resp) # p = 0.0007938 R2 = 0.1662
qplot(bog$LAP,resp)

##### Model selection using step AIC for resp predictors. These parameters are for the scenario where we include 

resp<-realbog$D5.CO2
MBC<-realbog$MBC.mg.C.g.1.dry.soil
# logMBC<-log(MBC) # Data is slightly more normally distributed... Not enough to matter!
# resp_enzyC<-realbog$enzyC
# resp_enzyNut<-realbog$enzyNut
sat<-realbog$per..saturation
pH<-realbog$pH
DOC<-realbog$DOC..mg.C.g.1.dry.soil.
CN<-realbog$peat.C.N
spruce<-realbog$Dist..Spruce..m.
blub<-realbog$Dist..Blueberry..m.

enzy_C_scale<-data.frame("AG_scale"=realbog$AG,"BG_scale"=realbog$BG,"CBH_scale"=realbog$CBH)
enzy_C_scale<-scale(enzy_C_scale,center=FALSE,scale=TRUE)
sum_C_scale<-enzy_C_scale[,1] + enzy_C_scale[,2] + enzy_C_scale[,3]

enzy_nut_scale<-data.frame("NAG_scale"=realbog$NAG,"PHOS_scale"=realbog$PHOS,"LAP_scale"=realbog$LAP)
enzy_nut_scale<-scale(enzy_nut_scale,center=FALSE,scale=TRUE)
sum_nut_scale<-enzy_nut_scale[,1] + enzy_nut_scale[,2] + enzy_nut_scale[,3]

AG<-realbog$AG
BG<-realbog$BG
CBH<-realbog$CBH
NAG<-realbog$NAG
LAP<-realbog$LAP
PHOS<-realbog$PHOS

# These parameters are on the full dataset (n=59) and excludes MBC and DOC
# resp<-bog$D5.CO2
# resp_enzyC<-bog$enzyC
# resp_enzyNut<-bog$enzyNut
# sat<-bog$per..saturation
# pH<-bog$pH
# CN<-bog$peat.C.N
# spruce<-bog$Dist..Spruce..m.
# blub<-bog$Dist..Blueberry..m.

# Figure out parameter order in initial model by comparing resp~paramater R2

# predict1<-data.frame(MBC,resp_enzyC,resp_enzyNut,sat,pH,DOC,CN,spruce,blub)
# predict1<-data.frame(sum_C_scale,sum_nut_scale,sat,pH,CN,spruce,blub,MBC,DOC)
predict1<-data.frame(AG,BG,CBH,NAG,PHOS,LAP,sat,pH,CN,spruce,blub,MBC,DOC)

param_R2<-function(predict,response){
  N<-ncol(predict)
  R2vect<-c()
  for (i in 1:N) {
    regression<-lm(response~predict[,i])
    R2vect[i]<-summary(regression)$adj.r.squared
  }
  results<-data.frame(names(predict),R2vect)
  colnames(results)<-c("Predictor","R2")
  results<-results[order(-results$R2),] 
  return(results)
}

resp_rs<-param_R2(predict=predict1,response=resp)
print(resp_rs)

# StepAIC model selection 

library(MASS)

mod001<-lm(resp~sat+ PHOS+ LAP+ MBC+ BG+ spruce+ CBH+ NAG+ DOC+ pH+ blub+ CN+ AG)
# mod001<-lm(resp~sat+ resp_enzyNut+ spruce+ resp_enzyC+ pH+ CN+ blub)
startmod<-lm(resp~1)
stepmod001<-stepAIC(startmod,direction="both",scope=resp~sat+ PHOS+ LAP+ MBC+ BG+ spruce+ CBH+ NAG+ DOC+ pH+ blub+ CN+ AG) # best model: resp ~ sat + BG + spruce + DOC, AIC = 19.5
summary(stepmod001)

mod1<-lm(resp~sat)
summary(mod1)
mod2<-lm(resp~sat+BG)
summary(mod2)
mod3<-lm(resp~sat+BG+spruce)
summary(mod3)
mod4<-lm(resp~sat+BG+spruce+DOC)
summary(mod4)

##### Linear regressions between parameters in best fit model and resp

resp_sat<-lm(resp~sat)
summary(resp_sat) # p < 0.001, R2 = 0.411

resp_sat_plot<-ggplot(data=realbog,aes(x=sat,y=resp)) +
  geom_smooth(method=lm,formula=y~x,colour="black",size=1.25) +
  geom_point() +
  ylim(0,8.1) +
  ylab(expression(bold(paste("Total Respiration (",mu,"mol"," ","CO"[2]," ","m"^-2," s"^-1,")")))) +
  xlab(label="Water Saturation (%)") +
  coord_cartesian(xlim=c(40,130), ylim=c(0,8.05)) + # Makes sure that CI and line fit extend to edge of graph
  annotate("text", x = 40, y = 7.9, label = "A", size=8, color="black") +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.text=element_text(colour="black",size=10),
        axis.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA)
  )

resp_spruce<-lm(resp~spruce)
summary(resp_spruce) # p = 0.007, R2 = 0.13

resp_spruce_plot<-ggplot(data=realbog,aes(x=spruce,y=resp)) +
  geom_smooth(method=lm,formula=y~x,fullrange=TRUE,colour="black",size=1.25) +
  geom_point() +
  ylab(label=NULL) +
  xlab(label="Distance from Nearest Spruce Tree (m)") +
  ylim(0,8.1) +
  xlim(0,2.45) +
  annotate("text", x = 0, y = 7.9, label = "B", size=8, color="black") +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.text=element_text(colour="black",size=10),
        axis.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA)
  )

resp_BG_mod<-lm(resp~BG) # p = 0.007, R2 = 0.129
summary(resp_BG_mod)

resp_BG_plot<-ggplot(data=realbog,aes(x=BG,y=resp)) +
  geom_smooth(method=lm,formula=y~x,colour="black",size=1.25) +
  geom_point() +
  ylab(expression(bold(paste("Total Respiration (",mu,"mol"," ","CO"[2]," ","m"^-2," s"^-1,")")))) +
  xlab(expression(bold(paste("Potential BG Activity"," (nmol"," g"^-1," h"^-1,")")))) +
  ylim(0,8.1) +
  annotate("text", x = 1750, y = 7.9, label = "C", size=8, color="black") +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.text=element_text(colour="black",size=10),
        axis.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA)
  )

resp_DOC<-lm(resp~DOC) # p = 0.16, R2 = 0.022
summary(resp_DOC)

resp_DOC_plot<-ggplot(data=realbog,aes(x=DOC,y=resp)) +
  geom_smooth(method=lm,formula=y~x,colour="black",size=1.25) +
  geom_point() +
  ylab(label=NULL) +
  xlab(expression(bold(paste("DOC (mg C ","g"^-1," dry soil)")))) +
  ylim(0,8) +
  annotate("text", x = 1.5, y = 7.9, label = "D", size=8, color="black") +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.text=element_text(colour="black",size=10),
        axis.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA)
  )

##### Panneled figure for best predictor regressions

library(gridExtra)
Figure2<-grid.arrange(resp_sat_plot,resp_spruce_plot,resp_BG_plot,resp_DOC_plot,nrow=2)
ggsave(filename="Figure2.pdf",plot=Figure2,dpi=300,width=8.75,height=7.5,units="in")

##### Model selection using step AIC for d13 predictors.

bog_starch<-filter(realbog,carbon=="label")

d13<-bog_starch$D5.d13
MBC<-bog_starch$MBC.mg.C.g.1.dry.soil
starch_enzyC<-bog_starch$enzyC
starch_enzyNut<-bog_starch$enzyNut
sat<-bog_starch$per..saturation
pH<-bog_starch$pH
DOC<-bog_starch$DOC..mg.C.g.1.dry.soil.
CN<-bog_starch$peat.C.N
spruce<-bog_starch$Dist..Spruce..m.
blub<-bog_starch$Dist..Blueberry..m.
resp<-bog_starch$D5.CO2
AG<-bog_starch$AG
BG<-bog_starch$BG
CBH<-bog_starch$CBH
NAG<-bog_starch$NAG
PHOS<-bog_starch$PHOS
LAP<-bog_starch$LAP

# Scale enzyme parameters
enzy_C_scale<-data.frame("AG_scale"=bog_starch$AG,"BG_scale"=bog_starch$BG,"CBH_scale"=bog_starch$CBH)
enzy_C_scale<-scale(enzy_C_scale,center=FALSE,scale=TRUE)
sum_C_scale<-enzy_C_scale[,1] + enzy_C_scale[,2] + enzy_C_scale[,3]

enzy_nut_scale<-data.frame("NAG_scale"=bog_starch$NAG,"PHOS_scale"=bog_starch$PHOS,"LAP_scale"=bog_starch$LAP)
enzy_nut_scale<-scale(enzy_nut_scale,center=FALSE,scale=TRUE)
sum_nut_scale<-enzy_nut_scale[,1] + enzy_nut_scale[,2] + enzy_nut_scale[,3]

# Figure out parameter order in initial model by comparing resp~paramater R2

predict2<-data.frame(MBC,AG,BG,CBH,NAG,LAP,PHOS,sat,pH,DOC,CN,spruce,blub)
# predict2<-data.frame(starch_enzyC,starch_enzyNut,sat,pH,CN,spruce,blub)

starch_rs<-param_R2(predict=predict2,response=d13)
print(starch_rs)

# StepAIC model selection
library(MASS)

mod002<-lm(d13~spruce+ MBC+ BG+ CBH+ sat+ LAP+ PHOS+ AG+ NAG+ DOC+ CN+ blub+ pH)
# mod002<-lm(d13~spruce+ starch_enzyC+ sat+ starch_enzyNut+ CN+ pH+ blub)
# Best model: resp~ spruce+ starch_enzyC + starch_enzyNut 
startmod<-lm(d13~1)
stepmod002<-stepAIC(startmod,direction="both",scope=d13~spruce+ MBC+ BG+ CBH+ sat+ LAP+ PHOS+ AG+ NAG+ DOC+ CN+ blub+ pH)
summary(stepmod002) # best model = spruce + BG + PHOS

mod01<-lm(d13~spruce)
summary(mod01)
mod02<-lm(d13~spruce+ BG)
summary(mod02)
mod03<-lm(d13~spruce + BG + PHOS)
summary(mod03)

##### Linear regressions between parameters in best fit model and d13

d13_spruce<-lm(d13~spruce)
summary(d13_spruce) # p = 0.002, R2 = 0.35

d13_spruce_plot<-ggplot(data=bog_starch,aes(x=spruce,y=d13)) +
  geom_smooth(method=lm,formula=y~x,fullrange=TRUE,colour="black",size=1.25) +
  geom_point() +
  ylab(label=expression(bold(paste("\u03B4"^{bold("13")}, "C Respired (\u2030)")))) +
  xlab(label="Distance from Nearest Spruce Tree (m)") +
  ylim(0,400) +
  xlim(0,2.45) +
  annotate("text", x = 0, y = 395, label = "A", size=8, color="black") +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.text=element_text(colour="black",size=10),
        axis.title.y=element_text(size=16,face="bold"),
        axis.title.x=element_text(size=16,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA)
  )

d13_BG<-lm(d13~BG) # p = 0.010, R2 = 0.2515
summary(d13_BG)

d13_BG_plot<-ggplot(data=bog_starch,aes(x=BG,y=d13)) +
  geom_smooth(method=lm,formula=y~x,colour="black",size=1.25) +
  geom_point() +
  ylab(label=NULL) +
  xlab(expression(bold(paste("Potential BG Activity"," (nmol"," g"^-1," h"^-1,")")))) +
  ylim(0,400) +
  annotate("text", x = 1750, y = 395, label = "B", size=8, color="black") +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.text=element_text(colour="black",size=10),
        axis.title.y=element_text(size=16,face="bold"),
        axis.title.x=element_text(size=16,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA)
  )

d13_PHOS<-lm(d13~PHOS) # p = 0.048, R2 = 0.14
summary(d13_PHOS)

d13_PHOS_plot<-ggplot(data=bog_starch,aes(x=PHOS,y=d13)) +
  geom_smooth(method=lm,formula=y~x,colour="black",size=1.25) +
  geom_point() +
  ylab(label=NULL) +
  xlab(expression(bold(paste("Potential PHOS Activity"," (nmol"," g"^-1," h"^-1,")")))) +
  ylim(0,400) +
  annotate("text", x = 1210, y = 395, label = "C", size=8, color="black") +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.text=element_text(colour="black",size=10),
        axis.title.y=element_text(size=16,face="bold"),
        axis.title.x=element_text(size=16,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA)
  )

##### Panneled figure for best predictor regressions

library(gridExtra)
Figure1<-grid.arrange(d13_spruce_plot,d13_BG_plot,d13_PHOS_plot,nrow=1)
# ggsave(filename="Figure1.jpeg",plot=Figure1,dpi=300,width=11,units="in")

###### Does respiration predict d13?
starch_resp<-lm(d13~resp)
summary(starch_resp) # no, p=0.88
qplot(x=resp,y=d13)

##### Mean bog properties table

# Use full dataset that includes NAs, bog1
bog1<-read.csv("bog_2015.csv")
names(bog1)

blub<-bog1$Dist..Blueberry..m.
spruce<-bog1$Dist..Spruce..m.
sat<-bog1$per..saturation
pH<-bog1$pH
peatCN<-bog1$peat.C.N
DOC<-bog1$DOC..mg.C.g.1.dry.soil.
DON<-bog1$DON..mg.N.g.1.dry.soil.
MBC<-bog1$MBC.mg.C.g.1.dry.soil
MBN<-bog1$MBN.mg.N.g.1.dry.soil
enzyC<-(bog1$AG + bog1$BG + bog1$CBH)
enzyN<-(bog1$NAG + bog1$PHOS + bog1$LAP)

tableDat<-data.frame(blub,spruce,sat,pH,peatCN,DOC,DON,MBC,MBN,enzyC,enzyN)

#########################################################################################################
# FUNCTION: bogMean
# Calculates site-level mean and standard deviation of all bog properties
# input: dat: data frame where columns are individual parameters and rows are observations across parameters for an individual mesocosm
# output: data frame listing mean and sd for each parameter
#--------------------------------------------------------------------------------------------------------

bogMean<-function(dat=NULL) {
  if(is.null(dat)) { # Minimalist code for default data frame to test function
    var1<-runif(10)
    var2<-runif(10)
    dat<-data.frame(var1,var2)
  }
  
  N<-ncol(dat)
  
  param<-c()
  meanDat<-c()
  sdDat<-c()
  sumDat<-data.frame(param,meanDat,sdDat)
  
  namesDat<-names(dat)
  
  for(i in 1:N) {
    meanVar<-mean(dat[,i], na.rm=TRUE) # removes NAs from a column before calculating mean
    sdVar<-sd(dat[,i], na.rm=TRUE)
    name<-namesDat[i]
    
    sumDat[i,1]<-name
    sumDat[i,2]<-meanVar
    sumDat[i,3]<-sdVar
  }
  
  names(sumDat)<-c("Parameter","Mean","Standard Deviation")
  return(sumDat)
}
#_________________________________________________________________

tableOutput<-bogMean(dat=tableDat)
write.csv(tableOutput,file="table_means.csv") # Creates a .csv file in working directory with parameters and their means and sds. This data goes into a .csv file that will be used as the skeleton for a kable table.

# See markdown file --> knit to pdf for using kable to make table look pretty for manuscript.



# SUPPLEMENTARY FIGURES ---------------------------------------------------

##### Mesh size does not affect microbial biomass (plot and stats)

mesh<-as.factor(bog$treatment)
MBC<-bog$MBC.mg.C.g.1.dry.soil

mesh_MBC<-lm(MBC~mesh)
anova(mesh_MBC) # p = 0.204, F = 1.646
summary(mesh_MBC)

mesh_cols<-c("1"="grey60","2"="grey32","3"="black")
mesh_scale<-c(1,2,3)

mesh_MBC_fig<-ggplot(data=bog,aes(x=mesh,y=MBC)) +
  geom_boxplot(lwd=1.5) +
  labs(x="Mesh Size",y=bquote(bold('Microbial Biomass C (mg C '*g^-1*' dry soil)'))) +
  # scale_colour_manual(values=mesh_cols,labels=c("1.45 mm","55 \u03BCm","5 \u03BCm"),guide=FALSE) + 
  # weird text is unicode for mu
  annotate("text",x=3.4,y=4.1,label="A",size=10) +
  scale_x_discrete(limits=mesh_scale, labels=c("1.45 mm","55 \u03BCm","5 \u03BCm")) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text=element_text(colour="black",size=10),
        axis.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA))

mesh_MBC_fig

##### Mesh size does not affect peat water saturation (plot and stats)

mesh<-as.factor(bog$treatment)
sat<-bog$per..saturation

mesh_sat<-lm(sat~mesh)
anova(mesh_sat) # p = 0.09372, F = 2.470
summary(mesh_sat)

mesh_cols<-c("1"="grey60","2"="grey32","3"="black")
mesh_scale<-c(1,2,3)

mesh_sat_fig<-ggplot(data=bog,aes(x=mesh,y=sat)) +
  geom_boxplot(lwd=1.5) +
  labs(x="Mesh Size",y="% Water Saturation") +
  # scale_colour_manual(values=mesh_cols,labels=c("1.45 mm","55 \u03BCm","5 \u03BCm"),guide=FALSE) + 
  # weird text is unicode for mu
  annotate("text",x=3.4,y=128,label="B",size=10) +
  scale_x_discrete(limits=mesh_scale, labels=c("1.45 mm","55 \u03BCm","5 \u03BCm")) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text=element_text(colour="black",size=10),
        axis.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA))

mesh_sat_fig

##### Paneled figure for supplemental Figure S1

library(gridExtra)
S1<-grid.arrange(mesh_MBC_fig,mesh_sat_fig,nrow=1)
ggsave(filename="S1.jpeg",plot=S1,dpi=300)

##### Labeled starch treatment was successful (increase in 13C CO2 respiration in labeled mesocosms)

carbon<-as.factor(bog$carbon)
d13<-bog$D5.d13

d13effect<-t.test(formula=d13~carbon) # p < 0.001, t = -9.5717
labeled<-filter(bog,carbon=="label")
water<-filter(bog,carbon=="control")
mean(labeled$D5.d13) # 134.38
sd(labeled$D5.d13) # 82.99
mean(water$D5.d13) # -7.82
sd(water$D5.d13) # 6.65

# Figure for S2
d13_scale<-c(1,2)

d13effect_fig<-ggplot(data=bog,aes(x=carbon,y=d13)) +
  geom_boxplot(lwd=1.5) +
  labs(x=expression(bold(paste(text=" "^{bold("13")},"C"," Starch Addition"))),y=expression(bold(paste("\u03B4"^{bold("13")}, "C Respired (\u2030)")))) +
  scale_x_discrete(labels=c("control"="Control","label"="Labeled")) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA))

## Save plot
ggsave(filename="S2.jpeg",plot=d13effect_fig,dpi=300)


#### Starch addition treatment doesn't affect: total resp, microbial biomass, enzyme activity (all 6), peat saturation, pH

carbon<-as.factor(bog$carbon)
resp<-bog$D5.CO2
MBC<-bog$MBC.mg.C.g.1.dry.soil
AG<-bog$AG
BG<-bog$BG
CBH<-bog$CBH
NAG<-bog$NAG
PHOS<-bog$PHOS
LAP<-bog$LAP
sat<-bog$per..saturation
pH<-bog$pH
DOC<-bog$DOC..mg.C.g.1.dry.soil.
CN<-bog$peat.C.N


# t-tests first
t.test(formula=resp~carbon) # t = -0.63, p = 0.53
t.test(formula=MBC~carbon) # t = -0.98, p = 0.33
t.test(formula=AG~carbon) # t = 0.43, p = 0.67
t.test(formula=BG~carbon) # check output summary below for t and p values
t.test(formula=CBH~carbon) 
t.test(formula=NAG~carbon) 
t.test(formula=PHOS~carbon) 
t.test(formula=LAP~carbon) 
starch_enzyC<-t.test(formula=enzyC~carbon)
starch_enzyNut<-t.test(formula=enzyNut~carbon)
starch_sat<-t.test(formula=sat~carbon)  
starch_pH<-t.test(formula=pH~carbon) 
starch_DOC<-t.test(formula=DOC~carbon) 
starch_CN<-t.test(formula=CN~carbon) 

# MIGHT need to build a supplemental figure displaying all of this

##### Respiration and d13 by DOC and MBC to validate removing from model selection

DOC_resp<-lm(resp~DOC)
summary(DOC_resp)

MBC_resp<-lm(resp~MBC)
summary(MBC_resp)

DOC_d13<-lm(d13~DOC)
summary(DOC_d13)

MBC_d13<-lm(d13~MBC)
summary(MBC_d13)

##### Atom percent 13CO2 across sampling days

five_day_atom<-c(bog_starch$D1.13C.atom..,bog_starch$D2.13C.atom..,bog_starch$D3.13C.atom..,bog_starch$D4.13C.atom..,bog_starch$D5.13C.atom..)

atom_days<-rep(1:5,each=22)

atom_dat<-data.frame("Day"=atom_days,"Atom_Percent"=five_day_atom)
atom_dat$Day<-as.factor(atom_dat$Day)

# I think boxplots are the correct way to do this...
atom_days_fig<-ggplot(data=atom_dat,aes(x=Day,y=Atom_Percent)) +
    geom_boxplot(lwd=1) +
    labs(y=expression(bold(paste(text=" "^{bold("13")},"C"," Atom % in ","CO"[bold(2)]))),x=expression(bold(paste(text="Days After"," "^{bold("13")},"C","-Starch Addition")))) +
    theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
          axis.text.y=element_text(colour="black",size=10),
          axis.text.x=element_text(colour="black",size=14),
          axis.title=element_text(size=14,face="bold"),
          panel.border=element_rect(fill=NA,colour="black",size=1.5),
          panel.background=element_rect(fill=NA))

atom_days_smooth<-ggplot(data=atom_dat,aes(x=as.numeric(atom_dat$Day),y=Atom_Percent)) +
  geom_point() +
  geom_smooth(col="black") +
  labs(y=expression(bold(paste(text=" "^{bold("13")},"C"," Atom % in ","CO"[bold(2)]))),x=expression(bold(paste(text="Days After"," "^{bold("13")},"C","-Starch Addition")))) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA))

##### Atom percent 13CO2 across sampling days

five_day_d13<-c(bog_starch$D1.d13,bog_starch$D2.d13,bog_starch$D3.d13,bog_starch$D4.d13,bog_starch$D5.d13)

d13_days<-rep(1:5,each=22)

d13_dat<-data.frame("Day"=d13_days,"d13"=five_day_d13)
d13_dat$Day<-as.factor(d13_dat$Day)

d13_days_fig<-ggplot(data=d13_dat,aes(x=Day,y=d13)) +
  geom_boxplot(lwd=1) +ylab(label=expression(bold(paste("\u03B4"^{bold("13")}, "C Respired (\u2030)")))) +
  labs(y=expression(bold(paste("\u03B4"^{bold("13")}, "C Respired (\u2030)"))),x=expression(bold(paste(text="Days After"," "^{bold("13")},"C","-Starch Addition")))) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA))

ggsave(filename="d13_5_days.jpg",plot=d13_days_fig)

# Distrubtion in mesocosm proximity to spruce, blueberry ------------------

Distance<-c(bog$Dist..Spruce..m.,bog$Dist..Blueberry..m.)
Tree<-c(rep("spruce",times=59),rep("blub",times=59))
dist_df<-data.frame("Tree"=Tree,"Distance"=Distance)

ggplot(data=dist_df, aes(dist_df$Distance,fill=dist_df$Tree)) + 
  geom_density(alpha=0.7) +
  labs(fill="Focal plant species",x="Mesocosm distance to nearest stem (m)",y="Density") +
  scale_fill_discrete(labels=c("spruce"="black spruce","blub"="blueberry")) +
  theme_classic()

ggsave(filename="Plant_dist.jpg")



# Is water-microbial activity relationship driven by pH? ------------------

#### Water saturation vs. pH
sat_pH<-ggplot(data=bog,aes(x=bog$per..saturation,y=bog$pH)) +
  geom_smooth(method=lm,formula=y~x,colour="black",size=1.25) +
  geom_point() +
  ylab(expression(bold(paste("pH")))) +
  xlab(label="Water Saturation (%)") +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.text=element_text(colour="black",size=10),
        axis.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA)
  )

sat_pH_mod<-lm(bog$pH~bog$per..saturation)
summary(sat_pH_mod) # Significant, p =0.02575, F = 5.244

resp_sat_pH<-lm(bog$D5.CO2~bog$per..saturation*bog$pH)
summary(resp_sat_pH)

resp_pH<-lm(bog$D5.CO2~bog$pH)
summary(resp_pH)
# While there is a negative relationship between peat saturation and pH (p =0.02575, F = 5.244), this relationship does not seem to confound the strong correlation between peat saturation and respiration, as there is no indication of an impact of pH alone (p = 0.08687, F = 3.035) or an interactive effect between pH and peat saturation on heterotrophic respiration (p = 0.983, F = 12.87).