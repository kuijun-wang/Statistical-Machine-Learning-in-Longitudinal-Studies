setwd("C:/Users/wangk/Desktop/Practicum")
load("Allpts5114.RData")
library(gam)
library(randomForest)
library(MASS)
library(caret)
library("e1071")
library(ggplot2)
library(tidyverse)
library(caret)
library(glmnet)

#imputation missing value for EXAMAGE
attach(Allpts)
AGE_new<-data.frame(EXAMAGE,EX2_AGE,EX3_AGE,EX4_AGE,EX5_AGE,EX6_AGE,EX7_AGE,EX8_AGE,EX9_AGE)
EX12_AGE<-AGE_new$EXAMAGE+2
AGE_new$EX2_AGE[is.na(AGE_new$EX2_AGE)]=EX12_AGE[is.na(AGE_new$EX2_AGE)]
EX23_AGE<-AGE_new$EX2_AGE+3
AGE_new$EX3_AGE[is.na(AGE_new$EX3_AGE)]=EX23_AGE[is.na(AGE_new$EX3_AGE)]
EX34_AGE<-AGE_new$EX3_AGE+2
AGE_new$EX4_AGE[is.na(AGE_new$EX4_AGE)]=EX34_AGE[is.na(AGE_new$EX4_AGE)]
EX45_AGE<-AGE_new$EX4_AGE+3
AGE_new$EX5_AGE[is.na(AGE_new$EX5_AGE)]=EX45_AGE[is.na(AGE_new$EX5_AGE)]
EX56_AGE<-AGE_new$EX5_AGE+5
AGE_new$EX6_AGE[is.na(AGE_new$EX6_AGE)]=EX56_AGE[is.na(AGE_new$EX6_AGE)]
EX67_AGE<-AGE_new$EX6_AGE+5
AGE_new$EX7_AGE[is.na(AGE_new$EX7_AGE)]=EX67_AGE[is.na(AGE_new$EX7_AGE)]
EX78_AGE<-AGE_new$EX7_AGE+5
AGE_new$EX8_AGE[is.na(AGE_new$EX8_AGE)]=EX78_AGE[is.na(AGE_new$EX8_AGE)]
EX89_AGE<-AGE_new$EX8_AGE+5
AGE_new$EX9_AGE[is.na(AGE_new$EX9_AGE)]=EX89_AGE[is.na(AGE_new$EX9_AGE)]
################################################
#Y15 covariates for Y15 CAC>0
################################################
CARDIA6 =data.frame(Male=Allpts$male, Race=Allpts$aarace, Education=Allpts$F03ED, NeverMeetPA=Allpts$PA6,
                    Creatinine=Allpts$FL7CCREAT,CHOL=Allpts$FL1CHOL,PAtotalscore=Allpts$F18TOTAL,
                    EverSmoked=Allpts$F10SMOKE,dislipedima=Allpts$dys_lip6,HTN=Allpts$htn6,BMI=Allpts$F20BMI,
                    GLUCOSE=Allpts$FL7CGLU,SBP=Allpts$F02sbp, DBP=Allpts$F02dbp, EGFR=Allpts$GFR15_CKDEPI,
                    NTRIG= Allpts$FL1NTRIG, HDL= Allpts$FL1HDL, LDL=  Allpts$FL1LDL,DM=  Allpts$diab6  )
CARDIA6$Y15cac=ifelse(Allpts$Y15cactot>0,1,0)#Create binary variables for Y15cactot>0
CARDIA6=CARDIA6[Allpts$CVDyr.byage>=AGE_new$EX6_AGE,]#deleted the patients who had CVD event before exam6
CARDIA6$Y15cac<-as.factor(CARDIA6$Y15cac)

CARDIA6<-CARDIA6[!is.na(CARDIA6$Y15cac),]
nrow(CARDIA6)#number of the sample has Y15 cac value
nrow(CARDIA6[CARDIA6$Y15cac==0,])#number of the sample with Y15 cac=0
nrow(CARDIA6[CARDIA6$Y15cac==1,])#number of the sample with Y15 cac>0


##full model to see variable importance
#delete the data with the missing value
set.seed(111)
Y15_rf1 = randomForest(Y15cac~.,data=CARDIA6, ntree=100, importance=TRUE,
                      proximity=TRUE,na.action=na.omit)
Y15_rf1_im<-as.data.frame(importance(Y15_rf1))
Y15_rf1_im<-Y15_rf1_im[order(-Y15_rf1_im$MeanDecreaseAccuracy),]
varImpPlot(Y15_rf1)
#Rough Imputation of Missing Values
set.seed(121)
Y15_rf2 = randomForest(Y15cac~.,data=CARDIA6, ntree=100, importance=TRUE,
                      proximity=TRUE,na.action=na.roughfix)
importance(Y15_rf2)
Y15_rf2_im<-as.data.frame(importance(Y15_rf2))
Y15_rf2_im<-Y15_rf2_im[order(-Y15_rf2_im$MeanDecreaseAccuracy),]
varImpPlot(Y15_rf2)


################################################
#Y15 and Y5 covariates for Y15 CAC>0
################################################
#prepare the data
CARDIA<-data.frame(Male=Allpts$male, Race=Allpts$aarace, Education15=Allpts$F03ED, 
            NeverMeetPA15=Allpts$PA6, Creatinine15=Allpts$FL7CCREAT,CHOL15=Allpts$FL1CHOL,
            PAtotalscore15=Allpts$F18TOTAL, EverSmoked15=Allpts$F10SMOKE,
            dislipedima15=Allpts$dys_lip6,HTN15=Allpts$htn6,BMI15=Allpts$F20BMI,
            GLUCOSE15=Allpts$FL7CGLU,SBP15=Allpts$F02sbp, DBP15=Allpts$F02dbp, 
            EGFR15=Allpts$GFR15_CKDEPI,NTRIG15= Allpts$FL1NTRIG, HDL15= Allpts$FL1HDL, 
            LDL15= Allpts$FL1LDL, DM15= Allpts$diab6, Education5=Allpts$C03ED, 
            NeverMeetPA5=Allpts$PA3,CHOL5=Allpts$CL1CHOL,PAtotalscore5=Allpts$C18TOTAL,
            EverSmoked5=Allpts$C10SMOKE,dislipedima5=Allpts$dys_lip3,HTN5=Allpts$htn3,
            BMI5=Allpts$C20BMI, SBP5=Allpts$c02sbp, DBP5=Allpts$c02dbp, 
            NTRIG5= Allpts$CL1NTRIG, HDL5= Allpts$CL1HDL, LDL5=Allpts$CL1LDL,DM5=Allpts$diab3)
CARDIA$Y15cac=ifelse(Allpts$Y15cactot>0,1,0)#Create binary variables for Y15cactot>0
CARDIA=CARDIA[Allpts$CVDyr.byage>=AGE_new$EX6_AGE,]#deleted the patients who had CVD event before exam6
CARDIA$Y15cac<-as.factor(CARDIA$Y15cac)
CARDIA$DBP5<-as.numeric(as.character(CARDIA$DBP5))
CARDIA$SBP5<-as.numeric(as.character(CARDIA$SBP5))

CARDIA<-CARDIA[!is.na(CARDIA$Y15cac),]
nrow(CARDIA)#number of the sample has Y15 cac value
nrow(CARDIA[CARDIA$Y15cac==0,])#number of the sample with Y15 cac=0
nrow(CARDIA[CARDIA$Y15cac==1,])#number of the sample with Y15 cac>0

##full model to see variable importance
#delete the data with the missing value
set.seed(111)
Y15_rf31 = randomForest(Y15cac~.,data=CARDIA, ntree=100, importance=TRUE,
                       proximity=TRUE,na.action=na.omit)
importance(Y15_rf31)
Y15_rf31_im<-as.data.frame(importance(Y15_rf31))
Y15_rf31_im<-Y15_rf31_im[order(-Y15_rf31_im$MeanDecreaseAccuracy),]
varImpPlot(Y15_rf31)

#Rough Imputation of Missing Values
set.seed(121)
Y15_rf32 = randomForest(Y15cac~.,data=CARDIA, ntree=100, importance=TRUE,
                       proximity=TRUE,na.action=na.roughfix)
importance(Y15_rf32)
Y15_rf32_im<-as.data.frame(importance(Y15_rf32))
Y15_rf32_im<-Y15_rf32_im[order(-Y15_rf32_im$MeanDecreaseAccuracy),]
varImpPlot(Y15_rf32)


#########################################
############logistic regression##########

############prepare the data##########
#Y15 covariates for Y15 CAC>0 
#First ten
Rdata15<-data.frame(Male=CARDIA6$Male,CHOL=CARDIA6$CHOL,LDL=CARDIA6$LDL,HDL=CARDIA6$HDL,
                    Creatinine=CARDIA6$Creatinine, DBP=CARDIA6$DBP,NTRIG=CARDIA6$NTRIG,
                    SBP=CARDIA6$SBP, BMI=CARDIA6$BMI,dislipedima=CARDIA6$dislipedima)
Rdata15$Y15cac<-as.factor(CARDIA6$Y15cac)
###logistic regression##
#Apply logistic regression to Y15cac > 0, and show the covariate effects;
Y15_logit <- glm(Y15cac ~ ., data=Rdata15, family=binomial(link="logit"))
summary(Y15_logit)



###LASSO###
x = model.matrix(Y15cac ~ . -1, data = Rdata15)
y = Rdata15$Y15cac[-nas]
nas=c()
for (i in 1:3017){if (sum(is.na(Rdata15[i,]))>0){nas=c(nas,i)}}
fit.lasso = glmnet(x, y)
plot(fit.lasso, xvar = "lambda", label = TRUE)
cv.lasso = cv.glmnet(x, y)
plot(cv.lasso)
coef(cv.lasso)
###stepwise###



############prepare the data##########
#Y15 and Y5 covariates for Y15 CAC>0 
#First ten
Rdata515<-data.frame(Male=CARDIA$Male,CHOL15=CARDIA$CHOL15,LDL15=CARDIA$LDL15,
                     CHOL5=CARDIA$CHOL5,SBP5=CARDIA$SBP5,LDL5=CARDIA$LDL5,
                     HDL15=CARDIA$HDL15,HDL5=CARDIA$HDL5,BMI15=CARDIA$BMI15,
                     DBP5=CARDIA$DBP5)
Rdata515$Y15cac<-as.factor(CARDIA$Y15cac)
###logistic regression##
#Apply logistic regression to Y15cac > 0, and show the covariate effects;
Y515_logit <- glm(Y15cac ~ ., data=Rdata515, family=binomial(link="logit"))
summary(Y515_logit)
