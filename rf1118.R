setwd("C:/Users/wangk/Desktop/Practicum")
load("Allpts5114.RData")
library(gam)
library(randomForest)
library(MASS)
library(caret)
library("e1071")


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

#training and testing
train = sample(1:nrow(CARDIA6), 4000)
set.seed(100)
Y15_rf01 = randomForest(Y15cac~.,data=CARDIA6,subset = train, ntree=100, importance=TRUE,
                       proximity=TRUE,na.action=na.omit)
Y15_rf01#delete the data with the missing value
prediction<-predict(Y15_rf01, CARDIA6[-train, ])
table(prediction,CARDIA6[-train, ]$Y15cac)

set.seed(100)
Y15_rf02 = randomForest(Y15cac~.,data=CARDIA6,subset = train, ntree=100, importance=TRUE,
                        proximity=TRUE,na.action=na.roughfix)
Y15_rf02#Rough Imputation of Missing Values
prediction<-predict(Y15_rf02, CARDIA6[-train, ])
table(prediction,CARDIA6[-train, ]$Y15cac)

##full model to see variable importance
#delete the data with the missing value
set.seed(111)
Y15_rf1 = randomForest(Y15cac~.,data=CARDIA6, ntree=100, importance=TRUE,
                      proximity=TRUE,na.action=na.omit)
importance(Y15_rf1)
Y15_rf1_im<-as.data.frame(importance(Y15_rf1))
plot(Y15_rf1_im$MeanDecreaseAccuracy,Y15_rf1_im$MeanDecreaseGini,
     col="lightblue",pch=19, cex=2, 
     xlab = "MeanDecreaseAccuracy",
     ylab = "MeanDecreaseGini",
     main="Variable Importance for Y15 covariates for Y15 CAC>0
     (delete the data with the missing value)")
text(MeanDecreaseGini~MeanDecreaseAccuracy, labels=rownames(Y15_rf1_im),
     data=Y15_rf1_im, cex=0.8, font=0.5)

#Rough Imputation of Missing Values
set.seed(121)
Y15_rf2 = randomForest(Y15cac~.,data=CARDIA6, ntree=100, importance=TRUE,
                      proximity=TRUE,na.action=na.roughfix)
importance(Y15_rf2)
Y15_rf2_im<-as.data.frame(importance(Y15_rf2))
plot(Y15_rf2_im$MeanDecreaseAccuracy,Y15_rf2_im$MeanDecreaseGini,
     col="lightblue",pch=19, cex=2, 
     xlab = "MeanDecreaseAccuracy",
     ylab = "MeanDecreaseGini",
     main="Variable Importance for Y15 covariates for Y15 CAC>0
     (Rough Imputation of Missing Values)")
text(MeanDecreaseGini~MeanDecreaseAccuracy, labels=rownames(Y15_rf2_im),
     data=Y15_rf2_im, cex=0.8, font=0.5)


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

#training and testing
train = sample(1:nrow(CARDIA), 4000)
set.seed(100)
Y15_rf301 = randomForest(Y15cac~.,data=CARDIA,subset = train, ntree=100, importance=TRUE,
                        proximity=TRUE,na.action=na.omit)
Y15_rf301#delete the data with the missing value
prediction<-predict(Y15_rf301, CARDIA[-train, ])
table(prediction,CARDIA[-train, ]$Y15cac)

set.seed(100)
Y15_rf302 = randomForest(Y15cac~.,data=CARDIA,subset = train, ntree=100, importance=TRUE,
                        proximity=TRUE,na.action=na.roughfix)
Y15_rf302#Rough Imputation of Missing Values
prediction<-predict(Y15_rf302, CARDIA[-train, ])
table(prediction,CARDIA[-train, ]$Y15cac)

##full model to see variable importance
#delete the data with the missing value
set.seed(111)
Y15_rf31 = randomForest(Y15cac~.,data=CARDIA, ntree=100, importance=TRUE,
                       proximity=TRUE,na.action=na.omit)
importance(Y15_rf31)
Y15_rf31_im<-as.data.frame(importance(Y15_rf31))
plot(Y15_rf31_im$MeanDecreaseAccuracy,Y15_rf31_im$MeanDecreaseGini,
     col="lightblue",pch=19, cex=2, 
     xlab = "MeanDecreaseAccuracy",
     ylab = "MeanDecreaseGini",
     main="Variable Importance for Y15 covariates for Y15 CAC>0
     (delete the data with the missing value)")
text(MeanDecreaseGini~MeanDecreaseAccuracy, labels=rownames(Y15_rf31_im),
     data=Y15_rf31_im, cex=0.7, font=0.5)

#Rough Imputation of Missing Values
set.seed(121)
Y15_rf32 = randomForest(Y15cac~.,data=CARDIA, ntree=100, importance=TRUE,
                       proximity=TRUE,na.action=na.roughfix)
importance(Y15_rf32)
Y15_rf32_im<-as.data.frame(importance(Y15_rf32))
plot(Y15_rf32_im$MeanDecreaseAccuracy,Y15_rf32_im$MeanDecreaseGini,
     col="lightblue",pch=19, cex=2, 
     xlab = "MeanDecreaseAccuracy",
     ylab = "MeanDecreaseGini",
     main="Variable Importance for Y15 covariates for Y15 CAC>0
     (Rough Imputation of Missing Values)")
text(MeanDecreaseGini~MeanDecreaseAccuracy, labels=rownames(Y15_rf32_im),
     data=Y15_rf32_im, cex=0.7, font=0.5)


################################################
# GAM plot of prob of Y15 CAC>0 vs Y5 LDL.
################################################
library(mgcv)
GAMdata<-data.frame(Y15cac=ifelse(Allpts$Y15cactot>0,1,0),Y5LDL=Allpts$CL1LDL)
GAMdata<-GAMdata[Allpts$CVDyr.byage>=AGE_new$EX6_AGE&Allpts$DeadFU.byage>=AGE_new$EX6_AGE,]#deleted the patients who had CVD event before exam6,and death before Y15
GAMdata$Y15cac<-as.factor(GAMdata$Y15cac)
gam1<-gam(Y15cac~s(Y5LDL),data=GAMdata,family=binomial(link="logit"))
summary(gam1)
#trans <- function(x)
#  + binomial()$linkinv(x)
plot(gam1,shade = TRUE, trans = plogis)
