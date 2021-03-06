setwd("C:/Users/wangk/Desktop/Practicum")
load("Allpts5114.RData")
#install.packages("rpart")
#install.packages("rpart.plot")
#install.packages("gam")
library(rpart)
library(rpart.plot)
library(gam)

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

#Y15 covariates for Y15 CAC>0
CARDIA6 =data.frame(Male=Allpts$male, Race=Allpts$aarace, Education=Allpts$F03ED, NeverMeetPA=Allpts$PA6,
                    Creatinine=Allpts$FL7CCREAT,CHOL=Allpts$FL1CHOL,PAtotalscore=Allpts$F18TOTAL,
                    EverSmoked=Allpts$F10SMOKE,dislipedima=Allpts$dys_lip6,HTN=Allpts$htn6,BMI=Allpts$F20BMI,
                    GLUCOSE=Allpts$FL7CGLU,SBP=Allpts$F02sbp, DBP=Allpts$F02dbp, EGFR=Allpts$GFR15_CKDEPI,
                    NTRIG= Allpts$FL1NTRIG, HDL= Allpts$FL1HDL, LDL=  Allpts$FL1LDL,DM=  Allpts$diab6  )
CARDIA6$Y15cac=ifelse(Allpts$Y15cactot>0,1,0)#Create binary variables for Y15cactot>0
CARDIA6=CARDIA6[Allpts$CVDyr.byage>=AGE_new$EX6_AGE,]#deleted the patients who had CVD event before exam6
tree1=rpart(Y15cac~Male+Race+Education+NeverMeetPA+EverSmoked+DM+Creatinine+CHOL+PAtotalscore+dislipedima
            +HTN+BMI+GLUCOSE+SBP+DBP+EGFR+HDL+NTRIG+LDL,data=CARDIA6)
rpart.plot(tree1,roundint=FALSE,type=2,extra=1, under = TRUE,cex=1.5)
summary(tree1)
#Y15 covariates for Y15 CAC>0 male and women
tree1m=rpart(Y15cac~Race+Education+NeverMeetPA+EverSmoked+DM+Creatinine+CHOL+PAtotalscore+dislipedima
            +HTN+BMI+GLUCOSE+SBP+DBP+EGFR+HDL+NTRIG+LDL,data=CARDIA6[CARDIA6$Male==1,])
rpart.plot(tree1m,roundint=FALSE,type=2,extra=0, under = TRUE,cex=1)
tree1f=rpart(Y15cac~Race+Education+NeverMeetPA+EverSmoked+DM+Creatinine+CHOL+PAtotalscore+dislipedima
             +HTN+BMI+GLUCOSE+SBP+DBP+EGFR+HDL+NTRIG+LDL,data=CARDIA6[CARDIA6$Male==0,])
rpart.plot(tree1f,roundint=FALSE,type=2,extra=0, under = TRUE,cex=1)

#Y15 covariates for Y20 CAC>0
CARDIA6 =data.frame(Male=Allpts$male, Race=Allpts$aarace, Education=Allpts$F03ED, NeverMeetPA=Allpts$PA6,
                    Creatinine=Allpts$FL7CCREAT,CHOL=Allpts$FL1CHOL,PAtotalscore=Allpts$F18TOTAL,
                    EverSmoked=Allpts$F10SMOKE,dislipedima=Allpts$dys_lip6,HTN=Allpts$htn6,BMI=Allpts$F20BMI,
                    GLUCOSE=Allpts$FL7CGLU,SBP=Allpts$F02sbp, DBP=Allpts$F02dbp, EGFR=Allpts$GFR15_CKDEPI,
                    NTRIG= Allpts$FL1NTRIG, HDL= Allpts$FL1HDL, LDL=  Allpts$FL1LDL,DM=  Allpts$diab6  )
CARDIA6$Y20cac=ifelse(Allpts$Y20cactot>0,1,0)#Create binary variables for Y20cactot>0
CARDIA6=CARDIA6[Allpts$CVDyr.byage>=AGE_new$EX7_AGE&Allpts$DeadFU.byage>=AGE_new$EX7_AGE,]#deleted the patients who had CVD event before exam7,and death before Y20

tree2=rpart(Y20cac~Male+Race+Education+NeverMeetPA+EverSmoked+DM+Creatinine+CHOL+PAtotalscore+dislipedima
            +HTN+BMI+GLUCOSE+SBP+DBP+EGFR+HDL+NTRIG+LDL,data=CARDIA6)
rpart.plot(tree2,roundint=FALSE,type=2,extra=1, under = TRUE,cex=1.5)
summary(tree2)
#Y15 covariates for Y20 CAC>0 male and women
tree2m=rpart(Y20cac~Race+Education+NeverMeetPA+EverSmoked+DM+Creatinine+CHOL+PAtotalscore+dislipedima
             +HTN+BMI+GLUCOSE+SBP+DBP+EGFR+HDL+NTRIG+LDL,data=CARDIA6[CARDIA6$Male==1,])
rpart.plot(tree2m,roundint=FALSE,type=2,extra=0, under = TRUE,cex=1)
tree2f=rpart(Y20cac~Race+Education+NeverMeetPA+EverSmoked+DM+Creatinine+CHOL+PAtotalscore+dislipedima
             +HTN+BMI+GLUCOSE+SBP+DBP+EGFR+HDL+NTRIG+LDL,data=CARDIA6[CARDIA6$Male==0,])
rpart.plot(tree2f,roundint=FALSE,type=2,extra=1, under = TRUE,cex=1)


#Y15 covariates for Y25 CAC>0
CARDIA6 =data.frame(Male=Allpts$male, Race=Allpts$aarace, Education=Allpts$F03ED, NeverMeetPA=Allpts$PA6,
                    Creatinine=Allpts$FL7CCREAT,CHOL=Allpts$FL1CHOL,PAtotalscore=Allpts$F18TOTAL,
                    EverSmoked=Allpts$F10SMOKE,dislipedima=Allpts$dys_lip6,HTN=Allpts$htn6,BMI=Allpts$F20BMI,
                    GLUCOSE=Allpts$FL7CGLU,SBP=Allpts$F02sbp, DBP=Allpts$F02dbp, EGFR=Allpts$GFR15_CKDEPI,
                    NTRIG= Allpts$FL1NTRIG, HDL= Allpts$FL1HDL, LDL=  Allpts$FL1LDL,DM=  Allpts$diab6  )
CARDIA6$Y25cac=ifelse(Allpts$Y25cactot>0,1,0)#Create binary variables for Y25cactot>0
CARDIA6=CARDIA6[Allpts$CVDyr.byage>=AGE_new$EX8_AGE&Allpts$DeadFU.byage>=AGE_new$EX8_AGE,]#deleted the patients who had CVD event before exam8,and death before Y25

tree3=rpart(Y25cac~Male+Race+Education+NeverMeetPA+EverSmoked+DM+Creatinine+CHOL+PAtotalscore+dislipedima
            +HTN+BMI+GLUCOSE+SBP+DBP+EGFR+HDL+NTRIG+LDL,data=CARDIA6)
rpart.plot(tree3,roundint=FALSE,type=2,extra=1, under = TRUE,cex=1.5)
#Y15 covariates for Y25 CAC>0 male and women
tree3m=rpart(Y25cac~Race+Education+NeverMeetPA+EverSmoked+DM+Creatinine+CHOL+PAtotalscore+dislipedima
             +HTN+BMI+GLUCOSE+SBP+DBP+EGFR+HDL+NTRIG+LDL,data=CARDIA6[CARDIA6$Male==1,])
rpart.plot(tree3m,roundint=FALSE,type=2,extra=0, under = TRUE,cex=1)
tree3f=rpart(Y25cac~Race+Education+NeverMeetPA+EverSmoked+DM+Creatinine+CHOL+PAtotalscore+dislipedima
             +HTN+BMI+GLUCOSE+SBP+DBP+EGFR+HDL+NTRIG+LDL,data=CARDIA6[CARDIA6$Male==0,])
rpart.plot(tree3f,roundint=FALSE,type=2,extra=1, under = TRUE,cex=1)




#Y5 covariates for Y15 CAC>0?
CARDIA3 =data.frame(Male=Allpts$male, Race=Allpts$aarace, Education=Allpts$C03ED, NeverMeetPA=Allpts$PA3,
                    CHOL=Allpts$CL1CHOL,PAtotalscore=Allpts$C18TOTAL,
                    EverSmoked=Allpts$C10SMOKE,dislipedima=Allpts$dys_lip3,HTN=Allpts$htn3,BMI=Allpts$C20BMI,
                    SBP=Allpts$c02sbp, DBP=Allpts$c02dbp,
                    NTRIG= Allpts$CL1NTRIG, HDL= Allpts$CL1HDL, LDL=Allpts$CL1LDL,DM=Allpts$diab3 )
CARDIA3$DBP<-as.numeric(as.character(CARDIA3$DBP))
CARDIA3$SBP<-as.numeric(as.character(CARDIA3$SBP))
CARDIA3$Y15cac=ifelse(Allpts$Y15cactot>0,1,0)#Create binary variables for Y15cactot>0
CARDIA3=CARDIA3[Allpts$CVDyr.byage>=AGE_new$EX6_AGE&Allpts$DeadFU.byage>=AGE_new$EX6_AGE,]#deleted the patients who had CVD event before exam6,and death before Y15
tree4=rpart(Y15cac~Male+Race+Education+NeverMeetPA+EverSmoked+DM+CHOL+PAtotalscore+dislipedima
            +HTN+BMI+SBP+DBP+HDL+NTRIG+LDL,data=CARDIA3)
rpart.plot(tree4,roundint=FALSE,type=2,extra=1, under = TRUE,cex=1.5)
#Y5 covariates for Y15 CAC>0 male and women
tree4m=rpart(Y15cac~Male+Race+Education+NeverMeetPA+EverSmoked+DM+CHOL+PAtotalscore+dislipedima
             +HTN+BMI+SBP+DBP+HDL+NTRIG+LDL,data=CARDIA3[CARDIA3$Male==1,])
rpart.plot(tree4m,roundint=FALSE,type=2,extra=1, under = TRUE,cex=1.5)
tree4f=rpart(Y15cac~Male+Race+Education+NeverMeetPA+EverSmoked+DM+CHOL+PAtotalscore+dislipedima
             +HTN+BMI+SBP+DBP+HDL+NTRIG+LDL,data=CARDIA3[CARDIA3$Male==0,])
rpart.plot(tree4f,roundint=FALSE,type=2,extra=1, under = TRUE,cex=1.5)

# GAM plot of prob of Y15 CAC>0 vs Y5 LDL.
GAMdata<-data.frame(Y15cac=ifelse(Allpts$Y15cactot>0,1,0),Y5LDL=Allpts$CL1LDL)
GAMdata<-GAMdata[Allpts$CVDyr.byage>=AGE_new$EX6_AGE&Allpts$DeadFU.byage>=AGE_new$EX6_AGE,]#deleted the patients who had CVD event before exam6,and death before Y15
gam1<-gam(Y15cac~s(Y5LDL),data=GAMdata)
summary(gam1)
plot(gam1,se=TRUE)