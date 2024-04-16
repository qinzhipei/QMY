library(ggplot2)
library(rms) 
library(tuple)
library(XLConnect)
library(broom)
library(survminer)
library(survival) 
library(ggplot2)
library(pheatmap)

data1<-read.table(file=r'{E:\QMY\result\Datamerge0504.csv}', 
                  header=TRUE, sep=",")
dd <- datadist(data1) #数据分布：dd
options(datadist='dd') #为后续程序设定数据环境

fit1<- cph(Surv(Duration_ukb, incident_ASCVD) ~ sex + rcs(F_pseudo_avg,3)+age + smoke + TDI + 
               ethnic_selected+ BMI + WHR
             + metabolic_syd,data=data1)
fit1


HR1<-Predict(fit1,F_pseudo_avg,fun=exp,ref.zero = FALSE)
ggplot(HR1)
