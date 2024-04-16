install.packages("broom")
install.packages("forestplot") #森林图
install.packages("XLConnect")
install.packages("survminer")
install.packages("survival")
library(XLConnect)
library(broom)
library(survminer)
library(survival) 


data1<-read.table(file=r'{E:\QMY\Result\Data_0322.csv}', 
                  header=TRUE, sep=",")
avg_step_count_y = data1$avg_step_count_y
quantile(avg_step_count_y)

data1_q1 = data1[which(data1$avg_step_count_y <8800.3),] 
data1_q2 = data1[which(8800.3<= data1$avg_step_count_y & data1$avg_step_count_y <10913.8),] 
data1_q3 = data1[which(10913.8<= data1$avg_step_count_y & data1$avg_step_count_y <13164.2),] 
data1_q4 = data1[which(13164.2<= data1$avg_step_count_y),] 

table(data1_q1$termi)
table(data1_q2$termi)
table(data1_q3$termi)
table(data1_q4$termi)

res.Cox <- coxph(Surv(Duration_ukb, termi)~
                   avg_step_count_y+
                   sex+age, data = data1)
summary(res.Cox)


HR<-read.table(file=r'{E:\QMY\Result\HR0319.csv}', 
               header=TRUE, sep=",")


#绘制森林图
hr=sprintf("%.3f",HR$"HR")#获取HR列取小数点后3位
hrLow=sprintf("%.3f",HR$"HR.95L")#获取95%置信区间取HR.95L列小数点后3位
hrHigh=sprintf("%.3f",HR$"HR.95H")#获取95%置信区间HR.95H列小数点后3位
pVal=ifelse(HR$PValue<0.001, "<0.001", sprintf("%.3f", HR$PValue)) #获取P值
Hazard.ratio=paste0(hr,"(",hrLow,"-",hrHigh,")")

split.screen(c(1,2))
screen(1)
xlim = c(0,3)
ylim=c(1,6)
par(mar=c(4,2,1.5,1.5))
print(plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab=""))
text.cex=0.8 
text(0,5:1,HR$X,adj=0,cex=text.cex)
text(1.5-0.5*0.2,5:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,6,'pvalue',
                                                   cex=text.cex,font=2,adj=1)
text(3,5:1,Hazard.ratio,adj=1,cex=text.cex);text(3,6,'Hazard ratio',cex=text.cex,font=2,adj=1,)

screen(2)
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
print(plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio"))#设置x轴的标题
arrows(as.numeric(hrLow),5:1,as.numeric(hrHigh),5:1,angle=90,code=3,length=0.03,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2) #添加中线，设置中线的位置，颜色，类型，宽度
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'blue')#设置中线的取值
points(as.numeric(hr), 5:1, pch = 15, col = boxcolor, cex=1.5)#pch表示点的样式，设置点的大小，颜色
axis(1)
dev.off()

#代谢分组
data_0 = data1[which(data1$metabolic_syd == 0),] 
data_1 = data1[which(data1$metabolic_syd == 1),] 


res.Cox <- coxph(Surv(Duration_ukb, termi)~
                   avg_step_count_y+
                   sex+age+smoke+alcohol+TDI+ethnic_x+season_sin+season_cos, data = data_0)
summary(res.Cox)


HR2<-read.table(file=r'{E:\QMY\Result\HR0322_3.csv}', 
               header=TRUE, sep=",",encoding="UTF8")

hr2=sprintf("%.3f",HR2$"HR")#获取HR列取小数点后3位
hrLow2=sprintf("%.3f",HR2$"HR.95L")#获取95%置信区间取HR.95L列小数点后3位
hrHigh2=sprintf("%.3f",HR2$"HR.95H")#获取95%置信区间HR.95H列小数点后3位
pVal=ifelse(HR2$PValue<0.001, "<0.005", sprintf("%.3f", HR2$PValue)) #获取P值
Hazard.ratio2=paste0(hr2,"(",hrLow2,"-",hrHigh2,")")
n=4
nRow=n+1
ylim=c(1,nRow)
split.screen(c(1,2))
screen(1)
xlim = c(0,4)
par(mar=c(6,3,1,1))
print(plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab=""))
text.cex=0.8
text(0,n:1,HR2$X,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',
                                           cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio2,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',
                                                 cex=text.cex,font=2,adj=1,)

screen(2)
ylim=c(1,nRow)
xlim = c(0.5,1.7)
print(plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio"))
arrows(as.numeric(hrLow2),n:1,as.numeric(hrHigh2),n:1,angle=90,code=3,length=0.03,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr2) > 1, 'red', 'blue')
points(as.numeric(hr2), n:1, pch = 15, col = boxcolor, cex=1.5)
axis(1)


#HR变化曲线
#age
age <- seq(min(data_0$age), max(data_0$age), 0.1) 
dt1 <- data.frame(age = age)
mm1 <- model.matrix(~ -1 + age, data = dt1)
dt1$betaX <- coefficients(res.Cox)[names(coefficients(res.Cox)) == "age"] * (dt1$age-mean(dt1$age))
dt1$betaX <- coefficients(res.Cox)[names(coefficients(res.Cox)) == "age"] * dt1$age
dt1$HR <- exp(dt1$betaX)/100

ggplot(data = dt1, aes(x = age)) + 
  geom_line(aes(y = HR, color = "HR"), linetype = 2) + 
  labs(color = "") + 
  theme(legend.position = c(0.8, 0.8))


#step
avg_step_count_y <- seq(min(data_0$avg_step_count_y),max(data_0$avg_step_count_y), 1)
dt2 <- data.frame(avg_step_count_y=avg_step_count_y)
mm2 <- model.matrix(~ -1 + avg_step_count_y, data = dt2)
dt2$betaX <- coefficients(res.Cox)[names(coefficients(res.Cox)) == "avg_step_count_y"] * (dt2$avg_step_count_y-mean(dt2$avg_step_count_y))

dt2$HR <- exp(dt2$betaX)

ggplot(data = dt2, aes(x = avg_step_count_y)) + 
  geom_line(aes(y = HR, color = "HR"), linetype = 2) + 
  labs(color = "") + 
  theme(legend.position = c(0.8, 0.8))

#alcohol
alcohol <- seq(min(data_0$alcohol),max(data_0$alcohol), 0.1)
dt3 <- data.frame(alcohol=alcohol)
mm3 <- model.matrix(~ -1+alcohol, data = dt3)

dt3$betaX <- coefficients(res.Cox)[names(coefficients(res.Cox)) == "alcohol"] * (dt3$alcohol-mean(dt3$alcohol))
dt3$se <- sqrt(diag(mm3 %*% vcov(res.Cox)[1, 1] %*% t(mm3)))
dt3$HR <- exp(dt3$betaX)

ggplot(data = dt3, aes(x = alcohol)) + 
  geom_line(aes(y = HR, color = "HR"), linetype = 2) + 
  labs(color = "") + 
  theme(legend.position = c(0.8, 0.8))

#TDI
TDI <- seq(min(data_1$TDI),max(data_1$TDI), 0.01)
dt3 <- data.frame(TDI=TDI)
dt3$betaX <- coefficients(res.Cox)[names(coefficients(res.Cox)) == "TDI"] * (dt3$TDI-mean(dt3$TDI))
dt3$HR <- exp(dt3$betaX)

ggplot(data = dt3, aes(x = TDI)) + 
  geom_line(aes(y = HR, color = "HR"), linetype = 2) + 
  labs(color = "") + 
  theme(legend.position = c(0.8, 0.8))











json_walmsley <- load("E:/QMY/data/json_walmsley.RData") 
write.csv(df_json_walmsley, "E:/QMY/data/json_walmsley.csv")
data1<-data.frame(lapply(df_json_walmsley, as.character), stringsAsFactors=FALSE)
write.csv(data1, "E:/QMY/data/json_walmsley.csv")






