install.packages("broom")
install.packages("forestplot") #森林图
install.packages("XLConnect")
install.packages("survminer")
install.packages("survival")
install.packages("tuple")
install.packages("patchwork")
library("tuple")
library(XLConnect)
library(broom)
library(survminer)
library(survival) 
library(ggplot2)


data1<-read.table(file=r'{E:\QMY\Result\Data_0322.csv}', 
                  header=TRUE, sep=",")
avg_step_count_y = data1$avg_step_count_y
quantile(avg_step_count_y)



#对整个人群作cox回归，计算出Hazard Ratio
res.Cox <- coxph(Surv(Duration_ukb, termi)~
                   avg_step_count_y+
                   sex+age, data = data1)
summary(res.Cox)

#将每个变量的Hazard Ratio维护到一个表里，并导入
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







#根据代谢病分组
data_0 = data1[which(data1$metabolic_syd == 0),]  #无代谢病人群
data_1 = data1[which(data1$metabolic_syd == 1),]  #有代谢病人群

#无代谢病人群作Cox回归
res.Cox <- coxph(Surv(Duration_ukb, termi)~
                   avg_step_count_y+
                   sex+age+smoke+alcohol+TDI+ethnic_x+season_sin+season_cos, data = data_0)
summary(res.Cox)

#绘制森林图
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
data1<-read.table(file=r'{E:\QMY\Result\Data0411.csv}', 
                  header=TRUE, sep=",")
#代谢病分组
data_0 = data1[which(data1$metabolic_syd == 0),]  #无代谢病人群
data_1 = data1[which(data1$metabolic_syd == 1),]  #有代谢病人群


#吸烟分组
data_0 = data1[which(data1$smoke == 0),]  #无吸烟史
data_1 = data1[which(data1$smoke == 1),]  #曾吸烟
data_2 = data1[which(data1$smoke == 2),]  #在吸烟

#年龄
data_0 = data1[which(data1$age <= 50),]  
data_1 = data1[which(50 < data1$age & data1$age <= 60),]  
data_2 = data1[which(data1$age >= 60),] 

#性别
data_0 = data1[which(data1$sex == 0),] 
data_1 = data1[which(data1$sex == 1),]  

#饮酒
data_0 = data1[which(data1$alcohol <14),] 
data_1 = data1[which(data1$alcohol >=14),]  



#Data1的cox回归
res.Cox <- coxph(Surv(Duration_ukb, termi)~
                   max_1_minute_avg+max_30_minute_avg_avg+metabolic_syd+
                   low_cadence_step_sum_avg+low_cadence_walk_minute_avg+mid_cadence_step_sum_avg+mid_cadence_walk_minute_avg+high_cadence_step_sum_avg+high_cadence_walk_minute_avg+
                   bicycling_overall_avg+cutPointMVPA_overall_avg+cutPointVPA_overall_avg+
                   sex+age+smoke+alcohol+TDI+ethnic_x, data = data1)
summary(res.Cox)

#分组的cox回归
res.Cox<- coxph(Surv(Duration_ukb, termi)~
                   step_sum_avg+max_1_minute_avg+max_30_minute_avg_avg+
                   sex+age+smoke+alcohol+TDI+ethnic_x, data = data_0)
summary(res.Cox0)

res.Cox1 <- coxph(Surv(Duration_ukb, termi)~
                    step_sum_avg+max_1_minute_avg+max_30_minute_avg_avg+
                    sex+age+smoke+alcohol+TDI+ethnic_x, data = data_1)
summary(res.Cox1)

res.Cox2 <- coxph(Surv(Duration_ukb, termi)~
                    step_sum_avg+max_1_minute_avg+max_30_minute_avg_avg+
                    sex+age+smoke+alcohol+TDI+ethnic_x, data = data_2)
summary(res.Cox2)

#定义HR计算函数
HR_fuc <- function(Data,X,xbin,Y,ybin){
                        x <- seq(min(Data[[X]]),max(Data[[X]]),xbin)
                        y <- seq(min(Data[[Y]]),max(Data[[Y]]), ybin)
                        result <- data.frame(matrix(NA, nrow = length(x), ncol = length(y))) # 创建一个空的数据框
                        dt1 <- data.frame(x = x)
                        dt2 <- data.frame(y=y)
                        x_HR <- exp(coefficients(res.Cox)[names(coefficients(res.Cox)) == X] * (dt1$x-mean(dt1$x)))
                        y_HR <- exp(coefficients(res.Cox)[names(coefficients(res.Cox)) == Y] * (dt2$y-mean(dt2$y)))
                        for (i in 1:nrow(result)) {
                             for (j in 1:ncol(result)) {
                                  result[i, j] <- x_HR[i]*y_HR[j]
                                     }
                        }
                        list(x = x, y = y, x_HR = x_HR,y_HR =y_HR,result = result)
}


#HR结果
#代谢病 sex 饮酒
HR_result_1 = HR_fuc(data_0,"max_30_minute_avg_avg",0.9,'step_sum_avg',200)
HR_result_2 = HR_fuc(data_1,"max_30_minute_avg_avg",0.9,'step_sum_avg',200)
write.csv(HR_result_1$y,"E:/QMY/Result/step_sum_avg_0_0414.csv", row.names = FALSE)
write.csv(HR_result_1$x,"E:/QMY/Result/max_30_minute_avg_avg0_0414.csv", row.names = FALSE)
write.csv(HR_result_1$result,"E:/QMY/Result/HR_result_1_0414.csv", row.names = FALSE)
write.csv(HR_result_2$y,"E:/QMY/Result/step_sum_avg_1_0414.csv", row.names = FALSE)
write.csv(HR_result_2$x,"E:/QMY/Result/max_30_minute_avg_avg1_0414.csv", row.names = FALSE)
write.csv(HR_result_2$result,"E:/QMY/Result/HR_result_2_0414.csv", row.names = FALSE)
#吸烟 年龄
HR_result_3 = HR_fuc(data_0,"max_30_minute_avg_avg",0.9,'step_sum_avg',200)
HR_result_4 = HR_fuc(data_1,"max_30_minute_avg_avg",0.9,'step_sum_avg',200)
HR_result_5 = HR_fuc(data_2,"max_30_minute_avg_avg",0.9,'step_sum_avg',200)
write.csv(HR_result_3$y,"E:/QMY/Result/step_sum_avg_0_0414.csv", row.names = FALSE)
write.csv(HR_result_3$x,"E:/QMY/Result/max_30_minute_avg_avg0_0414.csv", row.names = FALSE)
write.csv(HR_result_3$result,"E:/QMY/Result/HR_result_3_0414.csv", row.names = FALSE)
write.csv(HR_result_4$y,"E:/QMY/Result/step_sum_avg_1_0414.csv", row.names = FALSE)
write.csv(HR_result_4$x,"E:/QMY/Result/max_30_minute_avg_avg1_0414.csv", row.names = FALSE)
write.csv(HR_result_4$result,"E:/QMY/Result/HR_result_4_0414.csv", row.names = FALSE)
write.csv(HR_result_5$y,"E:/QMY/Result/step_sum_avg_2_0414.csv", row.names = FALSE)
write.csv(HR_result_5$x,"E:/QMY/Result/max_30_minute_avg_avg2_0414.csv", row.names = FALSE)
write.csv(HR_result_5$result,"E:/QMY/Result/HR_result_5_0414.csv", row.names = FALSE)

df<- expand.grid(x1 = HR_result_1$x,y1 = HR_result_1$y)
vec <- unlist(HR_result_1$result, use.names = FALSE)
Df<- cbind(df,vec)
df2<- expand.grid(x1 = HR_result_2$x,y1 = HR_result_2$y)
vec <- unlist(HR_result_2$result, use.names = FALSE)
Df2<- cbind(df2,vec)

library(patchwork)
p1 <- ggplot(Df, 
       aes(x1,y1, z = vec)) +
  geom_contour_filled() +
  labs(x = "Cadence_avg",
       y = "Steps_avg")
p2 <- ggplot(Df2, 
       aes(x1,y1, z = vec)) +
  geom_contour_filled() +
  labs(x = "Cadence_avg",
       y = "Steps_avg")
p1+p2+ plot_layout(nrow = 2)



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






data1<-read.table(file=r'{E:\QMY\Result\Data_0322.csv}', 
                  header=TRUE, sep=",")
res.Cox <- coxph(Surv(Duration_ukb, termi)~
                   avg_step_count_y+metabolic_syd+
                   sex+age+smoke+alcohol+TDI+ethnic_x+season_sin+season_cos, data = data1)
summary(res.Cox)


#将人按照步数分组
data1_q1 = data1[which(data1$avg_step_count_y <8800.3),] 
data1_q2 = data1[which(8800.3<= data1$avg_step_count_y & data1$avg_step_count_y <10913.8),] 
data1_q3 = data1[which(10913.8<= data1$avg_step_count_y & data1$avg_step_count_y <13164.2),] 
data1_q4 = data1[which(13164.2<= data1$avg_step_count_y),] 
#每组人的患病统计
table(data1_q1$termi)
table(data1_q2$termi)
table(data1_q3$termi)
table(data1_q4$termi)

#各个步数组与整体对比
res.Cox1 <- coxph(Surv(Duration_ukb, termi)~
                    Q1+metabolic_syd+
                    sex+age+smoke, data = data1)
summary(res.Cox1)
res.Cox2 <- coxph(Surv(Duration_ukb, termi)~
                    Q2+metabolic_syd+
                    sex+age+smoke, data = data1)
summary(res.Cox2)
res.Cox3 <- coxph(Surv(Duration_ukb, termi)~
                    Q3+metabolic_syd+
                    sex+age+smoke, data = data1)
summary(res.Cox3)
res.Cox4 <- coxph(Surv(Duration_ukb, termi)~
                    Q4+metabolic_syd+
                    sex+age+smoke, data = data1)
summary(res.Cox4)



#各个步数组与Q1对比
data1_q2= data1[which(data1$Q1 == 1 | data1$Q2 == 1 ),] #包含Q1和Q2的数据
data1_q3= data1[which(data1$Q1 == 1 | data1$Q3 == 1 ),] #包含Q1和Q3的数据
data1_q4= data1[which(data1$Q1 == 1 | data1$Q4 == 1 ),] #包含Q1和Q4的数据

res.Cox2_2 <- coxph(Surv(Duration_ukb, termi)~
                      Q2+metabolic_syd+
                      sex+age+smoke, data = data1_q2)
summary(res.Cox2_2)
res.Cox3_2 <- coxph(Surv(Duration_ukb, termi)~
                      Q3+metabolic_syd+
                      sex+age+smoke, data = data1_q3)
summary(res.Cox3_2)
res.Cox4_2 <- coxph(Surv(Duration_ukb, termi)~
                      Q4+metabolic_syd+
                      sex+age+smoke, data = data1_q4)
summary(res.Cox4_2)



#代谢分组
data_0 = data1[which(data1$metabolic_syd == 0),] 
data_1 = data1[which(data1$metabolic_syd == 1),] 

#无代谢病的cox回归
res.Cox <- coxph(Surv(Duration_ukb, termi)~
                    avg_step_count_y+
                    sex+age+smoke+alcohol+TDI+ethnic_x+season_sin+season_cos, data = data_0) #所有变量Cox回归
summary(res.Cox)

res.Cox <- coxph(Surv(Duration_ukb, termi)~
                   avg_step_count_y+
                   sex+age+alcohol, data = data_0) #显著变量Cox回归
summary(res.Cox)

#各个步数组与整体对比
res.Cox1 <- coxph(Surv(Duration_ukb, termi)~
                    Q1+
                    sex+age+alcohol, data = data_0)
summary(res.Cox1)
res.Cox2 <- coxph(Surv(Duration_ukb, termi)~
                    Q2+
                    sex+age+alcohol, data = data_0)
summary(res.Cox2)
res.Cox3 <- coxph(Surv(Duration_ukb, termi)~
                    Q3+
                    sex+age+alcohol, data = data_0)
summary(res.Cox3)
res.Cox4 <- coxph(Surv(Duration_ukb, termi)~
                    Q4+
                    sex+age+alcohol, data = data_0)
summary(res.Cox4)

#各个步数组与Q1对比
data_0_q2= data1[which(data1$Q1 == 1 | data1$Q2 == 1 ),] 
data_0_q3= data1[which(data1$Q1 == 1 | data1$Q3 == 1 ),] 
data_0_q4= data1[which(data1$Q1 == 1 | data1$Q4 == 1 ),] 

res.Cox2_2 <- coxph(Surv(Duration_ukb, termi)~
                      Q2+
                      sex+age+alcohol, data = data_0_q2)
summary(res.Cox2_2)
res.Cox3_2 <- coxph(Surv(Duration_ukb, termi)~
                      Q3+
                      sex+age+alcohol, data = data_0_q3)
summary(res.Cox3_2)
res.Cox4_2 <- coxph(Surv(Duration_ukb, termi)~
                      Q4+
                      sex+age+alcohol, data = data_0_q4)
summary(res.Cox4_2)










#0408
data1<-read.table(file=r'{E:\QMY\Result\Data0409.csv}', 
                  header=TRUE, sep=",")

res.Cox <- coxph(Surv(Duration_ukb, termi)~
                      feature1+feature2+feature3+feature4+
                   feature5+feature6+feature7+feature8+feature9, data = data1)
summary(res.Cox)

#0410
data1<-read.table(file=r'{E:\QMY\result\b_outside_0410.csv}', 
                  header=TRUE, sep=",")
res.Cox <- coxph(Surv(Duration_ukb, termi)~
                   a+b+c1+c2+c3+c4+c5+c6+c7+c8+c9+c10, data = data1)
summary(res.Cox)



















