install.packages("broom")
#森林图
install.packages("forestplot")
library(broom)
install.packages("XLConnect")
library("XLConnect")
install.packages("survminer") # 安装survminer包
install.packages("survival") # 安装survival包
library(survminer) # 加载包
library(survival) 




writeWorksheetToFile(r'{E:\QMY\Result\Result_Cox.xlsx}', 
                     data = res.Cox, 
                     sheet = "OLS", 
                     header = TRUE,
                     clearSheets = TRUE)


data1<-read.table(file=r'{E:\QMY\Result\Result01.csv}', 
                  header=TRUE, sep=",")
res.Cox <- coxph(Surv(min, termi)~
                 avg_step_count_y+
                 sex+age+smoke+alcohol+metabolic_syd, data = data1)
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


res.Cox <- coxph(Surv(min, termi)~
                   avg_step_count_y+
                   sex+age+smoke+alcohol+TDI+ethnic_x, data = data_0)
summary(res.Cox)


HR2<-read.table(file=r'{E:\QMY\Result\HR0319_2.csv}', 
               header=TRUE, sep=",",encoding="UTF8")

hr2=sprintf("%.3f",HR2$"HR")#获取HR列取小数点后3位
hrLow2=sprintf("%.3f",HR2$"HR.95L")#获取95%置信区间取HR.95L列小数点后3位
hrHigh2=sprintf("%.3f",HR2$"HR.95H")#获取95%置信区间HR.95H列小数点后3位
pVal=ifelse(HR2$PValue<0.001, "<0.001", sprintf("%.3f", HR2$PValue)) #获取P值
Hazard.ratio2=paste0(hr2,"(",hrLow2,"-",hrHigh2,")")


n=12
nRow=n+1
ylim=c(1,nRow)

split.screen(c(1,2))
screen(1)
xlim = c(0,4)
ylim=c(1,nRow)
par(mar=c(4,2,1.5,1.5))
print(plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab=""))
text.cex=0.8 
text(0,n:1,HR2$X,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',
                                           cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio2,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',
                                                 cex=text.cex,font=2,adj=1,)

screen(2)
xlim = c(0,2.5)
print(plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio"))
arrows(as.numeric(hrLow2),n:1,as.numeric(hrHigh2),n:1,angle=90,code=3,length=0.03,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr2) > 1, 'red', 'blue')
points(as.numeric(hr2), n:1, pch = 15, col = boxcolor, cex=1.5)
axis(1)





