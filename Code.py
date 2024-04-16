# -*- coding: utf-8 -*-

from lifelines.statistics import logrank_test, multivariate_logrank_test, pairwise_logrank_test
from lifelines import CoxPHFitter
from lifelines.utils import median_survival_times
from lifelines import KaplanMeierFitter
import pyreadr
import pandas as pd
import numpy as np
import pickle
import re
import matplotlib.pyplot as plt
from matplotlib import rcParams
from datetime import datetime
import math

plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False

#定义日期计算函数
def percentage(s):
    s = str(s)
    return int(s[:4])+int(s[4:len(s)])/1230.5



#导入表，字段处理
result = pyreadr.read_r(r'E:\QMY\data\select_attend_bak.RData') 
result2 = pyreadr.read_r(r'E:\QMY\data\acc_time.RData')
#result3 = pyreadr.read_r(r'C:\Users\Administrator\Documents\WeChat Files\wxid_1v7fkm60xfnu22\FileStorage\File\2023-03\avg_steps.RData')
result4 = pyreadr.read_r(r'E:\QMY\data\data_covariates.RData')
result5 = pyreadr.read_r(r'E:\QMY\data\death_time.RData')
result6 = pyreadr.read_r(r'E:\QMY\data\lost_follow_time.RData')
result7 = pyreadr.read_r(r'E:\QMY\data\cause_of_death.RData')
result8 = pyreadr.read_r(r'E:\QMY\data\ethnic_selected.RData')
result9 = pyreadr.read_r(r'E:\QMY\data\metabolic_new.RData')
result10 = pd.read_csv(r'E:\QMY\data\acc_subset.csv')
merged_result = pd.read_csv(r'E:\QMY\data\merged_result.csv')
merge2 = pd.read_csv(r'E:\QMY\data\merge2.csv')
result_0330 = pd.read_csv(r'E:\QMY\data\merged_result_v2(2).csv')
ffoutput = pd.read_csv(r'E:\QMY\data\ffoutput.csv')
metabolic_CO_T = pyreadr.read_r(r'E:\QMY\data\metabolic_CO_T.RData')
df_json_willetts_use = pd.read_csv(r'E:\QMY\data\df_json_willetts_use.csv') #加速度计数据
alcohol = pd.read_csv(r'E:\QMY\data\alcohol.csv')
data = pyreadr.read_r(r'E:\QMY\Data\complete_data_with_acc_circadian_0430.Rdata')
data2 = pyreadr.read_r(r'E:\QMY\Data\complete_data_0503.Rdata')



select_attend_bak = result["select_attend_bak"] #疾病发病时间
acc_time = result2["acc_time"] #加速度计时间
#avg_steps = result3["avg_steps"]
#avg_steps['eid'] = avg_steps['eid'].astype(int)
#avg_steps.rename(columns={'eid':'f.eid'},inplace=True)
data_covariates = result4["data_covariates"]
merged_result.rename(columns={'eid':'f.eid'},inplace=True)
death_time = result5['death_time']
lost_follow_time = result6['lost_follow_time']
cause_of_death = result7['cause_of_death']
ethnic_selected = result8['ethnic_selected']
#metabolic_new = result9['metabolic_new']
metabolic_CO_T = metabolic_CO_T['metabolic_CO_T']
metabolic_CO_T['eid'] = metabolic_CO_T['eid'].astype(int)
metabolic_CO_T.rename(columns={'eid':'f.eid'},inplace=True)
data = data['complete_data_with_acc_circadian_0430']
data2 = data2['complete_data_0503']


'''确定关键时间点'''
#疾病字段代码
illness = ['f.131296.0.0',
'f.131298.0.0', 'f.131300.0.0', 'f.131302.0.0', 'f.131304.0.0',
'f.131306.0.0', 'f.131360.0.0', 'f.131362.0.0',
'f.131364.0.0', 'f.131366.0.0', 'f.131368.0.0']

ill_code = ['I21','I22','I23','I24','I25','I60','I61','I62','I63','I64'] #疾病代码


#连接表
merge1 = pd.merge(select_attend_bak,acc_time,on='f.eid')
merge1 = pd.merge(merge1,lost_follow_time,on='f.eid')

#开始时间
merge1['f.90010.0.0_y'] = pd.to_datetime(merge1['f.90010.0.0_y']).dt.strftime('%Y%m%d')#加速度计时间
merge1['f.53.0.0'] = pd.to_datetime(merge1['f.53.0.0']).dt.strftime('%Y%m%d') #ukb时间

#各个疾病的时间转换
for i in illness:
    merge1[i] = pd.to_datetime(merge1[i]).dt.strftime('%Y%m%d')
merge1 = merge1.replace(np.nan,99999999)
merge1['f.90010.0.0_y'] = merge1['f.90010.0.0_y'].astype(int)
merge1['f.53.0.0'] = merge1['f.53.0.0'].astype(int)
merge1[illness] = merge1[illness].astype(int)
merge1['terminate'] = 20221022 #终点时间

#疾病第一次被诊断的时间
merge1['min'] = merge1[illness].min(axis=1)

#死亡时间
merge2 = pd.merge(merge1,death_time,on='f.eid')
merge2['f.40000.0.0'] = pd.to_datetime(merge2['f.40000.0.0']).dt.strftime('%Y%m%d')
merge2['f.40000.0.0'] = merge2['f.40000.0.0'].replace(np.nan,99999999)
merge2['f.40000.0.0'] = merge2['f.40000.0.0'].astype(int)



'''确定终点事件'''
merge3 = pd.merge(merge2,cause_of_death,on='f.eid')  
merge3[illness] = merge3[illness].astype(int)

#疾病第一次被诊断的时间
merge3['min'] = merge3[illness].min(axis=1)
merge3['f.40001.0.0'] = merge3['f.40001.0.0'].str[:3] #疾病类型：截取字符串
merge3['f.40001.0.0'].value_counts()
merge3['termi'] = np.zeros(len(merge3))

for i in range(len(merge3)):
    if merge3['min'][i] != 99999999: #有最早的病时间
        merge3['termi'][i] = 1
    if merge3['min'][i] == 99999999:
        if merge3['f.40001.0.0'][i] in ill_code:
            merge3['min'][i] = merge3['f.40000.0.0'][i] #最早得病时间=死亡时间
            merge3['termi'][i] = 1
        else:
            merge3['termi'][i] = 0



'''确定随访时长'''
merge3['随访时长'] =0
#merge1.query('min==20121018').index
#使用90010计算
for i in range(len(merge3)):
    if merge3.iloc[i]['f.90010.0.0_y'] == 99999999: #没参加
       merge3['随访时长'].iloc[i] = np.nan
    if merge3.iloc[i]['f.90010.0.0_y'] != 99999999 and merge3.iloc[i]['min'] != 99999999: #疾病被诊断
        merge3['随访时长'].iloc[i] = percentage(merge3.iloc[i]['min'])-percentage(merge3.iloc[i]['f.90010.0.0_y'])
    if merge3.iloc[i]['f.90010.0.0_y'] != 99999999 and merge3.iloc[i]['min'] == 99999999 and merge3.iloc[i]['f.40000.0.0'] != 99999999: #人员死亡
        merge3['随访时长'].iloc[i] = percentage(merge3.iloc[i]['f.40000.0.0'])-percentage(merge3.iloc[i]['f.90010.0.0_y'])
    if merge3.iloc[i]['f.90010.0.0_y'] != 99999999 and merge3.iloc[i]['min'] == 99999999 and merge3.iloc[i]['f.40000.0.0'] == 99999999:
        merge3['随访时长'].iloc[i] = percentage(merge3.iloc[i]['terminate'])-percentage(merge3.iloc[i]['f.90010.0.0_y'])
merge3['随访时长'].value_counts()

#使用53计算
merge3['Duration_ukb'] =0
for i in range(len(merge3)):
    if merge3.iloc[i]['f.53.0.0'] == 99999999: #没参加
       merge3['Duration_ukb'].iloc[i] = np.nan
    if merge3.iloc[i]['f.53.0.0'] != 99999999 and merge3.iloc[i]['min'] != 99999999: #疾病被诊断
        merge3['Duration_ukb'].iloc[i] = percentage(merge3.iloc[i]['min'])-percentage(merge3.iloc[i]['f.53.0.0'])
    if merge3.iloc[i]['f.53.0.0'] != 99999999 and merge3.iloc[i]['min'] == 99999999 and merge3.iloc[i]['f.40000.0.0'] != 99999999: #人员死亡
        merge3['Duration_ukb'].iloc[i] = percentage(merge3.iloc[i]['f.40000.0.0'])-percentage(merge3.iloc[i]['f.53.0.0'])
    if merge3.iloc[i]['f.53.0.0'] != 99999999 and merge3.iloc[i]['min'] == 99999999 and merge3.iloc[i]['f.40000.0.0'] == 99999999:
        merge3['Duration_ukb'].iloc[i] = percentage(merge3.iloc[i]['terminate'])-percentage(merge3.iloc[i]['f.53.0.0'])
merge3['Duration_ukb'].value_counts()

#随访时长分布
def distribution(column,Bin,Range,title):
    rcParams['axes.titlepad'] = 30 #标题与图间距
    fig,ax = plt.subplots(figsize=(15,10), dpi=300) 
    #plt.hist(normal['retire'],40,range=(-20,20))
    plt.hist(column,Bin,range=Range)
    ax.get_yaxis().get_major_formatter().set_scientific(False) #关掉y轴科学计数法
    plt.xticks(fontsize=26)
    plt.yticks(fontsize=26)
    plt.xlabel('随访时长（年）',fontsize=32)
    plt.ylabel('人数',fontsize=32)
    plt.title(title,fontsize=40)
#distribution(merge3['随访时长'],50,(0,11),'随访时长(f.90010.0.0)')
distribution(merge3['Duration_ukb'],50,(0,16),'随访时长(f.53.0.0)')

merge3.to_csv(r'E:\QMY\data\Terminate.csv',index=False,sep=',')


'''def save_variable(v,filename):
    f=open(filename,'wb')
    pickle.dump(v,f,0)
    f.close()
    return filename
filename = save_variable(merge3['随访时长'],r'E:\QMY\data\merge3.csv')'''











'''纳排标准'''
merge4 = pd.read_csv(r'E:\QMY\data\Terminate.csv')
merge4 = pd.merge(merge4,data_covariates,on='f.eid')
merge4 = pd.merge(merge4,ethnic_selected,on='f.eid') 
merge4 = pd.merge(merge4,acc_time,on='f.eid')
merge4 = pd.merge(merge4,metabolic_CO_T,on='f.eid') #48097
merge4 = pd.merge(merge4,alcohol,on='f.eid')

#排除参加访问前患卒中、冠心病的人群
'''merge4 = merge4[~merge4['f.53.0.0'].isin([99999999])]
merge4['f.53.0.0'] = pd.to_datetime(merge4['f.53.0.0']).dt.strftime('%Y%m%d')
merge4 = merge4.dropna(subset=['f.53.0.0'])
merge4['f.53.0.0'] = merge4['f.53.0.0'].astype(int)
merge4 = merge4[merge4['f.53.0.0'] < merge4['min']] '''

#排除异常值与缺失值
#merge4 = pd.merge(merge4,df_json_willetts_use,on='f.eid')
#merge4 = pd.merge(merge4,avg_steps,on='f.eid')
#merge4 = pd.merge(merge4,merged_result,on='f.eid') 
merge4.isnull().sum()
merge4 = merge4.dropna(subset=['ethnic_x','smoke','TDI']) #62365

#饮酒
merge4['alcohol'] = merge4['alcohol_weekly_volume'] + merge4['alcohol_monthly_volume']/4

#数据描述
rcParams['axes.titlepad'] = 15
plt.hist(merge4['avg_step_count_y'],100)
plt.title('日均步数',fontsize=20)
plt.boxplot(merge4['age']) #年龄
plt.title('年龄',fontsize=20)
plt.boxplot(merge4['alcohol']) #饮酒
plt.title('饮酒',fontsize=20)
plt.boxplot(merge4['TDI']) #汤森
plt.title('TDI',fontsize=20)
merge4['sex'].value_counts() #性别
merge4['sex'] = merge4['sex'].astype('float16')
merge4['smoke'].value_counts() #吸烟
merge4['smoke'] = merge4['smoke'].astype('float16') #类型转换，让corr有值
merge4['ethnic_selected'].value_counts() #民族
merge4['ethnic_x'] = merge4['ethnic_x'].astype('float16')
merge4['amp'].value_counts()
merge4['have_CO'].value_counts()

#排除协变量异常值
merge4 = merge4.query('smoke == 0 | smoke == 1 | smoke == 2') #排除抽烟 
merge4 = merge4.query('ethnic_selected != "others/unknown"') #排除民族 
merge4['alc'] = merge4['f.20117.0.0']
merge4 = merge4.query('alc != -3') #排除饮酒
merge4.isnull().sum()
#merge4 = merge4.dropna(subset = ['weekday_avg_step_count','weekend_avg_step_count',
                     #'MET-overall-sd','walking-overall-avg','walking-overall-sd']) #删掉缺失值 9879
#merge4.rename(columns={'quality-goodWearTime':'qualitygoodWearTime',
                                #'quality-goodCalibration':'qualitygoodCalibration'},inplace=True)
#merge4 = merge4.query('qualitygoodWearTime ==1 & qualitygoodCalibration ==1') #9776

walmsey = pd.read_csv(r'E:\QMY\data\walmsley_all\walmsley_1_2.csv') #anti-logic模型
merge4 = pd.merge(merge4,walmsey,on='f.eid')


#季节因素转化
merge4 = merge4.reset_index(drop=True)
#确定日期对应每年第几天
merge4['day'] = np.zeros(len(merge4))
for i in range(len(merge4)):
    merge4['day'][i] = (int(str(merge4['f.90010.0.0_y'][i])[4:6])-1)*30.5+int(str(merge4['f.90010.0.0_y'][i])[6:])-0.5
#确定季节协变量
merge4['season_sin'] = merge4['season_cos'] = np.zeros(len(merge4))
for i in range(len(merge4)):
    merge4['season_sin'][i] = math.sin(2*math.pi*merge4['day'][i]/365)
    merge4['season_cos'][i] = math.cos(2*math.pi*merge4['day'][i]/365)

merge4.to_csv(r'E:\QMY\result\Data428.csv',index=False,sep=',')

#步数与终点事件
plt.hist(merge4['avg_step_count_y'])
merge4.corr(method='spearman')['termi']
merge4['step_interval']= pd.cut(merge4['avg_step_count_y'],bins=50)
fig,ax = plt.subplots(figsize=(20, 10), dpi=150)
ax.plot(merge4.groupby('step_interval')['termi'].mean().values,'-ok', color='green',
         markersize=8, linewidth=6,
         markerfacecolor='black',
         markeredgewidth=2)
plt.xticks(size=30)
plt.yticks(size=30)
plt.xlim(4,30)
plt.xlabel('步数',fontsize=30)
plt.ylabel('产生终点事件比例',fontsize=30)
plt.title('步数与终点事件',fontsize=40)





cols1 = ['minimum_1', 'minimum_2', 'minimum_3', 'minimum_4', 'minimum_5', 'minimum_6', 'minimum_7']
cols2 = ['amp', 'amp_2', 'amp_3', 'amp_4', 'amp_5', 'amp_6', 'amp_7']
cols3 = ['alpha_1', 'alpha_2', 'alpha_3', 'alpha_4', 'alpha_5', 'alpha_6', 'alpha_7']
cols4 = ['beta_1', 'beta_2', 'beta_3', 'beta_4', 'beta_5', 'beta_6', 'beta_7']
cols5 = ['acrotime_1', 'acrotime_2', 'acrotime_3', 'acrotime_4', 'acrotime_5', 'acrotime_6', 'acrotime_7']
cols6 = ['F_pseudo_1', 'F_pseudo_2', 'F_pseudo_3', 'F_pseudo_4', 'F_pseudo_5', 'F_pseudo_6', 'F_pseudo_7']
cols7 = ['UpMesor_1', 'UpMesor_2', 'UpMesor_3', 'UpMesor_4', 'UpMesor_5', 'UpMesor_6', 'UpMesor_7']
cols8 = ['DownMesor_1', 'DownMesor_2', 'DownMesor_3', 'DownMesor_4', 'DownMesor_5', 'DownMesor_6', 'DownMesor_7']
cols9 = ['MESOR_1', 'MESOR_2', 'MESOR_3', 'MESOR_4', 'MESOR_5', 'MESOR_6', 'MESOR_7']
#Data424['minimum_sd7'] = np.std(walmsey[cols1],axis=1)






'''生存分析'''
merge4 = pd.read_csv(r'E:\QMY\Result\Data_0321.csv')
#生存函数：总体
def Kmf(table):
    kmf = KaplanMeierFitter() 
    kmf.fit(table['Duration_ukb'], event_observed=table['termi'])
    kmf.plot_survival_function()
    median_ = kmf.median_survival_time_
    median_confidence_interval_ = median_survival_times(kmf.confidence_interval_)
    print(median_confidence_interval_)

#分组生存函数
#代谢不作为分组
#步数分组：
def kmf_steps_4(table,ylim):
    kmf = KaplanMeierFitter()
    rcParams['axes.titlepad'] = 30 
    plt.rcParams.update({'font.size':25})
    plt.figure(figsize=(12,8), dpi=80) 
    ix1 = table.query('Q1 ==1').index
    ix2 = table.query('Q2 ==1').index
    ix3 = table.query('Q3 ==1').index
    ix4 = table.query('Q4 ==1').index
    kmf.fit(table['Duration_ukb'][ix1], table['termi'][ix1], label='step_Q1')
    ax = kmf.plot()
    kmf.fit(table['Duration_ukb'][ix2], table['termi'][ix2], label='step_Q2')
    ax = kmf.plot(ax=ax)
    kmf.fit(table['Duration_ukb'][ix3], table['termi'][ix3], label='step_Q3')
    ax = kmf.plot(ax=ax)
    kmf.fit(table['Duration_ukb'][ix4], table['termi'][ix4], label='step_Q4')
    ax = kmf.plot(ax=ax)
    plt.xlim(0,14)
    plt.ylim(ylim,1.0)
    plt.text(0.2,ylim+0.01,'P_value<0.005')
    plt.xlabel('随访时长（年）',rotation=0,fontsize=25)
    plt.ylabel('Survival Probability',rotation=90,fontsize=25)
    plt.title('生存函数',fontsize=32)


#不同代谢病分组
kmf = KaplanMeierFitter()
rcParams['axes.titlepad'] = 30 
plt.rcParams.update({'font.size':25})
plt.figure(figsize=(12,8), dpi=80) 
ix1 = merge4.query('alpha<0.2').index
ix2 = merge4.query('metabolic_syd ==1.0').index
kmf.fit(merge4['Duration_ukb'][ix1], merge4['termi'][ix1], label='无代谢综合征')
ax = kmf.plot()
kmf.fit(merge4['Duration_ukb'][ix2], merge4['termi'][ix2], label='有代谢综合征')
ax = kmf.plot(ax=ax)
plt.xlim(0,15)
plt.ylim(0.84,1.0)
plt.xlabel('随访时长（年）',rotation=0,fontsize=25)
plt.ylabel('Survival Probability',rotation=90,fontsize=25)
plt.title('生存函数',fontsize=32)
#log-rank test
results = logrank_test(merge4['Duration_ukb'][ix1], merge4['Duration_ukb'][ix2], merge4['termi'][ix1], 
                       merge4['termi'][ix2], alpha=.99)
results.print_summary()
print("p_value:{}".format(results.p_value))

#Sex
rcParams['axes.titlepad'] = 30 
plt.rcParams.update({'font.size':25})
plt.figure(figsize=(12,8), dpi=80) 
ix1 = merge4.query('sex ==0.0').index
ix2 = merge4.query('sex ==1.0').index
kmf.fit(merge4['Duration_ukb'][ix1], merge4['termi'][ix1], label='女性')
ax = kmf.plot()
kmf.fit(merge4['Duration_ukb'][ix2], merge4['termi'][ix2], label='男性')
ax = kmf.plot(ax=ax)
plt.xlim(0,15)
plt.ylim(0.86,1.0)
plt.xlabel('随访时长（年）',rotation=0,fontsize=25)
plt.ylabel('Survival Probability',rotation=90,fontsize=25)
plt.title('生存函数',fontsize=32)
#log-rank test
results = logrank_test(merge4['Duration_ukb'][ix1], merge4['Duration_ukb'][ix2], merge4['termi'][ix1], 
                       merge4['termi'][ix2], alpha=.99)
results.print_summary()
print("p_value:{}".format(results.p_value))




#加入代谢病分组后
rcParams['axes.titlepad'] = 30 
plt.rcParams.update({'font.size':25})
plt.figure(figsize=(12,8), dpi=80) 

ix1 = merge4.query('metabolic_syd == 0.0 & steps_over_1w ==0.0').index
ix2 = merge4.query('metabolic_syd == 0.0 & steps_over_1w ==1.0').index
ix3 = merge4.query('metabolic_syd == 1.0 & steps_over_1w ==0.0').index
ix4 = merge4.query('metabolic_syd == 1.0 & steps_over_1w ==1.0').index
'''kmf.fit(merge4['Duration_ukb'][ix1], merge4['termi'][ix1], label='无代谢综合征且日均步数不到1万步')
ax = kmf.plot()
kmf.fit(merge4['Duration_ukb'][ix2], merge4['termi'][ix2], label='无代谢综合征且日均步数超过1万步')
ax = kmf.plot()'''
kmf.fit(merge4['Duration_ukb'][ix3], merge4['termi'][ix3], label='有代谢综合征且日均步数不到1万步')
ax = kmf.plot()
kmf.fit(merge4['Duration_ukb'][ix4], merge4['termi'][ix4], label='有代谢综合征且日均步数超过1万步')
ax = kmf.plot(ax=ax)
plt.xlim(0,15)
plt.ylim(0.82,1.0)
plt.xlabel('随访时长（年）',rotation=0,fontsize=25)
plt.ylabel('Survival Probability',rotation=90,fontsize=25)
plt.title('生存函数',fontsize=32)

#log-rank test
merge4['class'] = np.zeros(len(merge4))
merge4['class'][ix1] = '无代谢综合征且日均步数不到1万步'
merge4['class'][ix2] = '无代谢综合征且日均步数超过1万步'
merge4['class'][ix3] = '有代谢综合征且日均步数不到1万步'
merge4['class'][ix4] = '有代谢综合征且日均步数超过1万步'
results = pairwise_logrank_test(
        event_durations = merge4['Duration_ukb'],
        groups=merge4['class'],
        event_observed=merge4['termi']
)
summary = results.print_summary()
print("p_value:{}".format(results.p_value))

#不同代谢病和步数情况饼图
plt.pie(x=[len(ix1),len(ix2),len(ix3),len(ix4)],labels=['无代谢综合征\n且日均步数不到1万步',
        '无代谢综合征\n且日均步数超过1万步','有代谢综合征\n且日均步数不到1万步','有代谢综合征\n且日均步数超过1万步'],
        shadow=True, autopct="%0.2f%%")






#年龄分箱
rcParams['axes.titlepad'] = 30 
plt.rcParams.update({'font.size':25})
plt.figure(figsize=(12,8), dpi=80) 

ix1 = merge4.query('40<=age<50').index
ix2 = merge4.query('50<=age<60').index
ix3 = merge4.query('60<=age<=70').index

kmf = KaplanMeierFitter()
kmf.fit(merge4['Duration_ukb'][ix1], merge4['termi'][ix1], label='40-50岁人群')
ax = kmf.plot()
kmf.fit(merge4['Duration_ukb'][ix2], merge4['termi'][ix2], label='50-60岁人群')
ax = kmf.plot()
kmf.fit(merge4['Duration_ukb'][ix3], merge4['termi'][ix3], label='60-70岁人群')
ax = kmf.plot(ax=ax)
plt.xlabel('随访时长（年）',rotation=0,fontsize=25)
plt.ylabel('Survival Probability',rotation=90,fontsize=25)
plt.title('生存函数',fontsize=32)

merge4['class'] = np.zeros(len(merge4))
merge4['class'][ix1] = '40-50岁人群'
merge4['class'][ix2] = '50-60岁人群'
merge4['class'][ix3] = '60-70岁人群'
results = multivariate_logrank_test(
        event_durations = merge4['Duration_ukb'],
        groups=merge4['class'],
        event_observed=merge4['termi']
)
summary = results.print_summary()
print("p_value:{}".format(results.p_value))


#alcohol分箱
merge4['alcohol'].mean()
ix1 = merge4.query('alcohol<14').index
ix2 = merge4.query('alcohol>=14').index
plt.pie(x=[len(ix1),len(ix2)],labels=['未过量饮酒','过量饮酒'],
        shadow=True, autopct="%0.2f%%")

rcParams['axes.titlepad'] = 30 
plt.rcParams.update({'font.size':25})
plt.figure(figsize=(12,8), dpi=80) 
kmf.fit(merge4['Duration_ukb'][ix1], merge4['termi'][ix1], label='未过量饮酒（>14杯/周）')
ax = kmf.plot()
kmf.fit(merge4['Duration_ukb'][ix2], merge4['termi'][ix2], label='过量饮酒（>14杯/周）')
ax = kmf.plot(ax=ax)
plt.xlim(0,15.5)
plt.ylim(0.86,1.0)
plt.title('生存函数',fontsize=25)
#log-rank test
results = logrank_test(merge4['Duration_ukb'][ix1], merge4['Duration_ukb'][ix2], merge4['termi'][ix1], 
                       merge4['termi'][ix2], alpha=.99)
results.print_summary()
print("p_value:{}".format(results.p_value))


#alcohol+代谢
ix1 = merge4.query('steps_over_1w == 0.0 & alcohol<14').index
ix2 = merge4.query('steps_over_1w == 0.0 &alcohol>=14').index
ix3 = merge4.query('metabolic_syd == 1.0 & alcohol<14').index
ix4 = merge4.query('metabolic_syd == 1.0 &alcohol>=14').index
merge4['class'] = np.zeros(len(merge4))
merge4['class'][ix1] = '日均步数不到1万步且未过量饮酒'
merge4['class'][ix2] = '日均步数不到1万步且过量饮酒'
merge4['class'][ix3] = '日均步数超过1万步且未过量饮酒'
merge4['class'][ix4] = '日均步数超过1万步且过量饮酒'

rcParams['axes.titlepad'] = 30 
plt.rcParams.update({'font.size':25})
plt.figure(figsize=(12,8), dpi=80) 
'''kmf.fit(merge4['Duration_ukb'][ix1], merge4['termi'][ix1], label='日均步数不到1万步且未过量饮酒')
ax = kmf.plot()
kmf.fit(merge4['Duration_ukb'][ix2], merge4['termi'][ix2], label='日均步数不到1万步且过量饮酒')
ax = kmf.plot()'''
kmf.fit(merge4['Duration_ukb'][ix3], merge4['termi'][ix3], label='日均步数超过1万步且未过量饮酒')
ax = kmf.plot()
kmf.fit(merge4['Duration_ukb'][ix4], merge4['termi'][ix4], label='日均步数超过1万步且过量饮酒')
ax = kmf.plot()
plt.xlim(0,15.5)
plt.ylim(0.85,1.0)
plt.title('生存函数',fontsize=25)

results = logrank_test(merge4['Duration_ukb'][ix3], merge4['Duration_ukb'][ix4], merge4['termi'][ix3], 
                       merge4['termi'][ix4], alpha=.99)
results.print_summary()
print("p_value:{}".format(results.p_value))



#smoke
rcParams['axes.titlepad'] = 30 
plt.rcParams.update({'font.size':25})
plt.figure(figsize=(12,8), dpi=80) 

ix1 = merge4.query('smoke==0.0').index
ix2 = merge4.query('smoke==1.0').index
ix3 = merge4.query('smoke==2.0').index

kmf.fit(merge4['Duration_ukb'][ix1], merge4['termi'][ix1], label='从不吸烟')
ax = kmf.plot()
kmf.fit(merge4['Duration_ukb'][ix2], merge4['termi'][ix2], label='曾经吸烟')
ax = kmf.plot()
kmf.fit(merge4['Duration_ukb'][ix3], merge4['termi'][ix3], label='现在吸烟')
ax = kmf.plot(ax=ax)
plt.xlim(0,15.5)
plt.ylim(0.86,1.0)
plt.xlabel('随访时长（年）',rotation=0,fontsize=25)
plt.ylabel('Survival Probability',rotation=90,fontsize=25)
plt.title('生存函数',fontsize=32)

merge4['class'] = np.zeros(len(merge4))
merge4['class'][ix1] = '从不吸烟'
merge4['class'][ix2] = '曾经吸烟'
merge4['class'][ix3] = '现在吸烟'
results = pairwise_logrank_test(
        event_durations = merge4['Duration_ukb'],
        groups=merge4['class'],
        event_observed=merge4['termi']
)
summary = results.print_summary()
print("p_value:{}".format(results.p_value))



#TDI
merge4['TDI'].mean()
merge4['TDI'].median()

rcParams['axes.titlepad'] = 30 
plt.rcParams.update({'font.size':25})
plt.figure(figsize=(12,8), dpi=80) 

ix1 = merge4.query('TDI<-2.5').index
ix2 = merge4.query('TDI>=-2.5').index

kmf.fit(merge4['Duration_ukb'][ix1], merge4['termi'][ix1], label='TDI<-2.5')
ax = kmf.plot()
kmf.fit(merge4['Duration_ukb'][ix2], merge4['termi'][ix2], label='TDI>=-2.5')
ax = kmf.plot(ax=ax)
plt.xlim(0,15.5)
plt.ylim(0.88,1.0)
plt.xlabel('随访时长（年）',rotation=0,fontsize=25)
plt.ylabel('Survival Probability',rotation=90,fontsize=25)
plt.title('生存函数',fontsize=32)

results = logrank_test(merge4['Duration_ukb'][ix1], merge4['Duration_ukb'][ix2], merge4['termi'][ix1], 
                       merge4['termi'][ix2], alpha=.99)
results.print_summary()
print("p_value:{}".format(results.p_value))



#种族
merge4['ethnic_selected'].value_counts()
rcParams['axes.titlepad'] = 30 
plt.rcParams.update({'font.size':25})
plt.figure(figsize=(12,8), dpi=80) 

ix1 = merge4.query('ethnic_selected=="White"').index
ix2 = merge4.query('ethnic_selected!="White"').index

kmf.fit(merge4['Duration_ukb'][ix1], merge4['termi'][ix1], label='白人')
ax = kmf.plot()
kmf.fit(merge4['Duration_ukb'][ix2], merge4['termi'][ix2], label='非白人')
ax = kmf.plot(ax=ax)
plt.xlim(0,15.5)
plt.ylim(0.88,1.0)
plt.xlabel('随访时长（年）',rotation=0,fontsize=25)
plt.ylabel('Survival Probability',rotation=90,fontsize=25)
plt.title('生存函数',fontsize=32)

results = logrank_test(merge4['Duration_ukb'][ix1], merge4['Duration_ukb'][ix2], merge4['termi'][ix1], 
                       merge4['termi'][ix2], alpha=.99)
results.print_summary()
print("p_value:{}".format(results.p_value))


merge4 = pd.read_csv(r'E:\QMY\Result\Data_0321.csv')
merge4['Q1'] = merge4['Q2'] = merge4['Q3'] = merge4['Q4'] = merge4['other'] = np.zeros(len(merge4))
ix1 = merge4.query('avg_step_count_y <8800.3').index
ix2 = merge4.query('8800.3<= avg_step_count_y<10913.8').index
ix3 = merge4.query('10913.8<= avg_step_count_y<13164.2').index
ix4 = merge4.query('13164.2<= avg_step_count_y').index
#other = merge4.query('13164.2<= avg_step_count_y<13165').index
merge4['Q1'][ix1] = 1
merge4['Q2'][ix2] = 1
merge4['Q3'][ix3] = 1
merge4['Q4'][ix4] = 1

merge4.to_csv(r'E:\QMY\result\Data_0322.csv',index=False,sep=',')


'''划分：有无代谢病'''
merge_0322 = pd.read_csv(r'E:\QMY\Result\Data_0322.csv')
merge_0322_0 = merge_0322.query('metabolic_syd == 0')
merge_0322_1 = merge_0322.query('metabolic_syd == 1')

print(Data428['alpha'].min(),
      merge_0322_0['avg_step_count_y'].quantile(0.25),
      merge_0322_0['avg_step_count_y'].median(),
      merge_0322_0['avg_step_count_y'].quantile(0.75),
      merge_0322_0['avg_step_count_y'].max())

merge_0322_0['Q1'] = merge_0322_0['Q2'] = merge_0322_0['Q3'] = merge_0322_0['Q4'] = merge_0322_0['other'] = np.zeros(len(merge_0322_0))
ix1 = merge_0322_0.query('avg_step_count_y <9275.893').index
ix2 = merge_0322_0.query('9275.893<= avg_step_count_y<=11389.017').index
ix3 = merge_0322_0.query('11389.017< avg_step_count_y<13565.807855').index
ix4 = merge_0322_0.query('13565.807855<= avg_step_count_y').index
#other = merge_0322_0.query('13164.2<= avg_step_count_y<13165').index
merge_0322_0['Q1'][ix1] = 1
merge_0322_0['Q2'][ix2] = 1
merge_0322_0['Q3'][ix3] = 1
merge_0322_0['Q4'][ix4] = 1

print(merge_0322_0['termi'][ix1].value_counts(),
merge_0322_0['termi'][ix2].value_counts(),
merge_0322_0['termi'][ix3].value_counts(),
merge_0322_0['termi'][ix4].value_counts())




'''0324 数据描述'''
merged_json_walmsley = pd.read_csv(r'E:\QMY\data\merged_json_walmsley.csv')
merged_json_willetts = pd.read_csv(r'E:\QMY\data\merged_json_willetts.csv')
walmsley_dictionary = pd.read_csv(r'E:\QMY\data\walmsley_dictionary.csv')
steps_dictionary = pd.read_csv(r'E:\QMY\data\steps_dictionary_0330.csv')
minute_steps = pd.read_csv(r'E:\QMY\data\minute_steps(2).csv') #一个人7天的分钟级步数情况
minute_steps = minute_steps.query('C1 !=0 ') #把步数为0的去掉后的步数分布
minute_steps = minute_steps.query('C2 !=0 ')
plt.hist(minute_steps['C2'],100)

#定义直方图函数
def hist(table,column):
    rcParams['axes.titlepad'] = 30 
    plt.rcParams.update({'font.size':25})
    plt.figure(figsize=(12,8), dpi=80) 
    plt.hist(table[column],50)
    plt.xlim(table[column].min()*0.95,table[column].max()*1.02)
    plt.title(column,fontsize=32)

#批量绘制直方图
for i in range(len(steps_dictionary['eid'])):
    hist(result_0330,steps_dictionary['eid'].values[i])
    plt.savefig('E:/QMY/result/Hist_result/Hist(steps,cadence)/{}.jpg'.format(steps_dictionary['eid'].values[i]))


plt.scatter(result_0330['daily_steps_avg'],result_0330['daily_max_30min_cadences_avg'])
#plt.title('步数与30分钟峰值节奏',fontsize=32)

plt.scatter(result_0330['daily_steps_avg'],result_0330['daily_minute_cadence_avgs_avg'])
plt.title('步数与分钟级步频',fontsize=32)

plt.scatter(result_0330['daily_max_30min_cadences_avg'],result_0330['daily_steps_avg']/60/24)
plt.title('步数与分钟级步频',fontsize=32)


'''4.8 频域特征分析'''
ffoutput = pd.read_csv(r'E:\QMY\data\ffoutput(1).csv')
merge_0322 = pd.read_csv(r'E:\QMY\Result\Data_0321.csv')

terminate = pd.read_csv(r'E:\QMY\data\terminate.csv')
terminate.fillna(-999, inplace=True)
terminate['f.53.0.0'] = pd.to_datetime(terminate['f.53.0.0']).dt.strftime('%Y%m%d')
terminate['f.53.0.0'] = terminate['f.53.0.0'].astype(int)
terminate['min'] = terminate['min'].astype(int)
terminate['terminate'] = 20221022
terminate['Duration_ukb'] = np.zeros(len(terminate))
for i in range(len(terminate)):
    if terminate['termi'][i] == 1: #有最早的病时间
        terminate['Duration_ukb'][i] = percentage(terminate['min'][i]) - percentage(terminate['f.53.0.0'][i])
    else:
        terminate['Duration_ukb'][i] = percentage(terminate['terminate'][i]) - percentage(terminate['f.53.0.0'][i])

terminate.to_csv(r'E:\QMY\data\Terminate.csv',index=False,sep=',')

ffoutput_0409 = pd.read_csv(r'E:\QMY\data\ffoutput0409.csv') #导入不同b的数据
ffoutput_0409[['c1','c2','c3',
               'c4','c5','c6','c7']] = ffoutput_0409[['c1','c2','c3',
                              'c4','c5','c6','c7']].abs() #将各个c取绝对值

merge_0409 = pd.merge(ffoutput_0409,terminate,on='f.eid')
merge_0409.to_csv(r'E:\QMY\result\Data0409.csv',index=False,sep=',')




'''4.10 步数各个指标的表 作人员关于佩戴时间的筛选
级数数据cox回归'''
merged_result_0409 = pd.read_csv(r'E:\QMY\data\merged_result_v2(2).csv')
merged_json_willetts = pd.read_csv(r'E:\QMY\data\merged_json_willetts.csv')
#处理f.eid
def clean_feid(df):
    # 遍历 f.eid 列的所有值，并进行处理
    rows_to_drop = []  # 需要删除的行索引列表
    for i in range(len(df)):
        feid = df.loc[i, 'f.eid']
        # 使用正则表达式匹配数字部分
        match = re.search(r'\d{7}$', feid)
        if match:
            # 如果本身就是7位数字，则保留原状
            feid_cleaned = match.group()
        else:
            # 否则，提取出最后的7位数字
            match = re.search(r'\d{7}$', feid.split('/')[-1])
            feid_cleaned = match.group() if match else ''
        
        # 判断是否以"cwa"结尾，如果是，则删除整行
        if feid_cleaned.endswith('cwa'):
            rows_to_drop.append(i)
        else:
            # 更新回 merged_json_willetts 表中对应的记录
            df.loc[i, 'f.eid'] = feid_cleaned

    # 删除需要删除的行
    df.drop(index=rows_to_drop, inplace=True)

    return df

# 清洗 f.eid 列的数据
merged_json_willetts = clean_feid(merged_json_willetts)
merged_json_willetts['f.eid'] = pd.to_numeric(merged_json_willetts['f.eid'], errors='coerce', downcast='integer')
# 删除包含非整数值的行
merged_json_willetts.dropna(subset=['f.eid'], inplace=True)
#df_json_willetts_use = pd.read_csv(r'E:\QMY\data\df_json_willetts_use.csv')
medi = pd.merge(merged_result_0409,merged_json_willetts,on='f.eid')
medi = medi.query('wearTime_overall>6.5')
medi.to_csv(r'E:\QMY\result\Data0424.csv', index=False,sep=',')

#级数数据
ft_multi_b =  pd.read_csv(r'E:\QMY\data\ft_multi_b.csv')
ft_one_b_inside = pd.read_csv(r'E:\QMY\data\ft_one_b_inside.csv')
ft_one_b_outside = pd.read_csv(r'E:\QMY\data\ft_one_b_outside.csv')

ft_multi_b[['c1','c2','c3','c4','c5','c6',
         'c7','c8','c9','c10']] = ft_multi_b[['c1','c2','c3',
                              'c4','c5','c6','c7','c8','c9','c10']].abs() #将各个c取绝对值
ft_one_b_inside[['c1','c2','c3','c4','c5','c6',
           'c7','c8','c9','c10']] = ft_one_b_inside[['c1','c2','c3',
                                'c4','c5','c6','c7','c8','c9','c10']].abs()
ft_one_b_outside[['c1','c2','c3','c4','c5','c6',
           'c7','c8','c9','c10']] = ft_one_b_outside[['c1','c2','c3',
                                'c4','c5','c6','c7','c8','c9','c10']].abs()
                                                
multi_b_0410 = pd.merge(ft_multi_b,terminate,on='f.eid')
b_inside_0410 = pd.merge(ft_one_b_inside,terminate,on='f.eid')
b_outside_0410 = pd.merge(ft_one_b_outside,terminate,on='f.eid')
multi_b_0410.to_csv(r'E:\QMY\result\multi_b_0410.csv',index=False,sep=',')
b_inside_0410.to_csv(r'E:\QMY\result\b_inside_0410.csv',index=False,sep=',')
b_outside_0410.to_csv(r'E:\QMY\result\b_outside_0410.csv',index=False,sep=',')





'''4.11 16464有步数人的生存分析'''
var_all_0411 =  pd.read_csv(r'E:\QMY\result\var_all_0411.csv')
step_0410 = pd.read_csv(r'E:\QMY\result\Data0411.csv')
merge_0411 = pd.merge(var_all_0411,step_0410,on='f.eid')
merge_0411.columns[:50]
merge_0411.columns[50:100]
merge_0411.columns[100:150]
merge_0411.columns[150:200]
merge_0411.columns[200:250]
merge_0411.columns[250:300]

def quantile(table,column): #分位数
    print(table[column].min(),
      table[column].quantile(0.25),
      table[column].median(),
      table[column].quantile(0.75),
      table[column].max())

'步数'
quantile(merge_0411,'step_sum_avg')
plt.hist(merge_0411['step_sum_avg'],50)
merge_0411['Q1'] = merge_0411['Q2'] = merge_0411['Q3'] = merge_0411['Q4'] = merge_0411['other'] = np.zeros(len(merge_0411))
ix1 = merge_0411.query('step_sum_avg <5719.24999975').index
ix2 = merge_0411.query('5719.24999975<= step_sum_avg<=6997.9285715').index
ix3 = merge_0411.query('6997.9285715< step_sum_avg <8413.535714').index
ix4 = merge_0411.query('8413.535714<= step_sum_avg').index
merge_0411['Q1'][ix1] = 1
merge_0411['Q2'][ix2] = 1
merge_0411['Q3'][ix3] = 1
merge_0411['Q4'][ix4] = 1
merge_0411.to_csv(r'E:\QMY\result\Data0411.csv', index=False,sep=',')
kmf_steps_4(merge_0411,0.88) #分组生存函数

'30 分钟峰值步频'
plt.hist(merge_0411['max_30_minute_avg_avg'],50)
quantile(merge_0411,'max_30_minute_avg_avg')
merge_0411['Q1_max30'] = merge_0411['Q2_max30'] = merge_0411['Q3_max30'] = merge_0411['Q4_max30']  = np.zeros(len(merge_0411))
ix1 = merge_0411.query('max_30_minute_avg_avg <40.160714285').index
ix2 = merge_0411.query('40.160714285<= max_30_minute_avg_avg <=48.47857143').index
ix3 = merge_0411.query('48.47857143< max_30_minute_avg_avg <58.37142857').index
ix4 = merge_0411.query('58.37142857<= max_30_minute_avg_avg').index
merge_0411['Q1_max30'][ix1] = 1
merge_0411['Q2_max30'][ix2] = 1
merge_0411['Q3_max30'][ix3] = 1
merge_0411['Q4_max30'][ix4] = 1


kmf = KaplanMeierFitter()
rcParams['axes.titlepad'] = 30 
plt.rcParams.update({'font.size':25})
plt.figure(figsize=(12,8), dpi=80) 
ix1 = merge_0411.query('Q1_max30 ==1').index
ix2 = merge_0411.query('Q2_max30 ==1').index
ix3 = merge_0411.query('Q3_max30 ==1').index
ix4 = merge_0411.query('Q4_max30 ==1').index
kmf.fit(merge_0411['Duration_ukb'][ix1], merge_0411['termi'][ix1], label='Max_cadence30_Q1')
ax = kmf.plot()
kmf.fit(merge_0411['Duration_ukb'][ix2], merge_0411['termi'][ix2], label='Max_cadence30_Q2')
ax = kmf.plot(ax=ax)
kmf.fit(merge_0411['Duration_ukb'][ix3], merge_0411['termi'][ix3], label='Max_cadence30_Q3')
ax = kmf.plot(ax=ax)
kmf.fit(merge_0411['Duration_ukb'][ix4], merge_0411['termi'][ix4], label='Max_cadence30_Q4')
ax = kmf.plot(ax=ax)
plt.xlim(0,15)
plt.ylim(0.88,1.0)
#plt.text(0.2,0.89,'P_value<0.005')
plt.xlabel('随访时长（年）',rotation=0,fontsize=25)
plt.ylabel('Survival Probability',rotation=90,fontsize=25)
plt.title('生存函数',fontsize=32)

def logrank():
    merge_0411['class'] = np.zeros(len(merge_0411))
    merge_0411['class'][ix1] = 'Q1'
    merge_0411['class'][ix2] = 'Q2'
    merge_0411['class'][ix3] = 'Q3'
    merge_0411['class'][ix4] = 'Q4'
    results = pairwise_logrank_test(
        event_durations = merge_0411['Duration_ukb'],
        groups=merge_0411['class'],
        event_observed=merge_0411['termi']
   )
    summary = results.print_summary()
    print("p_value:{}".format(results.p_value))

#步数与30 分钟峰值步频的关系
plt.scatter(merge_0411['step_sum_avg'],merge_0411['max_30_minute_avg_avg'])


'max_1_minute_avg'
plt.hist(merge_0411['max_1_minute_avg'],50)

'low_cadence'
plt.hist(merge_0411['low_cadence_step_sum_avg'],50)
plt.hist(merge_0411['low_cadence_walk_minute_avg'],50)

plt.figure(figsize=(12,8), dpi=80) 
plt.hist(merge_0411['low_cadence_step_sum_avg'],50)
plt.hist(merge_0411['mid_cadence_step_sum_avg'],50)
plt.hist(merge_0411['high_cadence_step_sum_avg'],50)
plt.hist(merge_0411['bicycling_overall_avg'],50)
plt.show()
print(merge_0411['step_sum_avg'].mean(),
      merge_0411['low_cadence_step_sum_avg'].mean(),
      merge_0411['mid_cadence_step_sum_avg'].mean(),
      merge_0411['high_cadence_step_sum_avg'].mean())

quantile(merge_0411,'low_cadence_step_sum_avg')
merge_0411['Q1_lowcad'] = merge_0411['Q2_lowcad'] = merge_0411['Q3_lowcad'] = merge_0411['Q4_lowcad']  = np.zeros(len(merge_0411))
ix1 = merge_0411.query('low_cadence_step_sum_avg <2893.285714').index
ix2 = merge_0411.query('2893.285714<= low_cadence_step_sum_avg <=3444.7857145').index
ix3 = merge_0411.query('3444.7857145< low_cadence_step_sum_avg <4015.142857').index
ix4 = merge_0411.query('4015.142857<= low_cadence_step_sum_avg').index
merge_0411['Q1_lowcad'][ix1] = 1
merge_0411['Q2_lowcad'][ix2] = 1
merge_0411['Q3_lowcad'][ix3] = 1
merge_0411['Q4_lowcad'][ix4] = 1

kmf = KaplanMeierFitter()
rcParams['axes.titlepad'] = 30 
plt.rcParams.update({'font.size':25})
plt.figure(figsize=(12,8), dpi=80) 
ix1 = merge_0411.query('Q1_lowcad ==1').index
ix2 = merge_0411.query('Q2_lowcad ==1').index
ix3 = merge_0411.query('Q3_lowcad ==1').index
ix4 = merge_0411.query('Q4_lowcad ==1').index
kmf.fit(merge_0411['Duration_ukb'][ix1], merge_0411['termi'][ix1], label='lowcad_Q1')
ax = kmf.plot()
kmf.fit(merge_0411['Duration_ukb'][ix2], merge_0411['termi'][ix2], label='lowcad_Q2')
ax = kmf.plot(ax=ax)
kmf.fit(merge_0411['Duration_ukb'][ix3], merge_0411['termi'][ix3], label='lowcad_Q3')
ax = kmf.plot(ax=ax)
kmf.fit(merge_0411['Duration_ukb'][ix4], merge_0411['termi'][ix4], label='lowcad_Q4')
ax = kmf.plot(ax=ax)
plt.xlim(0,15)
plt.ylim(0.88,1.0)
#plt.text(0.2,0.89,'P_value<0.005')
plt.xlabel('随访时长（年）',rotation=0,fontsize=25)
plt.ylabel('Survival Probability',rotation=90,fontsize=25)
plt.title('生存函数',fontsize=32)
logrank()



'high_cadence'
quantile(merge_0411,'high_cadence_step_sum_avg')
merge_0411['Q1_highcad'] = merge_0411['Q2_highcad'] = merge_0411['Q3_highcad'] = merge_0411['Q4_highcad']  = np.zeros(len(merge_0411))
ix1 = merge_0411.query('high_cadence_step_sum_avg <677.85714285').index
ix2 = merge_0411.query('677.85714285<= high_cadence_step_sum_avg <=1302.142857').index
ix3 = merge_0411.query('1302.142857< high_cadence_step_sum_avg <2242.7142855').index
ix4 = merge_0411.query('2242.7142855<= high_cadence_step_sum_avg').index
merge_0411['Q1_highcad'][ix1] = 1
merge_0411['Q2_highcad'][ix2] = 1
merge_0411['Q3_highcad'][ix3] = 1
merge_0411['Q4_highcad'][ix4] = 1

kmf = KaplanMeierFitter()
rcParams['axes.titlepad'] = 30 
plt.rcParams.update({'font.size':25})
plt.figure(figsize=(12,8), dpi=80) 
ix1 = merge_0411.query('Q1_highcad ==1').index
ix2 = merge_0411.query('Q2_highcad ==1').index
ix3 = merge_0411.query('Q3_highcad ==1').index
ix4 = merge_0411.query('Q4_highcad ==1').index
kmf.fit(merge_0411['Duration_ukb'][ix1], merge_0411['termi'][ix1], label='highcad_Q1')
ax = kmf.plot()
kmf.fit(merge_0411['Duration_ukb'][ix2], merge_0411['termi'][ix2], label='highcad_Q2')
ax = kmf.plot(ax=ax)
kmf.fit(merge_0411['Duration_ukb'][ix3], merge_0411['termi'][ix3], label='highcad_Q3')
ax = kmf.plot(ax=ax)
kmf.fit(merge_0411['Duration_ukb'][ix4], merge_0411['termi'][ix4], label='highcad_Q4')
ax = kmf.plot(ax=ax)
plt.xlim(0,15)
plt.ylim(0.9,1.0)
plt.xlabel('随访时长（年）',rotation=0,fontsize=25)
plt.ylabel('Survival Probability',rotation=90,fontsize=25)
plt.title('生存函数',fontsize=32)


'mid_cadence'
quantile(merge_0411,'mid_cadence_step_sum_avg')
merge_0411['Q1_midcad'] = merge_0411['Q2_midcad'] = merge_0411['Q3_midcad'] = merge_0411['Q4_midcad']  = np.zeros(len(merge_0411))
ix1 = merge_0411.query('mid_cadence_step_sum_avg <1368.142857').index
ix2 = merge_0411.query('1368.142857<= mid_cadence_step_sum_avg <=1898.428571').index
ix3 = merge_0411.query('1898.428571< mid_cadence_step_sum_avg <2528.6428575').index
ix4 = merge_0411.query('2528.6428575<= mid_cadence_step_sum_avg').index
merge_0411['Q1_midcad'][ix1] = 1
merge_0411['Q2_midcad'][ix2] = 1
merge_0411['Q3_midcad'][ix3] = 1
merge_0411['Q4_midcad'][ix4] = 1

kmf = KaplanMeierFitter()
rcParams['axes.titlepad'] = 30 
plt.rcParams.update({'font.size':25})
plt.figure(figsize=(12,8), dpi=80) 
ix1 = merge_0411.query('Q1_midcad ==1').index
ix2 = merge_0411.query('Q2_midcad ==1').index
ix3 = merge_0411.query('Q3_midcad ==1').index
ix4 = merge_0411.query('Q4_midcad ==1').index
kmf.fit(merge_0411['Duration_ukb'][ix1], merge_0411['termi'][ix1], label='midcad_Q1')
ax = kmf.plot()
kmf.fit(merge_0411['Duration_ukb'][ix2], merge_0411['termi'][ix2], label='midcad_Q2')
ax = kmf.plot(ax=ax)
kmf.fit(merge_0411['Duration_ukb'][ix3], merge_0411['termi'][ix3], label='midcad_Q3')
ax = kmf.plot(ax=ax)
kmf.fit(merge_0411['Duration_ukb'][ix4], merge_0411['termi'][ix4], label='midcad_Q4')
ax = kmf.plot(ax=ax)
plt.xlim(0,15)
plt.ylim(0.86,1.0)
plt.xlabel('随访时长（年）',rotation=0,fontsize=25)
plt.ylabel('Survival Probability',rotation=90,fontsize=25)
plt.title('生存函数',fontsize=32)
logrank()


#0504
plt.hist(data2['mean_minimum'],50,range=(0,100))
plt.boxplot(data2['mean_minimum'])

data = data[['f.eid', 
       'minimum_avg', 'amp_avg', 'alpha_avg', 'beta_avg', 'acrotime_avg',
       'F_pseudo_avg', 'UpMesor_avg', 'DownMesor_avg', 'MESOR_avg']]
datamerge =pd.merge(data,data2,on='f.eid')
datamerge.to_csv(r'E:\QMY\result\Datamerge0504.csv', index=False,sep=',')









#cox回归
merge_cox = merge4[['avg_step_count_y','ethnic_x','Duration_ukb','termi',
                    'sex','age','smoke','TDI','alcohol','metabolic_syd']]
cph = CoxPHFitter()
# 拟合 Cox 模型
cph.fit(merge_cox, duration_col='Duration_ukb', event_col='termi')
# 打印 Cox 模型的系数
np.set_printoptions(suppress=True)
cph.summary


import subprocess
import os
subprocess.call(["Rscript", 'E:/QMY/Code0322.R'])

with open('E:/QMY/Code0322.R', "r") as f:
    r_script_content = f.read()

 


