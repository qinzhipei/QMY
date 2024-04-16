# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 19:24:57 2023

@author: Administrator
"""

import rpy2.robjects as robjects
robjects.r['source'](r'E:\QMY\Code0322.R')

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import time

result_max30_1 = pd.read_csv(r'E:\QMY\result\HR_result_1_0414.csv')
result_max30_2 = pd.read_csv(r'E:\QMY\result\HR_result_2_0414.csv')
result_max30_3 = pd.read_csv(r'E:\QMY\result\HR_result_3_0414.csv')
result_max30_4 = pd.read_csv(r'E:\QMY\result\HR_result_4_0414.csv')
result_max30_5 = pd.read_csv(r'E:\QMY\result\HR_result_5_0414.csv')
max_30_minute_avg_0 = pd.read_csv(r'E:\QMY\result\max_30_minute_avg_avg0_0414.csv')
max_30_minute_avg_1 = pd.read_csv(r'E:\QMY\result\max_30_minute_avg_avg1_0414.csv')
max_30_minute_avg_2 = pd.read_csv(r'E:\QMY\result\max_30_minute_avg_avg2_0414.csv')
step_sum_0 = pd.read_csv(r'E:\QMY\result\step_sum_avg_0_0414.csv')
step_sum_1 = pd.read_csv(r'E:\QMY\result\step_sum_avg_1_0414.csv')
step_sum_2 = pd.read_csv(r'E:\QMY\result\step_sum_avg_2_0414.csv')


def Three(x,y,zlim,xlabel,ylabel):
    fig = plt.figure(figsize=(20, 10),dpi=100)             
    ax = fig.gca(fc='whitesmoke',projection='3d')
    #ax.plot_surface(x.iloc[i],y.iloc[m],result.iloc[i,m],rstride=1, cstride=2, cmap=plt.cm.Blues)
    for i in range(len(x)):
        for m in range(len(y)):
            #x,y = np.meshgrid(x,y)
            point = ax.scatter3D(xs=x.iloc[i],ys=y.iloc[m],zs=result.iloc[i,m],  
            zdir='z',color="blue", s=60)
    ax.set(xlabel=xlabel,ylabel=ylabel,zlabel='Hazard Ratio')
    ax.set_zlim(0,zlim)
    fig.colorbar(point, ax = ax, shrink = 0.6, aspect = 5)
    plt.show()

def Three3(x0,y0,z0,x1,y1,z1,x2,y2,z2,zmin,zlim,xlabel,ylabel,elev,azim):
    fig = plt.figure(figsize=(20, 10),dpi=100)             
    ax = plt.axes(projection='3d')
    xx0 = np.array(x0['x'])
    yy0 = np.array(y0['x'])
    xx1 = np.array(x1['x'])
    yy1 = np.array(y1['x'])
    xx2 = np.array(x2['x'])
    yy2 = np.array(y2['x'])
    X0, Y0 = np.meshgrid(xx0, yy0)
    X1, Y1 = np.meshgrid(xx1, yy1)
    X2, Y2 = np.meshgrid(xx2, yy2)
    zt0 =pd.DataFrame(z0.values.T)
    zt1 =pd.DataFrame(z1.values.T)
    zt2 =pd.DataFrame(z1.values.T)
    ax1 = ax.plot_surface(X0,Y0,np.array(zt0)*0.98,rstride=1, cstride=2, cmap=plt.cm.Blues)
    ax2 = ax.plot_surface(X1,Y1,np.array(zt1),rstride=1, cstride=2, cmap=plt.cm.Reds)
    ax3 = ax.plot_surface(X1,Y1,np.array(zt2)*1.05,rstride=1, cstride=2, cmap=plt.cm.Greens)
    ax.set_xlabel(xlabel,size=25, labelpad=15)
    ax.set_ylabel(ylabel,size=25, labelpad=15)
    ax.set_zlabel("Hazard Ratio",size=25, labelpad=10)
    ax.tick_params(labelsize=12)
                  #=xlabel,ylabel=ylabel,zlabel='Hazard Ratio')
    ax.set_zlim(zmin,zlim)
    ax.view_init(elev=elev,azim=azim) #仰角、方位角
    fig.colorbar(ax1, shrink=0.5, pad=-0.05,aspect=8)
    fig.colorbar(ax2, shrink=0.5, pad=-0.05,aspect=8)
    fig.colorbar(ax3, shrink=0.5, pad=0.1,aspect=8)
    ax1._facecolors2d = ax1._facecolor3d 
    ax2._facecolors2d = ax2._facecolor3d
    ax1._edgecolors2d = ax1._edgecolor3d 
    ax2._edgecolors2d = ax2._edgecolor3d
    ax3._edgecolors2d = ax3._edgecolor3d 
    color = ['cornflowerblue', 'r','limegreen']
    legend_lines = [mpl.lines.Line2D([0], [0], linestyle="none", marker='o', c=color[y]) for y in [0,1,2]]
    legend_labels = ['Age <50','Age 50-60','Age >60']
    ax.legend(legend_lines, legend_labels, numpoints=1,fontsize=20,loc='upper right')
    plt.show()

def Three2(x0,y0,z0,x1,y1,z1,zmin,zlim,xlabel,ylabel,elev,azim):
    fig = plt.figure(figsize=(20, 10),dpi=100)             
    ax = plt.axes(projection='3d')
    xx0 = np.array(x0['x'])
    yy0 = np.array(y0['x'])
    xx1 = np.array(x1['x'])
    yy1 = np.array(y1['x'])
    X0, Y0 = np.meshgrid(xx0, yy0)
    X1, Y1 = np.meshgrid(xx1, yy1)
    zt0 =pd.DataFrame(z0.values.T)
    zt1 =pd.DataFrame(z1.values.T)
    ax1 = ax.plot_surface(X0,Y0,np.array(zt0),rstride=1, cstride=2, cmap=plt.cm.Blues,label='metabolic_syd == 1')
    ax2 = ax.plot_surface(X1,Y1,np.array(zt1),rstride=1, cstride=2, cmap=plt.cm.Reds,label='metabolic_syd == 0')
    ax.set_xlabel(xlabel,size=25, labelpad=15)
    ax.set_ylabel(ylabel,size=25, labelpad=15)
    ax.set_zlabel("Hazard Ratio",size=25, labelpad=10)
    ax.tick_params(labelsize=12)
                  #=xlabel,ylabel=ylabel,zlabel='Hazard Ratio')
    ax.set_zlim(zmin,zlim)
    ax.view_init(elev=elev,azim=azim) #仰角、方位角
    fig.colorbar(ax1, shrink=0.5, pad=-0.05,aspect=8)
    fig.colorbar(ax2, shrink=0.5, pad=0.1,aspect=8)
    ax1._facecolors2d = ax1._facecolor3d 
    ax2._facecolors2d = ax2._facecolor3d
    ax1._edgecolors2d = ax1._edgecolor3d 
    ax2._edgecolors2d = ax2._edgecolor3d
    color = ['cornflowerblue', 'r']
    legend_lines = [mpl.lines.Line2D([0], [0], linestyle="none", marker='o', c=color[y]) for y in [0,1]]
    legend_labels = ['not excessive drinking','excessive drinking']
    ax.legend(legend_lines, legend_labels, numpoints=1,fontsize=20,loc='upper right')
    plt.show()


Three2(max_30_minute_avg_0,step_sum_0,result_max30_1,
         max_30_minute_avg_1,step_sum_1,result_max30_2,
         0.3,2.2,'max_30_minute_avg','avg_step_count',20,320) #代谢病

for i in range(0,360,20):
    Three3(max_30_minute_avg_0,step_sum_0,result_max30_3,
         max_30_minute_avg_1,step_sum_1,result_max30_4,
         max_30_minute_avg_2,step_sum_2,result_max30_5,
         0.3,2.1,'max_30_minute_avg','avg_step_count',30,i) #吸烟
    time.sleep(0.5)

         