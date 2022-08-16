# -*- coding: utf-8 -*-
"""
Created on Mon Jul 26 14:34:00 2021

@author: Yuanyuan
"""


import os 
from PIL import Image
import numpy as np
import math
import csv
import matplotlib.pyplot as plt
import porespy as ps
from scipy.optimize import curve_fit
import random 
import time
import copy

def tanh(x,a,b):    
    return 1-b*(np.exp(a*x) - np.exp(-a*x)) / (np.exp(a*x) + np.exp(-a*x)) 
#定义energy
def E(regions):   
    tpcf_regions = np.array(ps.metrics.two_point_correlation_fft(im=regions))
    delta=tpcf_regions[1,:]-tpcf_criterion
    E=np.sum(delta**2)
    return E

threshold = 100 
table = []
for i in range(256):
    if i < threshold:
        table.append(0)
    else:
        table.append(1)
path=r'C:\Users\Yuanyuan\Desktop\HM\20210718\reconstruction'
os.chdir(path)

sample_for_recon = ["1.jpg","2.jpg","3.jpg","4.jpg","5.jpg","6.jpg","7.jpg","8.jpg"] 
for name in sample_for_recon:
    img = Image.open(name)  #
    img = img.convert('L')
    bi_image = img.point(table, '1')
    
    tpcf = ps.metrics.two_point_correlation_fft(im=bi_image )
    tpcf = np.array(tpcf) #tpcf for one bi_pic
    popt, pcov = curve_fit(tanh, tpcf[0,:], tpcf[-1,:],p0=[0, 1])
    bi_image=np.array(bi_image)
    ratio=np.sum(bi_image==True)/(200*200)
    
    
    pixel=40000
    pixel_edge=200
    swap_num=1
    regions_init = np.zeros((pixel,1),bool)
    init1=np.random.permutation(pixel)
    regions_init[init1[:np.int(ratio*pixel)]]=True
    regions_init=np.reshape(regions_init, (pixel_edge,pixel_edge))
    regions_initimg=Image.fromarray(regions_init)
    regions_initimg
    tpcf_init = np.array(ps.metrics.two_point_correlation_fft(im=regions_initimg))
    tpcf_criterion=tanh(tpcf_init[0,:],popt[0],popt[1])#popt 是预测的tpcf的参数
    Energy_init=E(regions_init)
    
    strat=time.time()
    swap_num=1
    Energy=[]
    DeltaEnergy=[]
    DeltaE_temp=[]
    Fail_num=[]
    
    
    DeltaE_pre=-1
    DeltaEnergy.append(DeltaE_pre)
    regions_pre=regions_init
    E_pre=Energy_init
    Energy.append(E_pre)
    iter_re=0
    iter_com=0#更小
    #失败的次数
    T=1e-4 #接受概率-3 比 -2 接受概率更小
    alpha=0.95
    var = 1
    iter_fail=0
    p_instant=0.3
    
    while(E_pre>1e-6): 
        regions_new=copy.deepcopy(regions_pre)
        x_pre,y_pre=np.where(regions_pre==True)
        #print(E(regions_pre))
        order=random.sample(list(range(len(x_pre))),swap_num)
        regions_new[x_pre[order],y_pre[order]]=False
        #print(E(regions_pre))
        #print(E(regions_new))
        
        
        x_pre,y_pre=np.where(regions_pre==False)
        #print(E(regions_pre))
        #print(E(regions_new))
        order=random.sample(list(range(len(x_pre))),swap_num)
        regions_new[x_pre[order],y_pre[order]]=True
        #print(E(regions_pre))
        #print(E(regions_new))
        
    
        E_new=E(regions_new)
        DeltaE_new=E_new-E_pre
        iter_re+=1
        DeltaE_temp.append(DeltaE_new)
        #print(DeltaE_new)
        
       
        if (DeltaE_new<=0 or random.random()<np.exp(-DeltaE_new/T)):
            regions_pre=copy.deepcopy(regions_new)
            E_pre=copy.deepcopy(E_new)
            DeltaE_pre=copy.deepcopy(DeltaE_new)
            Energy.append(E_pre)
            DeltaEnergy.append(DeltaE_pre)
            T=alpha*T
            iter_com+=1
            print(iter_com)
            Fail_num.append(iter_fail)
            iter_fail=0
        else:
            iter_fail+=1
            
            
            
                
    ####show the recon       
    regions_initimg=Image.fromarray(regions_pre)
    regions_initimg
      
    np.sum(DeltaEnergy)
    end=time.time()
    

