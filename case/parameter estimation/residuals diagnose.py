# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 16:14:00 2022

@author: Yuanyuan
"""

#TPCF
#fitted TPCF
#noise distibution



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



###########100
threshold = 100 
table = []
for i in range(256):
    if i < threshold:
        table.append(1)
    else:
        table.append(0)


path=r'C:\Users\Yuanyuan\Desktop\HM\20220217\case\30case\recon\tpcf\TPCF-FITTED'
os.chdir(path)
num = [14,29,35,61]
for i in num:
    img = Image.open(str(i)+'.jpg')  #压缩后200*200
    img = img.convert('L')
    bi_image = img.point(table, '1')
    bi_image
    bi_image.save('bi_'+str(i)+'_100.jpg')
    
    tpcf_100 = ps.metrics.two_point_correlation_fft(im=bi_image )
    tpcf_100 = np.array(tpcf_100) #tpcf for one bi_pic
    popt_100, pcov_100 = curve_fit(tanh, tpcf_100[0,:], tpcf_100[-1,:],p0=[0, 1])
    print (popt_100)
    
    Y= tanh(tpcf_100[0,:], popt_100[0], popt_100[1])
    error=tanh(tpcf_100[0,:], popt_100[0], popt_100[1])-tpcf_100[1,:]
    error_mean = np.mean(error)
    error_std = np.std(error)
    print(error_mean)
    print(error_std)
    
    
    plt.figure()
    plt.plot(tpcf_100[0,:],tpcf_100[1,:])
    plt.plot(tpcf_100[0,:],Y)    
    plt.ylim(0,1)
    plt.savefig(str(i)+'-TPCF.jpg')













