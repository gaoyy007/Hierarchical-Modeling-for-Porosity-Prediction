# -*- coding: utf-8 -*-
"""
Created on Thu May 14 16:25:49 2020

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
def E(regions):   
    tpcf_regions = np.array(ps.metrics.two_point_correlation_fft(im=regions))
    delta=tpcf_regions[1,:]-tpcf_criterion
    E=np.sum(delta**2)
    return E

path=r'C:\Users\Yuanyuan\Desktop\recon'
os.chdir(path)
#tpcf_example=two_point_cor_fun('example.jpg')
img = Image.open('example.jpg')
threshold = 100 
table = []
for i in range(256):
    if i < threshold:
        table.append(1)
    else:
        table.append(0)
img = img.convert('L') #gray   
bi_image = img.point(table, '1')

name='bi_'+'example2.jpg'
bi_image.save(name) #save pic

tpcf = ps.metrics.two_point_correlation_fft(im=bi_image)
tpcf = np.array(tpcf) #tpcf for one bi_pic
popt, pcov = curve_fit(tanh, tpcf[0,:], tpcf[-1,:],p0=[0, 1])
bi_image=np.array(bi_image)
ratio=np.sum(bi_image==True)/(2208* 2752)


bi_image=bi_image+0

