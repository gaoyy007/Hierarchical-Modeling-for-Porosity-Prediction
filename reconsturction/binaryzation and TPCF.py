# -*- coding: utf-8 -*-
"""
Created on Tue Aug  3 10:26:29 2021

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

import os 
from PIL import Image
import numpy as np
import math
import csv
import matplotlib.pyplot as plt
import porespy as ps
from scipy.optimize import curve_fit
from parfor import parfor





def cal_entropy(img):#img array    
    tmp = []
    for i in range(256):
        tmp.append(0)
    val = 0
    k = 0
    entropy = 0
    #img=bi_image
    for i in range(len(img)):    
        for j in range(len(img[i])):        
            val = img[i][j]        
            tmp[val] = float(tmp[val] + 1)        
            k =  float(k + 1)
    for i in range(len(tmp)): 
        tmp[i] = float(tmp[i] / k)
    for i in range(len(tmp)):  
        if(tmp[i] == 0): 
            entropy = entropy
        else: entropy = float(entropy - tmp[i] * (math.log(tmp[i]) / math.log(2.0)))
    #print( entropy)
    return entropy


def two_point_cor_fun(picname):
    img = Image.open(picname)
    threshold = 100 
    table = []
    for i in range(256):
        if i < threshold:
            table.append(1)
        else:
            table.append(0)
    img = img.convert('L') #gray  
    boole_img = img.point(table, '1')
    
    bi_image=np.array(boole_img)+0
    square_percentage=np.sum(bi_image==1)/(262*210)
    
    tpcf = ps.metrics.two_point_correlation_fft(im=boole_img)
    tpcf = np.array(tpcf) #tpcf for one bi_pic
    return tpcf, square_percentage

def tanh(x,a,b):    
    return 1-b*(np.exp(a*x) - np.exp(-a*x)) / (np.exp(a*x) + np.exp(-a*x)) 




#二值化
def binary(img,picname):
    Img = img.convert('L') #gray        
    #Img.size 
    
    threshold = [40,60,80,100,120,140,160,180]
    table = []
    for i in range(256):
        if i < threshold[i]:
            table.append(0)
        else:
            table.append(1)
    bi_image = Img.point(table, '1')
    name="bi_"+picname
    bi_image.save(name) #save pic
    bi_image=np.array(bi_image)+0
    #square_percentage=np.sum(bi_image==1)/(2752*2208)
    return bi_image

intergration_path=r'C:\Users\Yuanyuan\Desktop\HM\20210718\reconstruction\example25'
new_folder=r'C:\Users\Yuanyuan\Desktop\HM\20210718\reconstruction\bi_example25'
#os.makedirs(new_folder)
file_list=os.listdir(intergration_path)
os.chdir(new_folder)
'''
for imagename in file_list:
    image=Image.open(imagename)

    print(imagename)
    print(image)
    binary(image,imagename)
'''
#计算TPCF
tpcf_file=None  
tpcf=None
para=None 
square_percentage=[]
path_offile=r'C:\Users\Yuanyuan\Desktop\HM\20210718\reconstruction\bi_example25_compressed'
for picname in os.listdir(path_offile):#read a pic           
    tpcf, square_percentage_value=two_point_cor_fun(picname)
    square_percentage_value=square_percentage_value.tolist()
    square_percentage.append(square_percentage_value)
    popt, pcov = curve_fit(tanh, tpcf[0,:], tpcf[-1,:],p0=[0, 1])
    if para is None:
        para=popt
    else:
        para=np.vstack((para,popt))
    
    if tpcf_file is None:
        tpcf_file  = tpcf
    else:
        tpcf_file = np.vstack((tpcf_file, tpcf[-1,:]))   


os.chdir(r'C:\Users\Yuanyuan\Desktop\HM\20210718\reconstruction')
tpcf_re, square_percentage_value=two_point_cor_fun('re1.png')

    
x = tpcf_file[0,:]
for y in range(16) :
    plt.plot(tpcf_re[0,:],tpcf_re[1,:])
    plt.plot(x,tpcf_file[y+1,:] )
plt.ylim(0, 1)
plt.xlabel('Pixel')
plt.ylabel('The two point correlation function')


        
