# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 21:45:37 2020

@author: GYY
MCMC_initial计算初始参数

"""


import os 
#############################sigmaE

def tanh(x,a,b):    
    return 1-np.multiply(b, np.divide(np.exp(np.multiply(a,x)) - np.exp(np.multiply(-a,x)),np.exp(np.multiply(a,x)) + np.exp(np.multiply(-a,x))))

error_e=[]
path=r"C:\Users\Yuanyuan\Desktop\HM\20210204\case matlab\new case"
os.chdir(path)
docnum=["01","02","03","04","05","06","10","11","12","13","14","15","16","17","18","19","20","21","22","23","25","26","27","28","29","31","32","33","34","35","36",
"37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","59","60"]
for i in docnum:

    Y_temp=np.loadtxt("tpcf_"+ i +".csv", delimiter=',')
    para_temp=np.loadtxt("para_"+ i +".csv", delimiter=',')
    for j in np.arange(16):
        error=tanh(Y_temp[0,:],para_temp[j,0],para_temp[j,1])-Y_temp[j+1,:]
        error_e.append( np.mean(error**2))
sigma_e=np.mean(error_e)      
    

import os 
from PIL import Image
import numpy as np
import math
import csv
import matplotlib.pyplot as plt
import porespy as ps
from scipy.optimize import curve_fit
import numpy
#from pyearth import Earth
from sklearn.datasets import make_friedman2
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import DotProduct, WhiteKernel
from scipy.optimize import curve_fit


#####################################2.有二次项拟合
path=r"C:\Users\Yuanyuan\Desktop\HM\20210204\case matlab\new case"
os.chdir(path)

data=np.loadtxt(r'C:\Users\Yuanyuan\Desktop\HM\20210204\case matlab\new case\h2.csv', delimiter=',')

xdata=data[:,0:2]
zdata=data[:,2]
def func(x,a,b,c,d,e,f):
    return a+b*x[0]+c*x[1]+d*x[0]*x[0]+e*x[1]*x[1]+f*x[0]*x[1]


popt, pcov = curve_fit(func, np.transpose(xdata), zdata,p0=[1,1,1,1,1,1])
###mse
z_hat = [func(i,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5]) for i in xdata]
z_hat_a=np.array(z_hat)
error=z_hat_a-zdata
plt.figure()
plt.scatter(np.arange(54 ),error)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.savefig('residual_a_beta.jpg')
np.savetxt('a_beta.csv', popt, delimiter = ',')  


xdata=data[:,0:2]
zdata=data[:,3]
def func(x,a,b,c,d,e,f):
    return a+b*x[0]+c*x[1]+d*x[0]*x[0]+e*x[1]*x[1]+f*x[0]*x[1]

popt, pcov = curve_fit(func, np.transpose(xdata), zdata,p0=[1,1,1,1,1,1])
###mse
z_hat = [func(i,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5]) for i in xdata]
z_hat_b=np.array(z_hat)
error=z_hat_b-zdata
plt.figure()
plt.scatter(np.arange(54 ),error)
plt.ylim(-0.2, 0.2)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.savefig('residual_b_beta.jpg')
np.savetxt('b_beta.csv', popt, delimiter = ',')  



xdata=data[:,0:2]
zdata=np.log(data[:,4])
def func(x,a,b,c,d,e,f):
    return a+b*x[0]+c*x[1]+d*x[0]*x[0]+e*x[1]*x[1]+f*x[0]*x[1]

popt, pcov = curve_fit(func, np.transpose(xdata), zdata,p0=[1,1,1,1,1,1])
###mse
z_hat = [func(i,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5]) for i in xdata]
z_hat_a_var=np.array(z_hat)
error=z_hat_a_var-zdata
sigmadelta_a_var=np.mean(error**2)
plt.figure()
plt.scatter(np.arange(54 ),error)
plt.ylim(-2, 2)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.savefig('residualz_a_var.jpg')
np.savetxt('a_var.csv', popt, delimiter = ',')  

xdata=data[:,0:2]
zdata=np.log(data[:,5])
def func(x,a,b,c,d,e,f):
    return a+b*x[0]+c*x[1]+d*x[0]*x[0]+e*x[1]*x[1]+f*x[0]*x[1]

popt, pcov = curve_fit(func, np.transpose(xdata), zdata,p0=[1,1,1,1,1,1])
###mse
z_hat = [func(i,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5]) for i in xdata]
z_hat_b_var=np.array(z_hat)
error=z_hat_b_var-zdata
sigmadelta_b_var=np.mean(error**2)
plt.figure()
plt.scatter(np.arange(54 ),error)
plt.ylim(-2, 2)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.savefig('residualz_b_var.jpg')
np.savetxt('b_var.csv', popt, delimiter = ',')  



##画出residual

##################################sigma与theta
sigma_a=[]
sigma_b=[]
ksi_a=np.empty([54,0])
ksi_b=np.empty([54,0])
docnum=["01","02","03","04","05","06","10","11","12","13","14","15","16","17","18","19","20","21","22","23","25","26","27","28","29","31","32","33","34","35","36",
"37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","59","60"]
i=0
for s in docnum:
    para=np.loadtxt("para_"+ s +".csv", delimiter=',')
    
    temp_a=z_hat_a[i]-para[:,0]
    temp_b=z_hat_b[i]-para[:,1]
    ksi_a=np.append(ksi_a,temp_a)
    ksi_b=np.append(ksi_b,temp_b)
    var_a=np.mean((z_hat_a[i]-para[:,0])**2)
    var_b=np.mean((z_hat_b[i]-para[:,1])**2)
    #sigma_a.append(var_a,axis=0)
    #sigma_b.append(var_b,axis=0)   
    sigma_a.append(var_a)
    sigma_b.append(var_b)   
    i=i+1
    
ksi_a=np.array(ksi_a)
ksi_b=np.array(ksi_b)
   
sigma_a=np.array(sigma_a)
sigma_b=np.array(sigma_b)

ksi_a=np.reshape(ksi_a,(54,-1))
ksi_b=np.reshape(ksi_a,(54,-1))

np.savetxt('ksi_a.csv', ksi_a, delimiter = ',')  
np.savetxt('ksi_b.csv', ksi_b, delimiter = ',')  

np.savetxt('sigma_a.csv', sigma_a, delimiter = ',')  
np.savetxt('sigma_b.csv', sigma_b, delimiter = ',')  
 


