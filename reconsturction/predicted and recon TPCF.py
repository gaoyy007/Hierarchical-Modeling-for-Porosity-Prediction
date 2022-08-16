# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 17:19:28 2021

@author: Yuanyuan
"""

import os
def tanh(x,a,b):    
    return 1-b*(np.exp(a*x) - np.exp(-a*x)) / (np.exp(a*x) + np.exp(-a*x)) 

os.chdir(r'C:\Users\Yuanyuan\Desktop\reconstruction\file')
img = Image.open('52.jpg')
threshold = 100 
table = []
for i in range(256):
    if i < threshold:
        table.append(1)
    else:
        table.append(0)
img = img.convert('L') #gray  
boole_img = img.point(table, '1')
tpcf_2c = ps.metrics.two_point_correlation_fft(im=boole_img)
tpcf_2c=np.array(tpcf_2c) #原始图


img = Image.open('bi_example1_recon.png')
threshold = 100 
table = []
for i in range(256):
    if i < threshold:
        table.append(1)
    else:
        table.append(0)
img = img.convert('L') #gray  
boole_img = img.point(table, '1')
tpcf_2cre = ps.metrics.two_point_correlation_fft(im=boole_img)
tpcf_2cre= np.array(tpcf_2cre) #重构图
popt, pcov = curve_fit(tanh, tpcf_2cre[0,:], tpcf_2cre[-1,:],p0=[0, 1])
y=tanh(tpcf_2cre[0,:],popt[0],popt[1])
  
    
    
    
##增加预测结果
x=tpcf_2c[0,:]




x_re=tpcf_2cre[0,:]
x_c=tpcf_2c[0,:]
plt.plot(x_c,tpcf_2c[1,:],'b',label='The TPCF of original image',linewidth=2)
plt.plot(x_re,y,'k',label='The predicted TPCF',linewidth=2)
plt.plot(x_re,tpcf_2cre[1,:],'r--',label='The TPCF of reconstructed image',linewidth=4)


plt.ylim(0, 1)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.legend(fontsize=13)




img = Image.open('57.jpg')
threshold = 100 
table = []
for i in range(256):
    if i < threshold:
        table.append(1)
    else:
        table.append(0)
img = img.convert('L') #gray  
boole_img = img.point(table, '1')
tpcf_2c = ps.metrics.two_point_correlation_fft(im=boole_img)
tpcf_2c=np.array(tpcf_2c) #原始图


img = Image.open('bi2_example2_recon.png')
threshold = 100 
table = []
for i in range(256):
    if i < threshold:
        table.append(1)
    else:
        table.append(0)
img = img.convert('L') #gray  
boole_img = img.point(table, '1')
tpcf_2cre = ps.metrics.two_point_correlation_fft(im=boole_img)
tpcf_2cre= np.array(tpcf_2cre) #重构图
popt, pcov = curve_fit(tanh, tpcf_2cre[0,:], tpcf_2cre[-1,:],p0=[0, 1])
y=tanh(tpcf_2cre[0,:],popt[0],popt[1])
  
    
    
    
##增加预测结果
x=tpcf_2c[0,:]




x_re=tpcf_2cre[0,:]
x_c=tpcf_2c[0,:]
plt.plot(x_c,tpcf_2c[1,:],'b',label='The TPCF of original image',linewidth=2)
plt.plot(x_re,y,'k',label='The predicted TPCF',linewidth=2)
plt.plot(x_re,tpcf_2cre[1,:],'r--',label='The TPCF of reconstructed image',linewidth=4)

plt.ylim(0, 1)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.legend(fontsize=13)






