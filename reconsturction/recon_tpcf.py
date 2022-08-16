# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 08:30:53 2020

@author: Yuanyuan
原始的tpcf 重建的tpcf
"""

os.chdir(r'C:\Users\Yuanyuan\Desktop\recon')
img = Image.open('example2_com.jpg')
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
tpcf_2c=np.array(tpcf_2c)


img = Image.open('bi_example2_recon.png')
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
tpcf_2cre= np.array(tpcf_2cre)


x_re=tpcf_2cre[0,:]
x_c=tpcf_2c[0,:]
plt.plot(x_re,tpcf_2cre[1,:],'r--',label='The TPCF of reference',linewidth=3)
plt.plot(x_c,tpcf_2c[1,:],label='The TPCF of reconstruction image',linewidth=3)
plt.ylim(0, 1)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.legend(fontsize=13)
plt.savefig('TPCF_2') 


img = Image.open('bi_example1_compress.jpg')
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
tpcf_2c=np.array(tpcf_2c)


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
tpcf_2cre= np.array(tpcf_2cre)


x_re=tpcf_2cre[0,:]
x_c=tpcf_2c[0,:]
plt.figure()
plt.plot(x_re,tpcf_2cre[1,:],'r--',label='The TPCF of reference',linewidth=3)
plt.plot(x_c,tpcf_2c[1,:],label='The TPCF of reconstruction image',linewidth=3)
plt.ylim(0, 1)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.legend(fontsize=13)
plt.savefig('TPCF_1') 


################################Energy作图 energy需要提前保存的
plt.plot(Energy,linewidth=2.5)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.savefig('Energy.jpg') 
