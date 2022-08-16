a=[200:200:1200];
b=[0.1:0.1:0.9];

i=6;
j=1;

load('freq_chosen')
filename = [ 'attenuation'  num2str(a(i)) '_' num2str(b(j)) '.mat'  ] ;
load(filename)

plot(freq_chosen,Atten,'b')
axis([2,2.5,0.5,1.2])
hold on;


