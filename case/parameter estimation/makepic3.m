%% 三维散点

clear 
close all
clc
load('meanb.mat')
load('Y.mat')
load('X.mat')
Y=flip(Y);

log_meanb=log(meanb./(1-meanb));


[x,y]=meshgrid(X,Y);
figure;
scatter3(x(:),y(:),meanb(:),'filled')
hold on
%% 拟合曲面图 下方等高线
V=100:10:1500;
P=400:-10:150;
[V,P]=meshgrid(V,P);
beta2=[-2.061,0.003572,-0.02933,5.001e-07,5.071e-05,-9.857e-06];
top=ones(26,141);
variable=[top;V;P;V.^2;P.^2;P.*V];
funmeanbofcase=1.835-0.0001963*V-0.0489*P+3.812e-06*V.^2+9.912e-05*P.^2-1.511e-05*P.*V;


%funmeanbofcase=-2.061+0.003472*V-0.02833*P+5.001e-07*V.^2+6.071e-05*P.^2+-9.857e-06*P.*V;
funmeanbofcase=exp(funmeanbofcase)./(1+exp(funmeanbofcase))
mesh(V,P,funmeanbofcase);


% [1.845,-0.0001989,-0.04953,3.833e-06,9.894e-05,-1.451e-05]

%% 等高线图
[x,y]=meshgrid(X,Y);
% coe_meanb=[0.147,2.09e-4,-1.36e-3,4.24e-8,2.88e-6,-7.61e-7];
% fitted_meanb=coe_meanb(1)+coe_meanb(2)*x+coe_meanb(3)*y+coe_meanb(4)*x.^2+coe_meanb(5)*y.^2+coe_meanb(6)*x.*y;
V=[100,200,300,500,800,1000,1300,1500];
P=[400,350,300,250,200,150];
[V,P]=meshgrid(V,P);
% funmeanbofcase=1.845-0.0001989*V-0.04953*P+3.833e-06*V.^2+9.894e-05*P.^2-1.451e-05*P.*V;
funmeanbofcase=1.835-0.0001963*V-0.0489*P+3.812e-06*V.^2+9.912e-05*P.^2-1.511e-05*P.*V;
%funmeanbofcase=1.845-0.0001989*V-0.04953*P+3.833e-06*V.^2+9.894e-05*P.^2-1.451e-05*P.*V;
funmeanbofcase=exp(funmeanbofcase)./(1+exp(funmeanbofcase))
xi=100:10:1500;
yi=400:-10:150;
[xi,yi]=meshgrid(xi,yi);
zi=interp2(x,y,funmeanbofcase,xi,yi,'spline');
xi=xi(:);
yi=yi(:);
zi=zi(:);
SI=scatteredInterpolant(xi,yi,zi);
xi=100:10:1500;
yi=400:-10:150;
[xi,yi]=meshgrid(xi,yi);
Ci = SI(xi,yi);
figure
imagesc(Ci);
% figure;
% xi=100:10:1500;
% yi=400:-10:150;
% [xi,yi]=meshgrid(xi,yi);
% zi=interp2(x,y,funmeanbofcase,xi,yi,'spline');
% xi=100:10:1500;
% yi=150:10:400;
% mesh(xi,yi,zi);
% shading interp 




%% 拟合误差图
[x,y]=meshgrid(X,Y);
coe_meanb=[0.147,2.09e-4,-1.36e-3,4.24e-8,2.88e-6,-7.61e-7];
fitted_meanb=coe_meanb(1)+coe_meanb(2)*x+coe_meanb(3)*y+coe_meanb(4)*x.^2+coe_meanb(5)*y.^2+coe_meanb(6)*x.*y;
xi=100:10:2000;
yi=400:-10:150;
[xi,yi]=meshgrid(xi,yi);
zi=interp2(x,y,error_meanb,xi,yi,'spline');
xi=xi(:);
yi=yi(:);
zi=zi(:);
SI=scatteredInterpolant(xi,yi,zi);
xi=100:10:2000;
yi=150:10:400;
[xi,yi]=meshgrid(xi,yi);
Ci = SI(xi,yi);
figure
imagesc(Ci);
figure;
xi=100:10:2000;
yi=400:-10:150;
[xi,yi]=meshgrid(xi,yi);
zi=interp2(x,y,error_meanb,xi,yi,'spline');
xi=100:10:2000;
yi=150:10:400;
mesh(xi,yi,zi);
shading interp 