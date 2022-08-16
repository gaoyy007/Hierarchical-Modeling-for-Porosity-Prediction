% error of meanb
coe_meanb=[0.147,2.09e-4,-1.36e-3,4.24e-8,2.88e-6,-7.61e-7];
fitted_meanb=coe_meanb(1)+coe_meanb(2)*x+coe_meanb(3)*y+coe_meanb(4)*x.^2+coe_meanb(5)*y.^2+coe_meanb(6)*x.*y;
error_meanb=meanb-fitted_meanb;


% error 光滑热力图；
xi=100:10:2000;
yi=150:10:400;
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

%绘制热图； 还有网格线
figure;
xi=100:10:2000;
yi=150:10:400;
[xi,yi]=meshgrid(xi,yi);
zi=interp2(x,y,error_meanb,xi,yi,'spline');
xi=100:10:2000;
yi=150:10:400;
heatmap(xi,yi,zi,'colormap',jet);

%绘制三维
figure;
xi=100:10:2000;
yi=150:10:400;
[xi,yi]=meshgrid(xi,yi);
zi=interp2(x,y,error_meanb,xi,yi,'spline');
xi=100:10:2000;
yi=150:10:400;
mesh(xi,yi,zi);
shading interp 

%meshc 样例
x1=43:237;x2=4.3:23.7;
[x1,x2]=meshgrid(x1,x2);
Y=(-0.1146+0.1473*x1+0.6493*x2-0.0052*x1.*x2+0.0001*x1.^2+0.0033*x2.^2)/9;
mesh(x1,x2,Y);
shading interp 
figure;
[c h]=contour(x1,x2,Y);


%shading interp  样例
subplot(3,1,1)
sphere(16)
axis square
shading flat
title('Flat Shading')

subplot(3,1,2)
sphere(16)
axis square
shading faceted
title('Faceted Shading')

subplot(3,1,3)
sphere(16)
axis square
shading interp
title('Interpolated Shading')
