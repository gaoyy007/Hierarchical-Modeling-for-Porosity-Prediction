%插值


%% b porosity
[x,y]=meshgrid(X,Y);
figure;
mesh(x,y,meanb); %无插值曲面

figure;
xi=100:10:2000;
yi=150:10:400;
[xi,yi]=meshgrid(xi,yi);
zi=interp2(x,y,meanb,xi,yi,'spline');
mesh(xi,yi,zi); %有插值曲面




figure;%热力网格图
heatmap(X,Y,meana,'colormap',jet);
saveas(gcf,'heatmap_meana.fig')
figure;
heatmap(X,Y,meanb,'colormap',jet);
saveas(gcf,'heatmap_meanb.fig')
figure;
heatmap(X,Y,stda,'colormap',jet);
saveas(gcf,'heatmap_stda.fig')
figure;
heatmap(X,Y,stdb,'colormap',jet);
saveas(gcf,'heatmap_stdb.fig')








%% a
[x,y]=meshgrid(X,Y);
figure;
mesh(x,y,meana);

figure;
xi=100:10:2000;
yi=150:10:400;
[xi,yi]=meshgrid(xi,yi);
zi=interp2(x,y,meana,xi,yi,'spline');
mesh(xi,yi,zi);


%cz=interp2(X,Y,meana,cx,cy,'method');