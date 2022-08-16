%颗粒圆度 
clear;
close all; 
%% 
%读取源图像 
I = imread('example2.PNG'); 
figure;
imshow(I); 
%% 
%灰度化、取反 
h = rgb2gray(I); 
figure;
imshow(h);%灰度图像 
h = imcomplement(h);%取反 
figure;
imshow(h); 
%% 
%中值滤波、二值化 
h = medfilt2(h,[4,4]); 
bw = im2bw(h,graythresh(h)); 
%% 
%消除噪点 
se = strel('disk',2); 
bw = imclose(bw,se); 
figure;
imshow(bw); 
%% 
%填补闭合图形，填充色为白色 
bw  = imfill(bw,'holes'); 
%% 
%边界寻找 
[B,L] = bwboundaries(bw,'noholes');
% 为每个闭合图形设置颜色显示 
figure;
imshow(label2rgb(L, @jet, [.5 .5 .5])) 
hold on 
for k = 1:length(B)   
    boundary = B{k};   
    plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2) 
end
%% 
%计算面积 
stats = regionprops(L,'Area','Centroid');

threshold = 0.94;

% 循环处理每个边界，length(B)是闭合图形的个数,即检测到的陶粒对象个数 
for k = 1:length(B)
    % 获取边界坐标'   
    boundary = B{k};
    % 计算周长   
    delta_sq = diff(boundary).^2;   
    perimeter = sum(sqrt(sum(delta_sq,2)));
    % 对标记为K的对象获取面积
    area = stats(k).Area;
    % 圆度计算公式4*PI*A/P^2
    metric = 4*pi*area/perimeter^2;
    % 结果显示   
    metric_string = sprintf('%2.2f',metric);
    % 用一个黑色小圆圈标记圆度大于threshold = 0.94 的对象   
    if metric > threshold     
        centroid = stats(k).Centroid;     
        plot(centroid(1),centroid(2),'ko');   
    end %设置显示字体   
    text(boundary(1,2)-35,boundary(1,1)+13,metric_string,'Color','y','FontSize',14,'FontWeight','bold');

end

title(['圆度识别结果，越圆越接近1，']);