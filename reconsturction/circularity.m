%����Բ�� 
clear;
close all; 
%% 
%��ȡԴͼ�� 
I = imread('example2.PNG'); 
figure;
imshow(I); 
%% 
%�ҶȻ���ȡ�� 
h = rgb2gray(I); 
figure;
imshow(h);%�Ҷ�ͼ�� 
h = imcomplement(h);%ȡ�� 
figure;
imshow(h); 
%% 
%��ֵ�˲�����ֵ�� 
h = medfilt2(h,[4,4]); 
bw = im2bw(h,graythresh(h)); 
%% 
%������� 
se = strel('disk',2); 
bw = imclose(bw,se); 
figure;
imshow(bw); 
%% 
%��պ�ͼ�Σ����ɫΪ��ɫ 
bw  = imfill(bw,'holes'); 
%% 
%�߽�Ѱ�� 
[B,L] = bwboundaries(bw,'noholes');
% Ϊÿ���պ�ͼ��������ɫ��ʾ 
figure;
imshow(label2rgb(L, @jet, [.5 .5 .5])) 
hold on 
for k = 1:length(B)   
    boundary = B{k};   
    plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2) 
end
%% 
%������� 
stats = regionprops(L,'Area','Centroid');

threshold = 0.94;

% ѭ������ÿ���߽磬length(B)�Ǳպ�ͼ�εĸ���,����⵽������������� 
for k = 1:length(B)
    % ��ȡ�߽�����'   
    boundary = B{k};
    % �����ܳ�   
    delta_sq = diff(boundary).^2;   
    perimeter = sum(sqrt(sum(delta_sq,2)));
    % �Ա��ΪK�Ķ����ȡ���
    area = stats(k).Area;
    % Բ�ȼ��㹫ʽ4*PI*A/P^2
    metric = 4*pi*area/perimeter^2;
    % �����ʾ   
    metric_string = sprintf('%2.2f',metric);
    % ��һ����ɫСԲȦ���Բ�ȴ���threshold = 0.94 �Ķ���   
    if metric > threshold     
        centroid = stats(k).Centroid;     
        plot(centroid(1),centroid(2),'ko');   
    end %������ʾ����   
    text(boundary(1,2)-35,boundary(1,1)+13,metric_string,'Color','y','FontSize',14,'FontWeight','bold');

end

title(['Բ��ʶ������ԽԲԽ�ӽ�1��']);