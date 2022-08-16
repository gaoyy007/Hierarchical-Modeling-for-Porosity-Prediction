close all
clear
m=6;
n=9;
Atten1=zeros(54,20,32);
index=0;
for i=1:m
%     for j=1:n
j=9;
        index=index+1;
        filename=['.\matrix data\','simulationmatrix',num2str(i*200),'_','0.',num2str(j),'.mat'];
        load(filename);
        for k=1:20
            [Atten1(index,k,:),freq]=calculateAtten(tempcell(:,k),0.001,8);
        end
        figure
        plot(freq,reshape(Atten1(index,:,:),20,32));
        axis([2,3,0,1])
%     end
end
        
        




