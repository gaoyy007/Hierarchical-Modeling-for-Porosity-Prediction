function [Atten,freq]=calculateAtten(signal,t_res,dist)
num_new=length(signal);
pixel_x=t_res;
xnew=(1:1:num_new)*t_res;
ynew=signal;
N=1000;
dividePoint=floor(num_new/2-100);
Energy=0;
for i=1:(num_new-N)
    x_fft=xnew(i:i+N-1);
    y_fft=ynew(i:i+N-1);
    Fn=fft(y_fft);
    Energy(i)=0;
    for j=1:N
        Energy(i)=Energy(i)+(abs(Fn(j)))^2;
    end
end
nn=1:1:(num_new-N);
maxE1=max(Energy(1:dividePoint));
start_point1=find(Energy==maxE1);

maxE2=max(Energy(dividePoint:num_new-N));
start_point2=find(Energy==maxE2);

y_fft=ynew(start_point1:start_point1+N-1);
figure(1)
plot(xnew,ynew,'r');
hold on
plot(xnew(start_point1:start_point1+N-1),ynew(start_point1:start_point1+N-1),'b');
hold on
plot(xnew(start_point2:start_point2+N-1),ynew(start_point2:start_point2+N-1),'g');
y_fft=[y_fft;zeros(20*N,1)];
Fn1=fft(y_fft);
T=(21*N-1)*pixel_x;
freq=0:1/T:(21*N-1)/T;
y_fft=ynew(start_point2:start_point2+N-1);
y_fft=[y_fft;zeros(20*N,1)];
Fn2=fft(y_fft);


start1=round(2/(1/T));
start2=round(3/(1/T));
% freq_chosen=0:1/T:floor((31*N-1)/28)/T;
% Fn1_chosen=Fn1(1:floor((31*N-1)/28)+1);
% Fn2_chosen=Fn2(1:floor((31*N-1)/28)+1);

freq_chosen=start1/T:1/T:start2/T;
Fn1_chosen=Fn1(start1:start2);
Fn2_chosen=Fn2(start1:start2);


figure(2)
plot(freq_chosen,abs(Fn1_chosen),'g');
hold on
plot(freq_chosen,abs(Fn2_chosen),'r');
axis([2,3,0,max(abs(Fn1_chosen))])


% figure(1);
% hold on
%  needF=freq_chosen(nn:mm);
atten=log(abs(Fn1_chosen./Fn2_chosen))/dist*1000/2;
atten=atten*8.6858896/1000;
% plot(freq_chosen,atten,'r');
% axis([2,3,0.1,0.9])
Atten=atten;
freq=freq_chosen;



end

