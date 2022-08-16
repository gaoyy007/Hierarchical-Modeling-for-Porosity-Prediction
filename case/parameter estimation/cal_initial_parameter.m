%% initial alpha &ksi &sigmaKSI
load('theta.mat')
load('a.mat')
load('b.mat')
m=54;
l=16;
n=2;
p=2;
q=6;
r=6;
alpha_initial=zeros(m,l,p);
alpha_initial(:,:,1)=a;
alpha_initial(:,:,2)=b;
ksi_initial=zeros(m,l,p);
beta1=[-3.733,-0.0005964,0.00441,-5.569e-08,-1.045e-05,1.245e-06];
beta2=[-2.061,0.003472,-0.02833,5.001e-07,6.071e-05,-9.857e-06];
gam1=[-8.641,-0.0005736,0.02684,-3.808e-07,-4.617e-05,3.2e-06];
gam2=[-3.365,0.0007443,-0.01262,7.584e-07,2.522e-05,-4.649e-06];
%Ô­Êý¾Ý

%initial value for Gibbs
gam(1,:)=[gam1,gam2];
beta(1,:)=[beta1,beta2];  
for i=1:m
    ksi_temp=phi(reshape(alpha_initial(i,:,:),2,16))-repmat(reshape(H_2(theta(i,:),p,q),p,p*q)*beta(1,:)',1,16);
    ksi_initial(i,:,:)=ksi_temp';
end

for i=1:m
    sigmaKSIinitial(i,1)=std(ksi_initial(i,:,1));
    sigmaKSIinitial(i,2)=std(ksi_initial(i,:,2));   
end


sigmaDelta_temp=zeros(2,54);
std_ab=zeros(2,54);
std_ab(1,:)=log(stda(:));
std_ab(2,:)=log(stdb(:));
sigmaDelta=zeros(1,2);

for i=1:m
    sigmaDelta_temp(:,i)=std_ab(:,i)-reshape(H_3(theta(i,:),p,r),p,p*r)*gam(1,:)';
end
sigmaDelta(1,1)=std(sigmaDelta_temp(1,:),0);
sigmaDelta(1,2)=std(sigmaDelta_temp(2,:),0);


save('sigmaDeltainitial.mat','sigmaDelta') 

save('sigmaKSIinitial.mat','sigmaKSIinitial') 
save('ksi_initial.mat','ksi_initial') 
save('alpha_initial.mat','alpha_initial') 

