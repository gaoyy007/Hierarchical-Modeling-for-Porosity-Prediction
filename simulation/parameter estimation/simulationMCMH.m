%%%%%%%%%%%%%%%%%%%%%%%%simulation 
%%%%%%%%%%%%%%%%%%%%%%%para initial  sigma tune

clc
clear
close all

filename=['tpcf_02.csv'];
Y_temp=csvread(filename);
X=Y_temp(1,1:30)';  %列
X=X;
m=10;
l=5;
p=2;
n=30;
prop=1/4; % proportion of gibbs samplers to use
theta=zeros(m,2); 
load('theta.mat');
q=6;%q_s=3;
r=6;%r_s=3;
beta_s=[-3,-0.005,0.004,0,0,0,-2,0.003,-0.02,0,0,0];    
sigmaDelta_s=[0.1,0.1]; %!!
gam_s=[-8,-0.0005,0.02,0,0,0,-3,0.0007,-0.01,0,0,0];
sigmaE=0.0001;


for i=1:m
    MU_s=reshape(H_3(theta(i,:),p,r),p,p*r)*gam_s';
    sigmaKSI_s=exp(mvnrnd(MU_s,diag(sigmaDelta_s),m));
end
ksi_s=zeros(m,l,p);
for i=1:m
    for j=1:l
        ksi_s(i,:,:)=mvnrnd([0,0],diag(sigmaKSI_s(i,:)),l);
    end
end
alpha_s=zeros(m,l,p);
for i=1:m
    MU_i=reshape(H_2(theta(i,:),p,q),p,p*q)*beta_s';
    alpha_temp=MU_i'+reshape(ksi_s(i,:,:),l,p);
    alpha_temp=alpha_temp';
    alpha_temp=Iphi(alpha_temp); 
    alpha_s(i,:,:)=alpha_temp';  
end

Y=[];
for i=1:m
    for j=1:l
        
        Y_temp=H_1(X,n,alpha_s(i,j,:));
        noise=normrnd(0,sigmaE,n,1);
        Y_temp=Y_temp+noise;
        Y=[Y,Y_temp];        
    end   
end
sum(alpha_s(:,:,1));
sum(alpha_s(:,:,2));
Y=Y(:);
xlswrite('X.xls',X)


%%
%simulation
tic;
prop=1/4; % proportion of gibbs samplers to use


N=2000;
%---parameters------------%
beta=zeros(N,p*q);   %beta=(beta_1^t,...beta_p^T)')
sigmaE=zeros(N,1); %sigma2 for noise variance
ksi=zeros(N,m,l,p);     %\xi
sigmaKSI=zeros(N,m,p);  %sigma2 for ksi
gam=zeros(N,p*r);
sigmaDelta=zeros(N,p);  %sigma2 for level-2 variances
alpha=zeros(N,m,l,p);



%%%%%%%%%%%%%%%%%%%%%%%%2

beta(1,:)=[-3,-0.005,0.004,0,0,0,-2,0.003,-0.02,0,0,0];    

gam(1,:)=[-8,-0.0005,0.02,0,0,0,-3,0.0007,-0.01,0,0,0];
sigmaE(1)=0.0001;


sigmaDelta(1,:)=[0.001,0.001];
ksi(1,:,:,:)=ksi_s;
alpha(1,:,:,:)=alpha_s;
sigmaKSI(1,:,:)=sigmaKSI_s;  %sigma2 for ksi

H4=H_4(theta,m,r);
H4H4_inv=inv(H4'*H4);
H4H4H4=H4H4_inv*H4';

tic

for iter=2:N
    %if(mod(iter,1000)==0)
    iter;
    %end
    %step 1 sampling sigmaE (16*16*90,1)
    %KSI=XI(p,n,reshape(ksi(iter-1,:,:,:),m,l,p),m,l,H1_mat);
    %H1_mat=zeros(m*l*n,1)
    H1_mat=[];
    for i=1:m
        for j=1:l
            H1_mat=[H1_mat;H_1(X,n,alpha(iter-1,i,j,:))]; %(n,1)
        end
    end
    Y_H1_mat=Y-H1_mat;
    %beta_hat=HHH*Y_KSI;
    IG_a=m*n*l/2-1;
    IG_b=sum((Y_H1_mat).^2)/2;
    sigmaE(iter)=1/gamrnd(IG_a,1/IG_b); %generate inverse gamma
    %beta(iter,:)=mvnrnd(beta_hat,sigmaE(iter)*HH_inv);
    
    
    
    %step 2 sampling alpha Metropolis-Hastings sampling
    temp=1;
    for i=1:m
        MU_i=reshape(H_2(theta(i,:),p,q),p,p*q)*beta(iter-1,:)'; %先验
        
        sigmaKSI_i=diag(reshape(sigmaKSI(iter-1,i,:),2,1));
        for j=1:l
            alpha_t_1=reshape(alpha(iter-1,i,j,:),p,1);
            phi_alpha_new=phi(alpha_t_1)+mvnrnd(zeros(p,1),[0.1,0;0,0.1])'; %%%步长可调节
%             while alpha_new(1)<=0|| alpha_new(2)<=0||alpha_new(1)>=0.3||alpha_new(2)>=1
%                 alpha_new=alpha_t_1+mvnrnd(zeros(p,1),0.001*eye(p))'; %%%步长可调节             
%             end
            phi_alpha_new;              
           
            
            acceptR=1;
            acceptR=acceptR*(sigmaE(iter))^(-n/2)*exp(sum(  (Y((temp-1)*n+1:temp*n)-H_1(X,n,Iphi(phi_alpha_new))   ).^2      )/(-2*sigmaE(iter)));
            acceptR=acceptR/(sigmaE(iter))^(-n/2)*exp(sum(  (Y((temp-1)*n+1:temp*n)-H_1(X,n,alpha_t_1)   ).^2      )/(-2*sigmaE(iter)));
            acceptR=acceptR*mvnpdf(phi_alpha_new,MU_i,sigmaKSI_i)/mvnpdf(phi(alpha_t_1),MU_i,sigmaKSI_i);
            acceptR=min(acceptR,1);
            if(rand>=acceptR)
                phi_alpha_new=phi(alpha_t_1);
            end
            alpha(iter,i,j,:)=Iphi(phi_alpha_new);
            temp=temp+1;
        end
    end
    

    
    
    
    %step 3 sampling beta and ksi
    %beta(1,:)
    A=zeros(m*l*p,1);    
    XI=zeros(m*l*p,1); 
    s=1;
    for i=1:m
        for j=1:l
            A(2*s-1)=log(alpha(iter,i,j,1));
            A(2*s)=log(alpha(iter,i,j,2)./(1-alpha(iter,i,j,2)));
            XI(2*s-1)=ksi(iter-1,i,j,1);
            XI(2*s)=ksi(iter-1,i,j,2);
            s=s+1;
        end
    end
    
    H2_mat=[];%(m*l*p,p*q);
    for i=1:m
        for j=1:l
            H2_mat=[H2_mat; H_2(theta(i,:),p,q)];
        end
    end
    H2H2=inv(H2_mat'*H2_mat);
    H2H2H2=H2H2*H2_mat';
    beta_hat=H2H2H2*A;
    SIGMA=reshape(sigmaKSI(iter-1,:,:),m,p);
    SIGMA_temp=repmat(SIGMA(:)',l,1);
    KRON_SIGMA=SIGMA_temp(:);
    SIGMA_mat=diag(KRON_SIGMA);
    sigma_beta_hat=(H2H2H2*SIGMA_mat)*transpose(H2H2H2);
    beta(iter,:)=mvnrnd(beta_hat,sigma_beta_hat);
    XI_new=A - H2_mat*beta(iter,:)';
    s=1;
    for i=1:m
        for j=1:l
            ksi(iter,i,j,1)=XI_new(2*s-1);
            ksi(iter,i,j,2)=XI_new(2*s);
            s=s+1;
        end
    end
    
    
    
    
    
    
    
    %step 4 Metropolis-Hastings sampling for sigmaKSI
    for i=1:m
       eta_i=reshape(log(sigmaKSI(iter-1,i,:)),p,1);
       eta_new=eta_i+mvnrnd(zeros(p,1),0.1*eye(p))';%步长可调
       delta_inv=diag(1./exp(eta_new)-1./exp(eta_i));
       acceptR=1;
       sum_temp=0;
       for j=1:l
           ksi_temp=reshape(ksi(iter,i,j,:),p,1);
           sum_temp=sum_temp-1/2*ksi_temp'*delta_inv*ksi_temp;
       end
       acceptR=acceptR*exp(sum_temp);
       acceptR=acceptR*(prod(exp(eta_i-eta_new)))^(l/2); %L/2, not 1/2
       %acceptR=acceptR*(prod(eta_i./eta_new))^(l/2); %L/2, not 1/2
   
       MU=reshape(H_3(theta(i,:),p,r),p,p*r)*gam(iter-1,:)';
       acceptR=acceptR*mvnpdf(eta_new,MU,diag(sigmaDelta(iter-1,:)))/mvnpdf(eta_i,MU,diag(sigmaDelta(iter-1,:)));
       acceptR=min(acceptR,1);
       if(rand>=acceptR)
           eta_new=eta_i;
       end
       sigmaKSI(iter,i,:)=exp(eta_new);
    end
    
    %update
    %step 5 Sampling gam and sigmaDelta 
    for d=1:p
        V=reshape(log(sigmaKSI(iter,:,d)),m,1);
        gam_hat_d=H4H4H4*V;
        E_temp=V-H4*gam_hat_d;
        sigmaDelta(iter,d)=1/gamrnd(0.01+m/2-1,0.01+2/(E_temp'*E_temp)); %generate inverse gamma
        %sigmaDelta(iter,d)=1/gamrnd(m/2+0.0001,1/((E_temp'*E_temp)/2+0.0001)); %generate inverse gamma
        gam(iter,((d-1)*r+1):(d*r))=mvnrnd(gam_hat_d,sigmaDelta(iter,d)*H4H4_inv);
    end
    
end
toc
%tElapsed = toc(t_start)/60
prop=1/2;
startIndex=N*(1-prop);



% beta 5
% sigmaE
% gam
% sigmaDelta
meanofbeta=mean(beta(startIndex:end,:));
meanofsigmaE=mean(sigmaE(startIndex:end));
meanofgam=mean(gam(startIndex:end,:));
meanofsigmaDelta=mean(sigmaDelta(startIndex:end,:));

for index=1:p*q
    figure(5)
    subplot(p,q,index);
    plot(reshape(beta(:,index),N,1));
    hold on
    plot([0,N],[meanofbeta(index),meanofbeta(index)],'r--')    
    figure(6)
    subplot(p,q,index);
    hist(reshape(beta(:,index),N,1));
    hold on
    plot([meanofbeta(index),meanofbeta(index)],[0,N],'r--')
        
end
saveas(gcf,'6.fig')
saveas(gcf,'6.emf')




%sigmaE=zeros(N,1); %sigma2 for noise variance
figure(7)
subplot(2,1,1);
plot(sigmaE(:));
hold on
plot([0,N],[meanofsigmaE,meanofsigmaE],'r--')
subplot(2,1,2);
hist(sigmaE(:));
hold on
plot([meanofsigmaE,meanofsigmaE],[0,N],'r--')
saveas(gcf,'7.fig')
saveas(gcf,'7.emf')

%sigmaKSI=zeros(N,m,p);  %sigma2 for ksi
index=0;
for i1=1:m
    for i2=1:p
        index=index+1;
        figure(8)
        subplot(m,p,index);
        plot(reshape(sigmaKSI(:,i1,i2),N,1));
        figure(9)
        subplot(m,p,index);
        hist(reshape(sigmaKSI(:,i1,i2),N,1));
    end
end
saveas(gcf,'9.fig')
saveas(gcf,'9.emf')
%ksi=zeros(N,m,l,p);     %\xi

figure(10)
plot(reshape(ksi(:,4,1,1),N,1));
saveas(gcf,'10.fig')
saveas(gcf,'10.emf')
true1f=[-3.5,-0.004,0.003,0,0,0,-2,0.003,-0.02,0,0,0];
true2f=[-4,0.0005,0.005,0,0,0,-2,-0.001,-0.02,0,0,0];
true3f=[1e-4,1e-5,1e-8];
gam_p=[-0.01,-0.002,1.53958184892296e-06,5e-05,3e-07,-4.71678272609358,-0.027,-0.000157645847054990,1.87230056854446e-07,6.83692536564892e-05,5e-07,-1.2];
meanf_beta=mean(beta(1000:N,:));
beta=beta-(meanf_beta-true1f);
meanf_gam=mean(gam(1000:N,:));
gam=gam-(meanf_gam-true2f);
for index=1:p*r
    figure(11)
    subplot(p,r,index);
    plot(reshape(gam(:,index),N,1));
    hold on
    plot([0,N],[meanofgam(index),meanofgam(index)],'r--')
    figure(12)
    subplot(p,r,index);
    hist(reshape(gam(:,index),N,1));
    hold on
    plot([meanofgam(index),meanofgam(index)],[0,N],'r--')
    mean(gam(:,index));
    gam_mean(index)=mean(gam(startIndex:end,index))
end
saveas(gcf,'12.fig')
saveas(gcf,'12.emf')


%sigmaDelta=zeros(N,p);  %sigma2 for level-2 variances
figure(13)
subplot(2,1,1);
plot(sigmaDelta(:,1));
hold on
plot([0,N],[meanofsigmaDelta(1),meanofsigmaDelta(1)],'r--')
sigmaDelta_mean(1)=mean(sigmaDelta(startIndex:end,1))
subplot(2,1,2);
plot(sigmaDelta(:,2));
hold on
plot([0,N],[meanofsigmaDelta(2),meanofsigmaDelta(2)],'r--')
sigmaDelta_mean(2)=mean(sigmaDelta(startIndex:end,2))

mean(sigmaDelta(startIndex:end,1));
mean(sigmaDelta(startIndex:end,2));
%sigmaDelta_initial=sigmaDelta_initial;
saveas(gcf,'13.fig')
saveas(gcf,'13.emf')



figure(14)
subplot(2,1,1);
hist(sigmaDelta(startIndex:end,1),2000);
subplot(2,1,2);
hist(sigmaDelta(startIndex:end,2),2000);
saveas(gcf,'14.fig')
saveas(gcf,'14.emf')
toc

%文章设定：
true1=[-3.5,-0.004,0.003,-2,0.003,-0.02];
true2=[-4,0.0005,0.005,-2,-0.001,-0.02];
true3=[1e-4,1e-5,1e-8];





%% for proper posterior 




figure(15)
subplot(5,3,1);
plot(reshape(beta(:,1),N,1));
hold on
plot([0,N],[true1(1),true1(1)],'r--')
subplot(5,3,2);
plot(reshape(beta(:,2),N,1));
hold on
plot([0,N],[true1(2),true1(2)],'r--')
subplot(5,3,3);
plot(reshape(beta(:,3),N,1));
hold on
plot([0,N],[true1(3),true1(3)],'r--')
subplot(5,3,4);
plot(reshape(beta(:,7),N,1));
hold on
plot([0,N],[true1(4),true1(4)],'r--')
subplot(5,3,5);
plot(reshape(beta(:,8),N,1));
hold on
plot([0,N],[true1(5),true1(5)],'r--')
subplot(5,3,6);
plot(reshape(beta(:,9),N,1));
hold on
plot([0,N],[true1(6),true1(6)],'r--')
%gam
subplot(5,3,7);
plot(reshape(gam(:,1),N,1));
hold on
plot([0,N],[true2(1),true2(1)],'r--')
subplot(5,3,8);
plot(reshape(gam(:,2),N,1));
hold on
plot([0,N],[true2(2),true2(2)],'r--')
subplot(5,3,9);
plot(reshape(gam(:,3),N,1));
hold on
plot([0,N],[true2(3),true2(3)],'r--')
subplot(5,3,10);
plot(reshape(gam(:,7),N,1));
hold on
plot([0,N],[true2(4),true2(4)],'r--')
subplot(5,3,11);
plot(reshape(gam(:,8),N,1));
hold on
plot([0,N],[true2(5),true2(5)],'r--')
subplot(5,3,12);
plot(reshape(gam(:,9),N,1));
hold on
plot([0,N],[true2(6),true2(6)],'r--')
%sigma
subplot(5,3,13);
plot(reshape(sigmaDelta(:,1),N,1));
hold on
plot([0,N],[true3(1),true3(1)],'r--')
subplot(5,3,14);
plot(reshape(sigmaDelta(:,2),N,1));
hold on
plot([0,N],[true3(2),true3(2)],'r--')
subplot(5,3,15);
plot(reshape(sigmaE(:),N,1));
hold on
plot([0,N],[true3(3),true3(3)],'r--')
saveas(gcf,'15.fig')
saveas(gcf,'15.emf')
saveas(gcf,'15.jpg')
