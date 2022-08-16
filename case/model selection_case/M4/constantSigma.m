%HLM with Variance Heterogeneity%
%------writen by yuanyuan (2020)-----%


clc
clear
close all

%linearFitting_only_beta;

clearvars -except b x AttenSelection beta1 beta2 gam1 gam2 sigma_error_initial ksi_initial sigmaKSI_initial sigmaDelta_initial
close all

t_start=tic;

gam1=gam1
gam2=gam2
beta1=beta1
beta2=beta2
sigma_error_initial=sigma_error_initial
sigmaDelta_initial=sigmaDelta_initial
prop=1/4; % proportion of gibbs samplers to use
%%%%%%%simulate data

m=length(b);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% linearFitting2

%--dimension setting-------%
m=16;  %for input material parameters  m=16 组参数 i
l=16;  %for locations or replication  l=16 16张图 j
% n=9;  %for x(frequency) n=很大 
[m,l,n]=size(AttenSelection);  %？？？？？？？？？？？？？                         
p=2;  %for alpha_ij
q=6;  %for beta dimension
r=6;  %for gamma 
theta=zeros(m,2);    %material input
theta=[100,800;200,300;200,800;200,1300;200,2000;200,3000;300,300;300,800;300,1300;300,2000;300,3000;400,300;400,800;400,1300;400,2000;400,3000]
%assign theta

%observations

docnum=[2,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20];
Y=[];
for i=1:m    
    filename=['tpcf_' num2str(docnum(i)) '.csv'];
    Y_temp=csvread(filename)
    Y=[Y;Y_temp(2:17,:)];
end
X=Y_temp(1,:)'  %列
Y=Y(:)' %列
    
%iteration number
N=40000;
%---parameters------------%

beta=zeros(N,p*q);   %beta=(beta_1^t,...beta_p^T)')
sigmaE=zeros(N,1); %sigma2 for noise variance
ksi=zeros(N,m,l,p);     %\xi
sigmaKSI=zeros(N,m,p);  %sigma2 for ksi

%initial value for Gibbs
beta1=[-0.000174646,-1.31423e-05,-1.48602e-08,3.50449e-07,3.83503e-09,0.0544238];
beta2=[0.00460975,-0.000200356,4.3674e-07,-7.70136e-06,-1.57011e-08,0.368937];
beta(1,:)=[beta1,beta2];   %beta=(beta_1^t,...beta_p^T)')
sigmaE(1)=0.0002149949371406595; %sigma2 for noise variance
load('ksi_initial.mat')
load('sigmaKSI_initial.mat')
load('alpha_initial.mat')
ksi(1,:,:,:)=ksi_initial;
sigmaKSI(1,:,:)=sigmaKSIinitial;  %sigma2 for ksi
%eta_temp=log(sigmaKSI(1,m,p));%rewrite to m*p



% h/H functions
H_mat=H(X,theta,m,l,n,p,q);
HH_inv=inv(H_mat'*H_mat);
HHH=HH_inv*H_mat';
H1_mat=H_1(X,p,n);
H1H1=H1_mat'*H1_mat;
%H2(theta_i) p*pq
H2_mat=zeros(m,p,p*q);
H3_mat=zeros(m,p,p*r);
for i=1:m
    H2_mat(i,:,:)=H_2(theta(i,:),p,q);
    H3_mat(i,:,:)=H_3(theta(i,:),p,r);
end
H4=H_4(theta,m,r);
H4H4_inv=inv(H4'*H4);
H4H4H4=H4H4_inv*H4';

for iter=2:N
    if(mod(iter,1000)==0)
    iter
    end
%step 1 sampling sigmaE and beta
    KSI=XI(p,n,reshape(ksi(iter-1,:,:,:),m,l,p),m,l,H1_mat);
    Y_KSI=Y-KSI;
    beta_hat=HHH*Y_KSI;
    IG_a=m*n*l/2-1;
    IG_b=sum((Y_KSI-H_mat*beta_hat).^2)/2;
    sigmaE(iter)=1/gamrnd(IG_a,1/IG_b); %generate inverse gamma
    beta(iter,:)=mvnrnd(beta_hat,sigmaE(iter)*HH_inv);
    
%step 2 sampling ksi
   
    for i=1:m
        for j=1:l
            Sigma_ksi_inv=inv(diag(reshape(sigmaKSI(iter-1,i,:),p,1)));
            ksi_temp=inv(H1H1+sigmaE(iter)*Sigma_ksi_inv);
            y=Y((((i-1)*l+j-1)*n+1):(((i-1)*l+j)*n));
            ksi_hat_ij= ksi_temp*(H1_mat'*(y-H1_mat*reshape(H2_mat(i,:,:),p,p*q)*beta(iter,:)'));
            Sigma_ksi_hat_ij=ksi_temp*sigmaE(iter);
            ksi(iter,i,j,:)=mvnrnd(ksi_hat_ij,Sigma_ksi_hat_ij);
        end
    end
    
%step 3  sigmaKSI

for k=1:p
    ksi_temp=reshape(ksi(iter,:,:,k),l*m,1);
    ksi_2=sum((ksi_temp-mean(ksi_temp)).^2)/2;
    ksi_1=(m*l-2)/2;
    temp=1/gamrnd(ksi_1,1/ksi_2);
    for i=1:m
     sigmaKSI(iter,i,k)=temp;
    end
end
       
end

%%
prop=1/4; % proportion of gibbs samplers to use
startIndex=round((1-prop)*N);


%beta=zeros(N,p*q);   %beta=(beta_1',...beta_p')'
for index=1:p*q
    figure(5)
    subplot(p,q,index);
    plot(reshape(beta(:,index),N,1));
    figure(6)
    subplot(p,q,index);
    hist(reshape(beta(:,index),N,1));
    mean(beta(N/2:end,index))
    beta_mean(index)=mean(beta(startIndex:end,index));
end
[beta1,beta2]

%sigmaE=zeros(N,1); %sigma2 for noise variance
figure(7)
subplot(2,1,1);
plot(sigmaE(:));
subplot(2,1,2);
hist(sigmaE(:));
mean_sigmaE=mean(sigmaE(N/2:end))
sigma_error_initial=sigma_error_initial
sigmaE_mean=mean(sigmaE(startIndex:end));

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
        hist(reshape(sigmaKSI(:,i1,i2),N,1),100);
        mean(sigmaKSI(startIndex:end,i1,i2))
    end
end

%ksi=zeros(N,m,l,p);     %\xi

figure(10)
plot(reshape(ksi(:,4,1,1),N,1));



%%
%inference on the input parameters
%%
%inference on the input parameters

x=x;
Atten_test=reshape(AttenSelection(8,:,:),l,n);
L=20;
n=length(x);
sigma2_ksi_test=zeros(1,2);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prop=1/4; % proportion of gibbs samplers to use
startIndex=round((1-prop)*N);

theta_s=0.6:0.01:0.9;

m_s=length(theta_s);
sigmaKSI_theta=zeros(round(prop*N),m_s,p);
for i=1:round(prop*N)
    for j=1:m_s
        for k=1:p
            sigmaKSI_theta(i,j,k)=sigmaKSI(startIndex+i,1,k);
        end
    end
end

%calculate the Likelihood
Likehood=zeros(1,m_s);
for i=1:m_s
    theta_s(i)
    for g=1:round(prop*N)
        f_g=0;
        for j=1:L
            mu_temp=H_1(x,p,n)*H_2(theta_s(i),p,q)*beta(startIndex+g,:)';
            sigma_temp=H_1(x,p,n)*diag(reshape(sigmaKSI_theta(g,i,:),p,1))*H_1(x,p,n)'+sigmaE(g+startIndex)*eye(n);
            f_g=f_g+log(mvnpdf(reshape(Atten_test(j,:),n,1),mu_temp,sigma_temp));
        end
        Likehood(i)=Likehood(i)+exp(f_g);
    end
    
end
 Likehood= Likehood/round(prop*N);
 Likelihood=Likehood;
 
 figure
 bar(theta_s,Likelihood/sum(Likelihood)*100)



%%
%get distribution for log-likelihood 
%%%%%%%simulate data

prop=1/4; % proportion of gibbs samplers to use
startIndex=round((1-prop)*N);
b_test=0.68;
theta_s=b_test;


sigmaKSI_theta2=zeros(p,1);

for k=1:p
    sigmaKSI_theta2(k)=exp(gam_mean(((k-1)*r+1):(k*r))*h3(theta_s,r)+normrnd(0,sqrt(sigmaDelta_mean(k))));
end



for NN=1:1000
    NN

    x=x;
    L=20;
    n=length(x);
    sigma2_ksi_test=zeros(1,2);
    A=zeros(1,L);
    B=zeros(1,L);
    
    Atten_test=zeros(L,n);
    for j=1:L
        A(j)=beta_mean(1)*b_test+beta_mean(2)+normrnd(0,sqrt(sigmaKSI_theta2(1)));
        B(j)=beta_mean(3)*b_test+beta_mean(4)+normrnd(0,sqrt(sigmaKSI_theta2(2)));
        for k=1:n
            Atten_test2(j,k)=x(k)*A(j)+B(j)+normrnd(0,sqrt(sigmaE_mean));
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %calculate the Likelihood
    
    f_g=0;
    for j=1:L
        temp=0;
        for NNN=1:1000
            for k=1:p
                sigmaKSI_theta2(k)=exp(gam_mean(((k-1)*r+1):(k*r))*h3(theta_s,r)+normrnd(0,sqrt(sigmaDelta_mean(k))));
            end
            mu_temp=H_1(x,p,n)*H_2(theta_s,p,q)*beta_mean';
            sigma_temp=H_1(x,p,n)*diag(sigmaKSI_theta2)*H_1(x,p,n)'+sigmaE_mean*eye(n);
            temp=temp+(mvnpdf(reshape(Atten_test2(j,:),n,1),mu_temp,sigma_temp));
        end
        f_g=f_g+log(temp/1000);
    end
    Likelihood2(NN)=f_g;
end
figure
hist(Likelihood2(1:1000),50)

%calculate the experimental data likelihood
f_g=0;
for j=1:L
    temp=0;
    for NNN=1:1000
        for k=1:p
            sigmaKSI_theta2(k)=exp(gam_mean(((k-1)*r+1):(k*r))*h3(theta_s,r)+normrnd(0,sqrt(sigmaDelta_mean(k))));
        end
        mu_temp=H_1(x,p,n)*H_2(theta_s,p,q)*beta_mean';
        sigma_temp=H_1(x,p,n)*diag(sigmaKSI_theta2)*H_1(x,p,n)'+sigmaE_mean*eye(n);
        temp=temp+(mvnpdf(reshape(Atten_test(j,:),n,1),mu_temp,sigma_temp));
    end
    f_g=f_g+log(temp/1000);
end



for i=1:1000
    ss(i)=normrnd(1,1);
    ll(i)=normpdf(ss(i),1,1);
end
bar(ll)

