%HLM with Variance Heterogeneity%
%------writen by Yuanyuan(2020)-----%

clc
clear
close all

%linearFitting_only_beta;

%clearvars -except b x AttenSelection beta1 beta2 gam1 gam2 sigma_error_initial ksi_initial sigmaKSI_initial sigmaDelta_initial
%close all

t_start=tic;


prop=1/4; % proportion of gibbs samplers to use


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% linearFitting2

%--dimension setting-------%
m=16;  %for input material parameters  m=16 组参数 i
l=16;  %for locations or replication  l=16 16张图 j
n=91;  %for x(frequency) n=很大 
%[m,l,n]=size(AttenSelection);  %？？？？？？？？？？？？？                         
p=2;  %for alpha_ij
q=4;  %for beta dimension
r=4;  %for gamma 
theta=zeros(m,2);    %material input
theta=[100,800;200,300;200,800;200,1300;200,2000;200,3000;300,300;300,800;300,1300;300,2000;300,3000;400,300;400,800;400,1300;400,2000;400,3000]
%assign theta

%observations

docnum=['02';'06';'07';'08';'09';'10';'11';'12';'13';'14';'15';'16';'17';'18';'19';'20'];
Y=[];
for i=1:m    
    filename=['tpcf_' num2str(docnum(i,:)) '.csv'];
    Y_temp=csvread(filename);
    Y=[Y;Y_temp(2:17,:)];
end
X=Y_temp(1,:)';  %列
Y=Y';
Y=Y(:); %列
    
%iteration number
%N=40000;
N=2;
%---parameters------------%
beta=zeros(N,p*q);   %beta=(beta_1^t,...beta_p^T)')
sigmaE=zeros(N,1); %sigma2 for noise variance
ksi=zeros(N,m,l,p);     %\xi
sigmaKSI=zeros(N,m,p);  %sigma2 for ksi
gam=zeros(N,p*r);
sigmaDelta=zeros(N,p);  %sigma2 for level-2 variances
alpha=zeros(N,m,l,p);



2.72E-02	8.20E-06	-3.61E-06	-5.05E-09
8.29E-01	6.21E-04	-1.77E-04	2.10E-07
3.38E-04	-6.65E-07	-1.74E-07	3.84E-10
6.74E-05	-1.38E-06	2.13E-06	-2.61E-09

beta1=[-0.000174646,-1.31423e-05,-1.48602e-08,3.50449e-07,3.83503e-09,0.0544238];
beta2=[0.00460975,-0.000200356,4.3674e-07,-7.70136e-06,-1.57011e-08,0.368937];
gam1=[-1.26135e-06,-4.52638e-07,-4.25883e-10,4.89297e-09,1.1734e-10,0.000599707];
gam2=[4.91669e-06,-3.2943e-06,1.75461e-09,-2.65689e-08,1.01832e-09,0.00488845];


%beta1=[-0.000174646,-1.31423e-05,-1.48602e-08,3.50449e-07,3.83503e-09,0.0544238];
%beta2=[0.00460975,-0.000200356,4.3674e-07,-7.70136e-06,-1.57011e-08,0.368937];
%gam1=[-1.26135e-06,-4.52638e-07,-4.25883e-10,4.89297e-09,1.1734e-10,0.000599707];
%gam2=[4.91669e-06,-3.2943e-06,1.75461e-09,-2.65689e-08,1.01832e-09,0.00488845];
%initial value for Gibbs

beta(1,:)=[beta1,beta2];   %beta=(beta_1^t,...beta_p^T)')
gam(1,:)=[gam1,gam2];
sigmaE(1)=0.0002149949371406595; %sigma2 for noise variance
sigmaDelta(1,:)=[9.173744750793996e-08,1.995555226461552e-06];
load('ksi_initial.mat')
load('sigmaKSIinitial.mat')
load('alpha_initial.mat')
ksi(1,:,:,:)=ksi_initial;
alpha(1,:,:,:)=alpha_initial;
sigmaKSI(1,:,:)=sigmaKSIinitial;  %sigma2 for ksi
%eta_temp=log(sigmaKSI(1,m,p));%rewrite to m*p




% h/H functions  H_1(X,n,a,b)(n,1)
% H_mat=H(X,theta,m,l,n,p,q);
% HH_inv=inv(H_mat'*H_mat);
% HHH=HH_inv*H_mat';
% H1_mat=H_1(X,n,a,b);  %??????
% H1H1=H1_mat'*H1_mat;
%H2(theta_i) p*pq



% H2_mat=zeros(m,p,p*q);
% H3_mat=zeros(m,p,p*r);
% for i=1:m
%     %H2_mat(i,:,:)=H_2(theta(i,:),p,q);
%     %H3_mat(i,:,:)=H_3(theta(i,:),p,r);
% end



H4=H_4(theta,m,r);
H4H4_inv=inv(H4'*H4);
H4H4H4=H4H4_inv*H4';






tic

for iter=2:N
    %if(mod(iter,1000)==0)
    iter
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
    sigmaE(iter)=1/gamrnd(IG_a,IG_b); %generate inverse gamma
    %beta(iter,:)=mvnrnd(beta_hat,sigmaE(iter)*HH_inv);
    


%step 2 sampling alpha Metropolis-Hastings sampling
   
    temp=1;
    for i=1:m
        MU_i=reshape(H_2(theta(i,:),p,q),p,p*q)*beta(iter-1,:)'; %先验
        
        sigmaKSI_i=diag(reshape(sigmaKSI(iter-1,i,:),2,1));
        for j=1:l
            alpha_t_1=reshape(alpha(iter-1,i,j,:),p,1);           
            alpha_new=alpha_t_1+mvnrnd(zeros(p,1),0.001*eye(p))'; %%%步长可调节
            while alpha_new(1)<0|| alpha_new(2)<0
                alpha_new=alpha_t_1+mvnrnd(zeros(p,1),0.001*eye(p))'; %%%步长可调节             
            end
            alpha_new;              
            
            
            
            
            
            acceptR=1;
            acceptR=acceptR*(sigmaE(iter))^(-n/2)*exp(sum(  (Y((temp-1)*n+1:temp*n)-H_1(X,n,alpha_new)   ).^2      )/(-2*sigmaE(iter)));
            acceptR=acceptR/(sigmaE(iter))^(-n/2)*exp(sum(  (Y((temp-1)*n+1:temp*n)-H_1(X,n,alpha_t_1)   ).^2      )/(-2*sigmaE(iter)));
            acceptR=acceptR*mvnpdf(alpha_new,MU_i,sigmaKSI_i)/mvnpdf(alpha_t_1,MU_i,sigmaKSI_i);
            acceptR=min(acceptR,1);
            if(rand>=acceptR)
                alpha_new=alpha_t_1;
            end
            alpha(iter,i,j,:)=alpha_new;
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
            A(2*s-1)=alpha(iter,i,j,1);
            A(2*s)=alpha(iter,i,j,2);
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
    
    
    %Metropolis-Hastings sampling for sigmaKSI
    
            
            
            
%             Sigma_ksi_inv=inv(diag(reshape(sigmaKSI(iter-1,i,:),p,1)));
%             ksi_temp=inv(H1H1+sigmaE(iter)*Sigma_ksi_inv);
%             y=Y((((i-1)*l+j-1)*n+1):(((i-1)*l+j)*n));
%             ksi_hat_ij= ksi_temp*(H1_mat'*(y-H1_mat*reshape(H2_mat(i,:,:),p,p*q)*beta(iter,:)'));
%             Sigma_ksi_hat_ij=ksi_temp*sigmaE(iter);
%             ksi(iter,i,j,:)=mvnrnd(ksi_hat_ij,Sigma_ksi_hat_ij);
            
%         end
%     end
    
    

    
    
    
    
    
%step 4 Metropolis-Hastings sampling for sigmaKSI
    for i=1:m
       eta_i=reshape(sigmaKSI(iter-1,i,:),p,1);
       eta_new=eta_i+mvnrnd(zeros(p,1),0.00001*eye(p))';%步长可调
       delta_inv=diag(1./eta_new-1./eta_i);
       acceptR=1;
       sum_temp=0;
       for j=1:l
           ksi_temp=reshape(ksi(iter,i,j,:),p,1);
           sum_temp=sum_temp-1/2*ksi_temp'*delta_inv*ksi_temp;
       end
       acceptR=acceptR*exp(sum_temp);
       %acceptR=acceptR*(prod(exp(eta_i-eta_new)))^(l/2); %L/2, not 1/2
       acceptR=acceptR*(prod(eta_i./eta_new))^(l/2); %L/2, not 1/2
   
       MU=reshape(H_3(theta(i,:),p,r),p,p*r)*gam(iter-1,:)';
       acceptR=acceptR*mvnpdf(eta_new,MU,diag(sigmaDelta(iter-1,:)))/mvnpdf(eta_i,MU,diag(sigmaDelta(iter-1,:)));
       acceptR=min(acceptR,1);
       if(rand>=acceptR)
           eta_new=eta_i;
       end
       sigmaKSI(iter,i,:)=eta_new;
    end
       
    %step 5 Sampling gam and sigmaDelta
    for d=1:p
        V=reshape(sigmaKSI(iter,:,d),m,1);
        gam_hat_d=H4H4H4*V;
        E_temp=V-H4*gam_hat_d;
        sigmaDelta(iter,d)=1/gamrnd(m/2-1,2/(E_temp'*E_temp)); %generate inverse gamma
        %sigmaDelta(iter,d)=1/gamrnd(m/2+0.0001,1/((E_temp'*E_temp)/2+0.0001)); %generate inverse gamma
        gam(iter,((d-1)*r+1):(d*r))=mvnrnd(gam_hat_d,sigmaDelta(iter,d)*H4H4_inv);
    end
end
toc



%%

%tElapsed = toc(t_start)/60
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
%sigma_error_initial=sigma_error_initial
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
        hist(reshape(sigmaKSI(:,i1,i2),N,1));
    end
end

%ksi=zeros(N,m,l,p);     %\xi

figure(10)
plot(reshape(ksi(:,4,1,1),N,1));

%gam=zeros(N,p*r);
for index=1:p*r
    figure(11)
    subplot(p,r,index);
    plot(reshape(gam(:,index),N,1));
    figure(12)
    subplot(p,r,index);
    hist(reshape(gam(:,index),N,1));
    mean(gam(:,index))
    gam_mean(index)=mean(gam(startIndex:end,index));
end
[gam1,gam2]


%sigmaDelta=zeros(N,p);  %sigma2 for level-2 variances
figure(13)
subplot(2,1,1);
plot(sigmaDelta(:,1));
sigmaDelta_mean(1)=mean(sigmaDelta(startIndex:end,1));
subplot(2,1,2);
plot(sigmaDelta(:,2));
sigmaDelta_mean(2)=mean(sigmaDelta(startIndex:end,2));

mean(sigmaDelta(startIndex:end,1))
mean(sigmaDelta(startIndex:end,2))
sigmaDelta_initial=sigmaDelta_initial



figure(14)
subplot(2,1,1);
hist(sigmaDelta(startIndex:end,1),2000);
subplot(2,1,2);
hist(sigmaDelta(startIndex:end,2),2000);


%%
%inference on the input parameters

% load 1-wt.mat
% x=XX1;
% x=x(1:5);
% Atten_test=YY1(setdiff(1:26,[4,9,11,16,21,25]),1:5);

load 5-wt.mat
x=XX_5_2(1:4);
Atten_test=YY_5_2([1,4,5,6,7,12,14,16,18,22,24,26,29,31,34,35,36,38,10,41],1:4)-0.1;
figure
plot(x,Atten_test)


L=20;
n=length(x);
sigma2_ksi_test=zeros(1,2);

% 

t_start2=tic;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prop=1/4; % proportion of gibbs samplers to use
startIndex=round((1-prop)*N);

theta_s=0.3:0.01:0.7;

m_s=length(theta_s);
sigmaKSI_theta=zeros(round(prop*N),m_s,p);
for i=1:round(prop*N)
    for j=1:m_s
        for k=1:p
            sigmaKSI_theta(i,j,k)=exp(reshape(gam(startIndex+i,((k-1)*r+1):(k*r)),1,r)*h3(theta_s(j),r)+normrnd(0,sqrt(sigmaDelta(startIndex+i,k))));
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
 
 tElapsed2 = toc(t_start2)/60
 
 
 figure
 theta_s=0.1:0.01:0.9;
 bar(theta_s,Likelihood/sum(Likelihood)*100)

 AA=(Likelihood/sum(Likelihood));
 BB=theta_s';
 AA*BB


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
        for k=1:p
            sigmaKSI_theta2(k)=exp(gam_mean(((k-1)*r+1):(k*r))*h3(theta_s,r)+normrnd(0,sqrt(sigmaDelta_mean(k))));
        end
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
hist(Likelihood2(1:179),50)

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



