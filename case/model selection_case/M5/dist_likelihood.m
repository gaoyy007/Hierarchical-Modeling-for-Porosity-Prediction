%HLM with Variance Heterogeneity%
%------writen by Jianguo Wu (2015)-----%

clc
clear
close all

linearFitting_only_beta;

clearvars -except x AttenSelection beta1 beta2 gam1 gam2 sigma_error_initial ksi_initial sigmaKSI_initial sigmaDelta_initial
close all
gam1=gam1
gam2=gam2
beta1=beta1
beta2=beta2
sigma_error_initial=sigma_error_initial
sigmaDelta_initial=sigmaDelta_initial
prop=1/4; % proportion of gibbs samplers to use
%%%%%%%simulate data
b=0.1:0.025:0.975;
x=2:0.2:3;
m=length(b);
l=20;
n=length(x);
sigma2_ksi=zeros(m,2);
A=zeros(m,l);
B=zeros(m,l);

AttenSelection=zeros(m,l,n);
for i=1:m
    sigma2_ksi(i,1)=exp(gam1(1)*b(i)+gam1(2)+normrnd(0,sqrt(sigmaDelta_initial(1))));
    sigma2_ksi(i,2)=exp(gam2(1)*b(i)+gam2(2)+normrnd(0,sqrt(sigmaDelta_initial(2))));
    for j=1:l
        A(i,j)=beta1(1)*b(i)+beta1(2)+normrnd(0,sqrt(sigma2_ksi(i,1)));
        B(i,j)=beta2(1)*b(i)+beta2(2)+normrnd(0,sqrt(sigma2_ksi(i,2)));
        for k=1:n
            AttenSelection(i,j,k)=x(k)*A(i,j)+B(i,j)+normrnd(0,sqrt(sigma_error_initial));
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
linearFitting2

%--dimension setting-------%
% m=9;  %for input material parameters 
% l=20;  %for locations or replication 
% n=9;  %for x(frequency)
[m,l,n]=size(AttenSelection);
p=2;  %for alpha_ij
q=2;  %for beta_p
r=2;  %for gamma
theta=zeros(m,1);    %material input

%assign theta
for i=1:m
    %theta(i,:)=0.1*i;
    theta(i,:)=b(i);
end
%observations

X=x;
Y=zeros(m*l*n,1);
index=0;
for i=1:m
    for j=1:l
        for k=1:n
            index=index+1;
            Y(index)=AttenSelection(i,j,k);
        end
    end
end
    
%iteration number
N=40000;
%---parameters------------%

beta=zeros(N,p*q);   %beta=(beta_1^t,...beta_p^T)')
sigmaE=zeros(N,1); %sigma2 for noise variance
ksi=zeros(N,m,l,p);     %\xi
sigmaKSI=zeros(N,m,p);  %sigma2 for ksi
gam=zeros(N,p*r);
sigmaDelta=zeros(N,p);  %sigma2 for level-2 variances

%initial value for Gibbs
%beta(1,:)=[0.2169,0.1322,-0.4046,0.5390];   %beta=(beta_1^t,...beta_p^T)')
beta(1,:)=[beta1,beta2];   %beta=(beta_1^t,...beta_p^T)')
ksi(1,:,:,:)=ksi_initial;
gam(1,:)=[gam1,gam2];
sigmaKSI(1,:,:)=sigmaKSI_initial;  %sigma2 for ksi
%eta_temp=log(sigmaKSI(1,m,p));%rewrite to m*p
sigmaE(1)=sigma_error_initial; %sigma2 for noise variance
sigmaDelta(1,:)=sigmaDelta_initial;


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
    
%step 3 Metropolis-Hastings sampling for sigmaKSI
    for i=1:m
       eta_i=reshape(log(sigmaKSI(iter-1,i,:)),p,1);
       eta_new=eta_i+mvnrnd(zeros(p,1),2*eye(p))';
       delta_inv=diag(1./exp(eta_new)-1./exp(eta_i));
       acceptR=1;
       sum_temp=0;
       for j=1:l
           ksi_temp=reshape(ksi(iter,i,j,:),p,1);
           sum_temp=sum_temp-1/2*ksi_temp'*delta_inv*ksi_temp;
       end
       acceptR=acceptR*exp(sum_temp);
       acceptR=acceptR*(prod(exp(eta_i-eta_new)))^(l/2); %L/2, not 1/2
       MU=reshape(H3_mat(i,:,:),p,p*r)*gam(iter-1,:)';
       acceptR=acceptR*mvnpdf(eta_new,MU,diag(sigmaDelta(iter-1,:)))/mvnpdf(eta_i,MU,diag(sigmaDelta(iter-1,:)));
       acceptR=min(acceptR,1);
       if(rand>=acceptR)
           eta_new=eta_i;
       end
       sigmaKSI(iter,i,:)=exp(eta_new);
    end
       
    %step 4 Sampling gam and sigmaDelta
    for d=1:p
        V=reshape(log(sigmaKSI(iter,:,d)),m,1);
        gam_hat_d=H4H4H4*V;
        E_temp=V-H4*gam_hat_d;
        sigmaDelta(iter,d)=1/gamrnd(m/2-1,2/(E_temp'*E_temp)); %generate inverse gamma
        %sigmaDelta(iter,d)=1/gamrnd(m/2+0.0001,1/((E_temp'*E_temp)/2+0.0001)); %generate inverse gamma
        gam(iter,((d-1)*r+1):(d*r))=mvnrnd(gam_hat_d,sigmaDelta(iter,d)*H4H4_inv);
    end
end

%%
%beta=zeros(N,p*q);   %beta=(beta_1',...beta_p')'
for index=1:p*q
    figure(5)
    subplot(p,q,index);
    plot(reshape(beta(:,index),N,1));
    figure(6)
    subplot(p,q,index);
    hist(reshape(beta(:,index),N,1));
    mean(beta(N/2:end,index))
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
end
[gam1,gam2]


%sigmaDelta=zeros(N,p);  %sigma2 for level-2 variances
figure(13)
subplot(2,1,1);
plot(sigmaDelta(:,1));
subplot(2,1,2);
plot(sigmaDelta(:,2));

mean(sigmaDelta(N/100:end,1))
mean(sigmaDelta(N/100:end,2))
sigmaDelta_initial=sigmaDelta_initial

figure(14)
subplot(2,1,1);
hist(sigmaDelta(:,1));
subplot(2,1,2);
hist(sigmaDelta(:,2));


%%
%get distribution for log-likelihood 
%%%%%%%simulate data
for NN=1:1000
    NN
    b_test=0.5;
    x=2:0.2:3;
    L=20;
    n=length(x);
    sigma2_ksi_test=zeros(1,2);
    A=zeros(1,L);
    B=zeros(1,L);
    
    Atten_test=zeros(L,n);
    sigma2_ksi(1)=exp(gam1(1)*b_test+gam1(2)+normrnd(0,sqrt(sigmaDelta_initial(1))));
    sigma2_ksi(2)=exp(gam2(1)*b_test+gam2(2)+normrnd(0,sqrt(sigmaDelta_initial(2))));
    for j=1:L
       % A(j)=beta1(1)*b_test+beta1(2)+normrnd(0,sqrt(sigma2_ksi_test(1)));
        A(j)=6*beta1(1)*b_test+beta1(2)+normrnd(0,sqrt(sigma2_ksi_test(1)));
        B(j)=beta2(1)*b_test+beta2(2)+normrnd(0,sqrt(sigma2_ksi_test(2)));
        for k=1:n
            Atten_test(j,k)=x(k)*A(j)+B(j)+normrnd(0,sqrt(sigma_error_initial));
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    prop=1/4; % proportion of gibbs samplers to use
    startIndex=round((1-prop)*N);
    
    theta_s=0.5;
    
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
        for g=1:round(prop*N)
            f_g=0;
            for j=1:L
                mu_temp=H_1(X,p,n)*H_2(theta_s(i),p,q)*beta(startIndex+g,:)';
                sigma_temp=H_1(X,p,n)*diag(reshape(sigmaKSI_theta(g,i,:),p,1))*H_1(X,p,n)'+sigmaE(g+startIndex)*eye(n);
                f_g=f_g+log(mvnpdf(reshape(Atten_test(j,:),n,1),mu_temp,sigma_temp));
            end
            Likehood(i)=Likehood(i)+exp(f_g);
        end
        
    end
    Likehood= Likehood/round(prop*N);
    Likelihood(NN)=log(Likehood);
end
hist(Likelihood)














