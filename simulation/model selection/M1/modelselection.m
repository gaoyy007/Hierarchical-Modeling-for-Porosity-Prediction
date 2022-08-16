%%%%%%%%%%%%%%%%%%%%%%%%model selection 
%M1 h_2 (θ)β=β_0+β_1 P+β_2 V+β_3 P^2+β_4 P*V

%%%%%%%%%%%%%%%%%%%%%%%


clc
clear
close all


%%
%simulation
tic;
prop=1/4; % proportion of gibbs samplers to use

filename=['tpcf_02.csv'];
Y_temp=csvread(filename);
X=Y_temp(1,1:30)';  %列
m=10; %10个样品
l=5; %5张图
p=2; %a b 两个参数
n=30;
prop=1/4; % proportion of gibbs samplers to use
theta=zeros(m,2); 
theta=csvread('theta.csv');
q=6;%q_s=3;
r=6;%r_s=3;

%model init
beta_s=[0.0526,-8.1706e-05,-2.4123e-05,3.5048e-08,4.6631e-09,1.9605e-08,0.3689,0.0046,-2.0036e-04,-7.7014e-06,-1.5701e-08,4.3674e-07];    
sigmaDelta_s=[0.1,0.1]; %!!
gam_s=[-13.4957,0.0338,-0.0031,-5.3064e-05,7.5698e-07,-1.7457e-06,-6.2032,-0.0289,0.0032,6.6612e-05,-2.2633e-07,-5.2715e-06];
sigmaE=0.0002149949371406595;


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

    alpha_s(i,:,:)=MU_i'+reshape(ksi_s(i,:,:),l,p);  

end

Y=[];

for i=1:m
    for j=1:l

        Y_temp=H_1(X,n,alpha_s(i,j,:));
        noise=normrnd(0,sigmaE,n,1);
        Y_temp=Y_temp+noise;
        Y=[Y,Y_temp];


    end
    subplot(2,5,m)
    plot(Y_temp)
end
Y=Y(:);
xlswrite('Y.xls',Y)
xlswrite('X.xls',X)

N=1000;
%---parameters------------%
beta=zeros(N,p*q);   %beta=(beta_1^t,...beta_p^T)')
sigmaE=zeros(N,1); %sigma2 for noise variance
ksi=zeros(N,m,l,p);     %\xi
sigmaKSI=zeros(N,m,p);  %sigma2 for ksi
gam=zeros(N,p*r);
sigmaDelta=zeros(N,p);  %sigma2 for level-2 variances
alpha=zeros(N,m,l,p);


% %initial value for Gibbs
% %%%%%%%%%%%%%%%%%1
beta(1,:)=beta_s;   %beta=(beta_1^t,...beta_p^T)')
gam(1,:)=gam_s;
sigmaE(1)=0.000001; %sigma2 for noise variance
%sigmaDelta(1,:)=[9.173744750793996e-08,1.995555226461552e-06];
sigmaDelta(1,:)=sigmaDelta_s;
load('ksi_initial.mat')
load('sigmaKSIinitial.mat')
load('alpha_initial.mat')
ksi(1,:,:,:)=ksi_s;
alpha(1,:,:,:)=alpha_s;
sigmaKSI(1,:,:)=sigmaKSI_s;  %sigma2 for ksi
%eta_temp=log(sigmaKSI(1,m,p));%rewrite to m*p



%%%%%%%%%%%%%%%%%%%%%%%%2

beta(1,:)=[0.01,-1e-05,-1e-05,1e-08,1e-09,1e-08,1,0.001,-1e-04,-1e-06,-1e-08,1e-07];    

gam(1,:)=[-10,0.01,-0.001,-1e-05,1e-07,-1e-06,-1,-0.01,0.001,1e-05,-1e-07,-1e-06];
sigmaE(1)=0.0002149949371406595;


%beta(1,:)=[0.001, -0.001,-1e-8, -1e-7, 1e-7,-8,-1e-2,1E-2,1e-7,-1e-6,2e-7,-5];   %beta=(beta_1^t,...beta_p^T)')
%beta(1,:)=zeros(1,12)+0.001;
%gam(1,:)=zeros(1,12)+0.001;
%sigmaE(1)=0.00001; %sigma2 for noise variance
%sigmaDelta(1,:)=[9.173744750793996e-08,1.995555226461552e-06];
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
            alpha_new=alpha_t_1+mvnrnd(zeros(p,1),[0.1,0;0,0.1])'; %%%步长可调节
            while alpha_new(1)<=0|| alpha_new(2)<=0||alpha_new(1)>=0.3||alpha_new(2)>=1
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
    
    
    %step 5 Sampling gam and sigmaDelta
    for d=1:p
        V=reshape(log(sigmaKSI(iter,:,d)),m,1);
        gam_hat_d=H4H4H4*V;
        E_temp=V-H4*gam_hat_d;
        sigmaDelta(iter,d)=1/gamrnd(m/2-1,2/(E_temp'*E_temp)); %generate inverse gamma
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

toc





%%%%%%%%%%%%%%%%%%%%%%%%simulation 
%%%%%%%%%%%%%%%%%%%%%%%para initial  sigma tune



close all




filename=['tpcf_02.csv'];
Y_temp=csvread(filename);
X=Y_temp(1,1:30)';  %列
m=10; %10个样品    %数据量过大导致inf 因此减少数据量 m=10 l=2 p=2 n=10 总共数据量 200
l=2; %5张图
p=2; %a b 两个参数
n=10;
bf=75.3003;
prop=1/4; % proportion of gibbs samplers to use
theta=zeros(m,2); 
theta=csvread('theta.csv');
q=6;%q_s=3;
r=6;%r_s=3;


beta_s=[0.0526,-8.1706e-05,-2.4123e-05,3.5048e-08,4.6631e-09,1.9605e-08,0.3689,0.0046,-2.0036e-04,-7.7014e-06,-1.5701e-08,4.3674e-07];    
sigmaDelta_s=[0.1,0.1]; %!!
gam_s=[-13.4957,0.0338,-0.0031,-5.3064e-05,7.5698e-07,-1.7457e-06,-6.2032,-0.0289,0.0032,6.6612e-05,-2.2633e-07,-5.2715e-06];
sigmaE=0.0002149949371406595;


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
    
    alpha_s(i,:,:)=MU_i'+reshape(ksi_s(i,:,:),l,p);  
    
end

Y=[];

for i=1:m
    for j=1:l
        
        Y_temp=H_1(X,n,alpha_s(i,j,:));
        noise=normrnd(0,sigmaE,n,1);
        Y_temp=Y_temp+noise;
        Y=[Y,Y_temp];
        
        
    end
    %subplot(2,5,m)
    %plot(Y_temp)
end
Y=Y(:);
%xlswrite('Y.xls',Y)
xlswrite('X.xls',X)







%%
%simulation
tic;
prop=1/4; % proportion of gibbs samplers to use


N=1000;
%---parameters------------%
beta=zeros(N,p*q);   %beta=(beta_1^t,...beta_p^T)')
sigmaE=zeros(N,1); %sigma2 for noise variance
ksi=zeros(N,m,l,p);     %\xi
sigmaKSI=zeros(N,m,p);  %sigma2 for ksi
gam=zeros(N,p*r);
sigmaDelta=zeros(N,p);  %sigma2 for level-2 variances
alpha=zeros(N,m,l,p);

%
% % %initial value for Gibbs
% % %%%%%%%%%%%%%%%%%1
% beta(1,:)=beta_s;   %beta=(beta_1^t,...beta_p^T)')
% gam(1,:)=gam_s;
% sigmaE(1)=0.000001; %sigma2 for noise variance
% %sigmaDelta(1,:)=[9.173744750793996e-08,1.995555226461552e-06];
% sigmaDelta(1,:)=sigmaDelta_s;
% load('ksi_initial.mat')
% load('sigmaKSIinitial.mat')
% load('alpha_initial.mat')
% ksi(1,:,:,:)=ksi_s;
% alpha(1,:,:,:)=alpha_s;
% sigmaKSI(1,:,:)=sigmaKSI_s;  %sigma2 for ksi
% %eta_temp=log(sigmaKSI(1,m,p));%rewrite to m*p



%%%%%%%%%%%%%%%%%%%%%%%%2

beta(1,:)=[0.01,-1e-05,-1e-05,1e-08,1e-09,1e-08,1,0.001,-1e-04,-1e-06,-1e-08,1e-07];    

gam(1,:)=[-10,0.01,-0.001,-1e-05,1e-07,-1e-06,-1,-0.01,0.001,1e-05,-1e-07,-1e-06];
sigmaE(1)=0.0002149949371406595;


%beta(1,:)=[0.001, -0.001,-1e-8, -1e-7, 1e-7,-8,-1e-2,1E-2,1e-7,-1e-6,2e-7,-5];   %beta=(beta_1^t,...beta_p^T)')
%beta(1,:)=zeros(1,12)+0.001;
%gam(1,:)=zeros(1,12)+0.001;
%sigmaE(1)=0.00001; %sigma2 for noise variance
%sigmaDelta(1,:)=[9.173744750793996e-08,1.995555226461552e-06];
sigmaDelta(1,:)=[0.001,0.001];
ksi(1,:,:,:)=ksi_s;
alpha(1,:,:,:)=alpha_s;
sigmaKSI(1,:,:)=sigmaKSI_s;  %sigma2 for ksi

H4=H_4(theta,m,r);
H4H4_inv=inv(H4'*H4);
H4H4H4=H4H4_inv*H4';



% define ML
tic
prod_pdf1=0;
prod_pdf2=0;
prod_pdf3=0;
pospdf1=0;
pospdf2=0;
pospdf3=0;
pospdf4=0;

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
    
    if iter>750
        pdf1=mvnpdf(Y,H1_mat,0.00021499);   %H1_mat 1500*1 
        prod_pdf1=prod_pdf1+log(prod(pdf1));
    else  
        prod_pdf1=prod_pdf1;
    end
    
  
    Y_H1_mat=Y-H1_mat;
    %beta_hat=HHH*Y_KSI;
    IG_a=m*n*l/2-1;
    IG_b=sum((Y_H1_mat).^2)/2;
    sigmaE(iter)=1/gamrnd(IG_a,1/IG_b); %generate inverse gamma
    %beta(iter,:)=mvnrnd(beta_hat,sigmaE(iter)*HH_inv);
    
    if iter>750
        %sigmaE(iter)=1/gamrnd(IG_a,1/IG_b);
        pospdf1_temp=gampdf(1/sigmaE(iter),IG_a,1/IG_b);
        %pospdf1_temp=gampdf(1/0.00021499,IG_a,1/IG_b)  %why 0？   pospdf1
        pospdf1=log(pospdf1_temp)+pospdf1;
    else
        pospdf1=pospdf1;
    end
    
    
    %step 2 sampling alpha Metropolis-Hastings sampling
    
    temp=1;
    for i=1:m
        
        MU_i=reshape(H_2(theta(i,:),p,q),p,p*q)*beta(iter-1,:)'; %先验
        
        sigmaKSI_i=diag(reshape(sigmaKSI(iter-1,i,:),2,1));
        for j=1:l
            
            alpha_t_1=reshape(alpha(iter-1,i,j,:),p,1);
            alpha_new=alpha_t_1+mvnrnd(zeros(p,1),[0.1,0;0,0.1])'; %%%步长可调节
            while alpha_new(1)<=0|| alpha_new(2)<=0||alpha_new(1)>=0.3||alpha_new(2)>=1
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
    
    if iter>750  % pospdf2
        pospdf2_temp=mvnpdf(meanofbeta',beta_hat,sigma_beta_hat);
        pospdf2=pospdf2+log(prod(pospdf2_temp));
    else
        pospdf2=pospdf2;
    end
    
    
    XI_new=A - H2_mat*beta(iter,:)';
    s=1;
    if iter>750
        
        
        %pdf2=mvnpdf(A,H2_mat*meanofbeta',SIGMA_mat);
        C=H2_mat*beta(iter,:)';
        pdf2=mvnpdf(A(1:10),C(1:10),SIGMA_mat(1:10,1:10));
        
        %0.422560777289687,0.421137203101778,5.851368886544799e-06
        %single 1.386993857266435e+02
        
        prod_pdf2=prod_pdf2+log(prod(pdf2));
    else
        prod_pdf2=prod_pdf2;
    end
    
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
    
    pdf3_temp=zeros(1,2);
    pospdf3_temp=zeros(1,2);
    pospdf4_temp=zeros(6,2);%这个6是有的变化的；
    
    %step 5 Sampling gam and sigmaDelta
    for d=1:p
        V=reshape(log(sigmaKSI(iter,:,d)),m,1);
        gam_hat_d=H4H4H4*V;
        E_temp=V-H4*gam_hat_d;
        sigmaDelta(iter,d)=1/gamrnd(m/2-1,2/(E_temp'*E_temp)); %generate inverse gamma
        %sigmaDelta(iter,d)=1/gamrnd(m/2+0.0001,1/((E_temp'*E_temp)/2+0.0001)); %generate inverse gamma
        gam(iter,((d-1)*r+1):(d*r))=mvnrnd(gam_hat_d,sigmaDelta(iter,d)*H4H4_inv);
        pdf3_temp(1,d)=prod(mvnpdf(V,H4*gam_hat_d,sigmaDelta(iter,d))); %log-likelihood
        pospdf3_temp(1,d)=gampdf(1/sigmaDelta(iter,d),m/2-1,2/(E_temp'*E_temp));%pos
        pospdf4_temp(1,d)=prod(mvnpdf(gam(iter,((d-1)*r+1):(d*r))',gam_hat_d,sigmaDelta(iter,d)*H4H4_inv) );
       

    
    end
    if iter>750
        pdf3=pdf3_temp(1,1)*pdf3_temp(1,2); 
        prod_pdf3=prod_pdf3+log(pdf3);
        
        pospdf3_ptemp=pospdf3_temp(1,1)*pospdf3_temp(1,2);  % temp 是两个相乘，ptemp则是前n个观测值的连乘
        pospdf3=pospdf3+log(pospdf3_ptemp);
        
        pospdf4_ptemp=pospdf4_temp(1,1)*pospdf4_temp(1,2);
        pospdf4=pospdf4+log(pospdf4_ptemp);
        
    else
        prod_pdf3=prod_pdf3;   
        pospdf3=pospdf3;   
        
         
        pospdf4=pospdf4; 
        
    end
    
 
    
    
end 


loglikelihood=prod_pdf1+prod_pdf2+prod_pdf3;
posterior=pospdf1+pospdf2+pospdf3+pospdf4;
%IBF=loglikelihood/posterior;
toc

