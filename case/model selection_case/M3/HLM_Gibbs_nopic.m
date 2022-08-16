
clc
clear
close all

%linearFitting_only_beta;

%clearvars -except b x AttenSelection beta1 beta2 gam1 gam2 sigma_error_initial ksi_initial sigmaKSI_initial sigmaDelta_initial
%close all

tic


prop=1/4; % proportion of gibbs samplers to use


%--dimension setting-------%
m=16;  %for input material parameters  m=16 组参数 i
l=16;  %for locations or replication  l=16 16张图 j
n=91;  %for x(frequency) n=很大 
%[m,l,n]=size(AttenSelection);  %？？？？？？？？？？？？？                         
p=2;  %for alpha_ij
q=6;  %for beta dimension
r=6;  %for gamma 
theta=zeros(m,2);    %material input
theta=[100,800;200,300;200,800;200,1300;200,2000;200,3000;300,300;300,800;300,1300;300,2000;300,3000;400,300;400,800;400,1300;400,2000;400,3000];

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
N=200;
%---parameters------------%
beta=zeros(N,p*q);   %beta=(beta_1^t,...beta_p^T)')
sigmaE=zeros(N,1); %sigma2 for noise variance
ksi=zeros(N,m,l,p);     %\xi
sigmaKSI=zeros(N,m,p);  %sigma2 for ksi
gam=zeros(N,p*r);
sigmaDelta=zeros(N,p);  %sigma2 for level-2 variances
alpha=zeros(N,m,l,p);

%% pre
% beta1=[-0.000174646,-1.31423e-05,-1.48602e-08,3.50449e-07,3.83503e-09,0.0544238];
% beta2=[0.00460975,-0.000200356,4.3674e-07,-7.70136e-06,-1.57011e-08,0.368937];
% %gam1=[-1.26135e-06,-4.52638e-07,-4.25883e-10,4.89297e-09,1.1734e-10,0.000599707];%nolog
% gam1=[0.00114621,-0.00206635,-3.28275e-06,1.36337e-05,5.19115e-07,-8.62905];
% %gam2=[4.91669e-06,-3.2943e-06,1.75461e-09,-2.65689e-08,1.01832e-09,0.00488845];%noloh
% gam2=[-0.00224754,-0.00102971,1.47094e-06,-4.89176e-06,2.60478e-07,-4.87428];

beta1=[5.26e-02,-8.17e-05,-2.41e-05,3.50e-08,4.66e-09,1.96e-08];
beta2=[3.69E-01,4.61E-03,-2.00E-04,-7.70E-06,-1.57E-08,4.37E-07];
%gam1=[-1.26135e-06,-4.52638e-07,-4.25883e-10,4.89297e-09,1.1734e-10,0.000599707];%nolog
gam1=[-1.35E+01,3.38E-02,-3.10E-03,-5.31E-05,7.57E-07,-1.75E-06];
%gam2=[4.91669e-06,-3.2943e-06,1.75461e-09,-2.65689e-08,1.01832e-09,0.00488845];%noloh
gam2=[-6.20E+00,-2.89E-02,3.19E-03,6.66E-05,-2.26E-07,-5.27E-06];




%initial value for Gibbs

beta(1,:)=[beta1,beta2];   %beta=(beta_1^t,...beta_p^T)')
gam(1,:)=[gam1,gam2];
sigmaE(1)=0.0002149949371406595; %sigma2 for noise variance
%sigmaDelta(1,:)=[9.173744750793996e-08,1.995555226461552e-06];
sigmaDelta(1,:)=[1.172017142009797,0.37757132424229733];
load('ksi_initial.mat')
load('sigmaKSIinitial.mat')
load('alpha_initial.mat')
ksi(1,:,:,:)=ksi_initial;
alpha(1,:,:,:)=alpha_initial;
sigmaKSI(1,:,:)=sigmaKSIinitial;  %sigma2 for ksi
%eta_temp=log(sigmaKSI(1,m,p));%rewrite to m*p



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
    sigmaE(iter)=1/gamrnd(IG_a,1/IG_b); %generate inverse gamma
    %beta(iter,:)=mvnrnd(beta_hat,sigmaE(iter)*HH_inv);
    


%step 2 sampling alpha Metropolis-Hastings sampling
   
    temp=1;
    for i=1:m
        MU_i=reshape(H_2(theta(i,:),p,q),p,p*q)*beta(iter-1,:)'; %先验
        
        sigmaKSI_i=diag(reshape(sigmaKSI(iter-1,i,:),2,1));
        for j=1:l
            alpha_t_1=reshape(alpha(iter-1,i,j,:),p,1);
            alpha_new=alpha_t_1+mvnrnd(zeros(p,1),[0.01,0;0,0.1])'; %%%步长可调节
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
prop=1/4; % proportion of gibbs samplers to use
startIndex=round((1-prop)*N);
%beta=zeros(N,p*q);   %beta=(beta_1',...beta_p')'



for index=1:p*q
    figure(5)
    subplot(p,q,index);
    plot(reshape(beta(:,index),N,1));
    %hold on
    %plot([0,N],[beta_p(index),beta_p(index)],'r--')
    figure(6)
    subplot(p,q,index);
    hist(reshape(beta(:,index),N,1));
    %hold on
    %plot([beta_p(index),beta_p(index)],[0,N],'r--')
    mean(beta(N/2:end,index));
    beta_mean(index)=mean(beta(startIndex:end,index))
    
end


%sigmaE=zeros(N,1); %sigma2 for noise variance
figure(7)
subplot(2,1,1);
plot(sigmaE(:));
%hold on
%plot([0,N],[0.14,0.14],'r--')
subplot(2,1,2);
hist(sigmaE(:));
%hold on
%plot([0.04,0.04],[0,N],'r--')
mean_sigmaE=mean(sigmaE(N/2:end));
%sigma_error_initial=sigma_error_initial
sigmaE_mean=mean(sigmaE(startIndex:end))


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
    %hold on
    %plot([0,N],[gam_p(index),gam_p(index)],'r--')
    figure(12)
    subplot(p,r,index);
    hist(reshape(gam(:,index),N,1));
    %hold on
    %plot([gam_p(index),gam_p(index)],[0,N],'r--')
    mean(gam(:,index));
    gam_mean(index)=mean(gam(startIndex:end,index))
end



%sigmaDelta=zeros(N,p);  %sigma2 for level-2 variances
figure(13)
subplot(2,1,1);
plot(sigmaDelta(:,1));
%hold on
%plot([0,N],[sigmaDelta_s(1),sigmaDelta_s(1)],'r--')
sigmaDelta_mean(1)=mean(sigmaDelta(startIndex:end,1))
subplot(2,1,2);
plot(sigmaDelta(:,2));
%hold on
%plot([0,N],[sigmaDelta_s(2),sigmaDelta_s(2)],'r--')
sigmaDelta_mean(2)=mean(sigmaDelta(startIndex:end,2))

mean(sigmaDelta(startIndex:end,1));
mean(sigmaDelta(startIndex:end,2));
%sigmaDelta_initial=sigmaDelta_initial;



figure(14)
subplot(2,1,1);
hist(sigmaDelta(startIndex:end,1),2000);
subplot(2,1,2);
hist(sigmaDelta(startIndex:end,2),2000);




