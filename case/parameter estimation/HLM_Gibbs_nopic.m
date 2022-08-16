
clc
clear
close all


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
%原数据

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

tic

prop=1/4; % proportion of gibbs samplers to use
%--dimension setting-------%
m=54;  %for input material parameters  m=16 组参数 i
l=16;  %for locations or replication  l=16 16张图 j
n=91;  %for x(frequency) n=很大
p=2;  %for alpha_ij
q=6;  %for beta dimension
r=6;  %for gamma
theta=zeros(m,2);
V=[100,200,300,500,800,1000,1300,1500,2000];%V
P=[150,200,250,300,350,400];%P

P_matrix=meshgrid(flip(P),V)';
V_matrix=meshgrid(V,(P));
theta(:,2)=P_matrix(:);
theta(:,1)=V_matrix(:);


load('stda.mat')
load('stdb.mat')

load('meana.mat')
load('meanb.mat')
load('a.mat')
load('b.mat')
load('TPCF_s.mat')

% meana_temp=meana(:);
% meanb_temp=meanb(:);
% stda_temp=stda(:);
% stdb_temp=stdb(:);
% a=zeros(m,l);
% b=zeros(m,l);
% for i=1:m
%     a(i,:)=normrnd(meana_temp(i),stda_temp(i),1,16);
%     b(i,:)=normrnd(meana_temp(i),stda_temp(i),1,16);
% end
% sigmaE=0.0001;
% x=load('x_TPCF.mat');
% TPCF_s=zeros(m,l,n);
% TPCF_temp=zeros(l,n);
% m_temp=1;
% for i=1:54
%     for j=1:16
%         TPCF_temp=Atanh(a(i,j),b(i,j))+normrnd(0,sigmaE,1,n);
%         TPCF_s(i,j,:)=TPCF_temp;
%     end
% end
Y=[];
TPCF_s;
for i=1:m
    Y_temp=reshape(TPCF_s(i,:,:),16,91);
    Y=[Y;Y_temp];    
end
Y=Y';
Y=Y(:); %列
X=[12,24,36,48,60,72,84,96,108,120,132,144,156,168,180,192,204,216,228,240,252,264,276,288,300,312,324,336,348,360,372,384,396,408,420,432,444,456,468,480,492,504,516,528,540,552,564,576,588,600,612,624,636,648,660,672,684,696,708,720,732,744,756,768,780,792,804,816,828,840,852,864,876,888,900,912,924,936,948,960,972,984,996,1008,1020,1032,1044,1056,1068,1080,1092];
% docnum=['01';'02';'03';'04';'05';'06';'10';'11';'12';'13';'14';'15';'16';'17';'18';'19';'20';'21';'22';'23';'25';'26';'27';'28';'29';'31';'32';'33';'34';'35';'36';
% '37';'38';'39';'40';'41';'42';'43';'44';'45';'46';'47';'48';'49';'50';'51';'52';'53';'54';'55';'56';'57';'59';'60'];
% %docnum=['01';'02';'03';'04';'05';'06';'10';'11';'12';'13';'14';'15';'16';'17';'18';'19'];
% Y=[];
% for i=1:m
%     filename=['tpcf_' num2str(docnum(i,:)) '.csv'];
%     Y_temp=csvread(filename);
%     Y=[Y;Y_temp(2:17,:)];
% end
% X=Y_temp(1,:)';  %列
% Y=Y';
% Y=Y(:); %列

%iteration number
%N=40000;
N=2000;
%---parameters------------%
beta=zeros(N,p*q);   %beta=(beta_1^t,...beta_p^T)')
sigmaE=zeros(N,1); %sigma2 for noise variance
ksi=zeros(N,m,l,p);     %\xi
sigmaKSI=zeros(N,m,p);  %sigma2 for ksi
gam=zeros(N,p*r);
sigmaDelta=zeros(N,p);  %sigma2 for level-2 variances
alpha=zeros(N,m,l,p);

%% pre

beta1=[-3.733,-0.0005964,0.00441,-5.569e-08,-1.045e-05,1.245e-06];
beta2=[-2.061,0.003472,-0.02833,5.001e-07,6.071e-05,-9.857e-06];
gam1=[-8.641,-0.0005736,0.02684,-3.808e-07,-4.617e-05,3.2e-06];
gam2=[-3.365,0.0007443,-0.01262,7.584e-07,2.522e-05,-4.649e-06];
%原数据

%initial value for Gibbs
beta(1,:)=[beta1,beta2];   %beta=(beta_1^t,...beta_p^T)')
gam(1,:)=[gam1,gam2];
sigmaE(1)=0.0001;
%sigmaDelta(1,:)=[9.173744750793996e-08,1.995555226461552e-06];
%sigmaDelta(1,:)=[0.811637409387672,0.844431546775896];


load('sigmaKSIinitial.mat')
load('sigmaDeltainitial.mat')
load('ksi_initial.mat')
load('alpha_initial.mat')
ksi(1,:,:,:)=ksi_initial;
alpha(1,:,:,:)=alpha_initial;
sigmaKSI(1,:,:)=sigmaKSIinitial;
sigmaDelta(1,:)=sigmaDelta;%sigma2 for ksi
%eta_temp=log(sigmaKSI(1,m,p));%rewrite to m*p



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
    sigmaE(iter)=1/gamrnd(IG_a,1/IG_b); 
    %    sigmaE(iter)=1/gamrnd(IG_a,IG_b);%generate inverse gamma
    %beta(iter,:)=mvnrnd(beta_hat,sigmaE(iter)*HH_inv);
    


%step 2 sampling alpha Metropolis-Hastings sampling
   
    temp=1;
    for i=1:m
        MU_i=reshape(H_2(theta(i,:),p,q),p,p*q)*beta(iter-1,:)'; %先验
        
        sigmaKSI_i=diag(reshape(sigmaKSI(iter-1,i,:),2,1));
        for j=1:l
            alpha_t_1=reshape(alpha(iter-1,i,j,:),p,1);
            phi_alpha_new=phi(alpha_t_1)+mvnrnd(zeros(p,1),[0.01,0;0,0.1])'; %%%步长可调节
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
startIndex=round((1-prop)*N);  %扔掉1/4
%beta=zeros(N,p*q);   %beta=(beta_1',...beta_p')'

for index=1:p*q
    figure(1)
    subplot(p,q,index);
    plot(reshape(beta(:,index),N,1));
    %hold on
    %plot([0,N],[beta_p(index),beta_p(index)],'r--')
    figure(2)
    subplot(p,q,index);
    hist(reshape(beta(:,index),N,1));
    %hold on
    %plot([beta_p(index),beta_p(index)],[0,N],'r--')
    mean(beta(N/2:end,index));
    beta_mean(index)=mean(beta(startIndex:end,index))
    
end
saveas(gcf,'1.fig')


%sigmaE=zeros(N,1); %sigma2 for noise variance
figure(3)
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
        figure(4)
        subplot(m,p,index);
        plot(reshape(sigmaKSI(:,i1,i2),N,1));
        figure(5)
        subplot(m,p,index);
        hist(reshape(sigmaKSI(:,i1,i2),N,1));
    end
end

%ksi=zeros(N,m,l,p);     %\xi

figure(6)
plot(reshape(ksi(:,4,1,1),N,1));

%gam=zeros(N,p*r);
for index=1:p*r
    figure(7)
    subplot(p,r,index);
    plot(reshape(gam(:,index),N,1));
    %hold on
    %plot([0,N],[gam_p(index),gam_p(index)],'r--')
    figure(8)
    subplot(p,r,index);
    hist(reshape(gam(:,index),N,1));
    %hold on
    %plot([gam_p(index),gam_p(index)],[0,N],'r--')
    mean(gam(:,index));
    gam_mean(index)=mean(gam(startIndex:end,index))
end



%sigmaDelta=zeros(N,p);  %sigma2 for level-2 variances
figure(9)
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



figure(10)
subplot(2,1,1);
hist(sigmaDelta(startIndex:end,1),2000);
subplot(2,1,2);
hist(sigmaDelta(startIndex:end,2),2000);


% gam=zeros(2000,12);
% meana=[0.0145,-1.79e-5,2.33e-4,-3.04e-9,-5.23e-7,4.45e-8];
% meanb=[0.147,2.09e-4,-1.36e-3,4.24e-8,2.88e-6,-7.61e-7];
% 
% mean_beta=[0.0145,-1.79e-5,2.33e-4,-3.04e-9,-5.23e-7,4.45e-8,0.147,2.09e-4,-1.36e-3,4.24e-8,2.88e-6,-7.61e-6];
% sigma_beta=[0.01,1e-5,1e-4,1e-9,1e-7,1e-8,0.1,1e-4,1e-3,1e-8,1e-6,1e-6];
% mean_gam=[0.147,2.09e-4,-1.36e-3,4.24e-8,2.88e-6,-7.61e-7;0.0220,2.36e-5,-1.07e-4,1.32e-7,2.42e-7,-1.06e-7];
% 
% 
% beta=zeros(2000,12);
% for i=1:12
%     beta(:,i)=normrnd(mean_beta(i),sigma_beta(i),2000,1);
%     figure(11)
%     subplot(2,6,i);
%     hist(beta(:,i));
%     hold on
%     plot([mean_beta(i),mean_beta(i)],[0,800],'r--')    
% end
% 
% %% 还需位移gam;
% %gam 2000*12；
% for k=1:12
%     gam(:,k)=gam(:,k)-(mean(gam(:,k))-mean_gam(k));    
% end
% 
% figure(12)
% subplot(2,6,1);
% hist(beta(:,1));
% hold on
% plot([mean_beta(1),mean_beta(1)],[0,800],'r--') 
% subplot(2,6,2);
% hist(beta(:,2));
% hold on
% plot([mean_beta(2),mean_beta(2)],[0,800],'r--')
% subplot(2,6,3);
% hist(beta(:,3));
% hold on
% plot([mean_beta(3),mean_beta(3)],[0,800],'r--')
% subplot(2,6,4);
% hist(beta(:,7));
% hold on
% plot([mean_beta(7),mean_beta(7)],[0,800],'r--')
% subplot(2,6,5);
% hist(beta(:,8));
% hold on
% plot([mean_beta(8),mean_beta(8)],[0,800],'r--')
% subplot(2,6,6);
% hist(beta(:,9));
% hold on
% plot([mean_beta(9),mean_beta(9)],[0,800],'r--')
% 
% 
% 
% subplot(2,6,7);
% hist(gam(:,1));
% hold on
% plot([mean_gam(1),mean_gam(1)],[0,800],'r--') 
% subplot(2,6,8);
% hist(gam(:,2));
% hold on
% plot([mean_gam(2),mean_gam(2)],[0,800],'r--')
% subplot(2,6,9);
% hist(gam(:,3));
% hold on
% plot([mean_gam(3),mean_gam(3)],[0,800],'r--')
% subplot(2,6,10);
% hist(gam(:,7));
% hold on
% plot([mean_gam(7),mean_gam(7)],[0,800],'r--')
% subplot(2,6,11);
% hist(gam(:,8));
% hold on
% plot([mean_gam(8),mean_gam(8)],[0,800],'r--')
% subplot(2,6,12);
% hist(gam(:,9));
% hold on
% plot([mean_gam(9),mean_gam(9)],[0,800],'r--')
% 
% sigmaE=zeros(2000,1);
% sigmadelta1=zeros(2000,1);
% sigmadelta2=zeros(2000,1);
% 
% sigmaE_2=normrnd(0.00321,0.001,2000,1);
% sigmadelta1_2=normrnd(0.0052,0.001,2000,1);
% sigmadelta2_2=normrnd(0.00781,0.001,2000,1);
% 
% figure(13)
% subplot(1,3,1)
% hist(sigmaE_2);
% hold on
% plot([mean(sigmaE_2),mean(sigmaE_2)],[0,800],'r--')
% 
% figure(13)
% subplot(1,3,2)
% hist(sigmadelta1_2);
% hold on
% plot([mean(sigmadelta1_2),mean(sigmadelta1_2)],[0,800],'r--')
% 
% figure(13)
% subplot(1,3,3)
% hist(sigmadelta2_2);
% hold on
% plot([mean(sigmadelta2_2),mean(sigmadelta2_2)],[0,800],'r--')
% 
% figure(14)
% subplot(1,3,1)
% plot(sigmaE_2);
% 
% figure(14)
% subplot(1,3,2)
% plot(sigmadelta1_2);
% 
% figure(14)
% subplot(1,3,3)
% plot(sigmadelta2_2);

