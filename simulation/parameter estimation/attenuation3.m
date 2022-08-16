

%%% theta is a vector 
clc
clear
close all

b=[0.4,0.425,0.45,0.475,0.5,0.525,0.55,0.575,0.6,0.625,0.65,0.675,0.7,0.725,0.75,0.775,0.8,0.825,0.85,0.875,0.9,0.925,0.95,0.975];

m=length(b);
p=2;
rep=20;
P=zeros(p,m,rep);
%%
%attenuation vs frequency fitting 
i=6;
figure(1);

AttenSelection=zeros(m,rep,11);

sigma_error_initial=0;
for j=1:m
    load('../attenuation/freq_chosen')
    filename = [ '../attenuation/attenuation3/attenuation' num2str(b(j)) 'thin.mat'  ] ;
    load(filename)
    x=freq_chosen(19:29);
    Atten=Atten(:,19:29);
    AttenSelection(j,:,:)=Atten;
    subplot(6,4,j);
    plot(x,Atten);
    axis([2.2,2.8,0.8,1.1])
    for k=1:rep
        P(:,j,k) = polyfit(x,Atten(k,:),p-1);
        sigma_error_initial=sigma_error_initial+sum((Atten(k,:)-(P(1,j,k)*x+P(2,j,k))).^2);
    end
end
sigma_error_initial=sigma_error_initial/(m*rep*length(x)-m*rep*p);
%%
%plot a(slope), b(intercept) vs beta(material parameter)
figure(2)
hold on
for j=1:m
    scatter(repmat(b(j),rep,1),reshape(P(1,j,:),rep,1))
end
beta1 = polyfit(b,mean(P(1,:,:),3),1);
hold on
plot(b,beta1(1)*b+beta1(2));

figure(3)
hold on
for j=1:m
    scatter(repmat(b(j),rep,1),reshape(P(2,j,:),rep,1));
end
beta2 = polyfit(b,mean(P(2,:,:),3),1);
hold on
plot(b,beta2(1)*b+beta2(2));

%%%% calculate the initial value ksi_ij
ksi_initial=zeros(m,rep,p);
sigmaKSI_initial=zeros(m,p);
for j=1:m
    for k=1:rep
        ksi_initial(j,k,1)=P(1,j,k)-(beta1(1)*b(j)+beta1(2));
        ksi_initial(j,k,2)=P(2,j,k)-(beta2(1)*b(j)+beta2(2));
    end
    sigmaKSI_initial(j,1)=var(ksi_initial(j,:,1));
    sigmaKSI_initial(j,2)=var(ksi_initial(j,:,2));
end
        

%% variance of a, b
sigmaDelta_initial=zeros(1,p);
V=zeros(m,p);
for j=1:m
    for k=1:p
        V(j,k)=var(P(k,j,:));
    end
end

%plot var(a), var(b) vs beta

figure(4)
scatter(b,reshape(log(V(:,1)),m,1))
gam1 = polyfit(b',reshape(log(V(:,1)),m,1),1);
hold on
plot(b,gam1(1)*b+gam1(2));
sigmaDelta_initial(1)=sum((reshape(log(V(:,1)),1,m)-(gam1(1)*b+gam1(2))).^2)/(m-p);


figure(5)
scatter(b,reshape(log(V(:,2)),m,1))
gam2 = polyfit(b',reshape(log(V(:,2)),m,1),1);
hold on
plot(b,gam2(1)*b+gam2(2));
sigmaDelta_initial(2)=sum((reshape(log(V(:,2)),1,m)-(gam2(1)*b+gam2(2))).^2)/(m-p);



  

