

%1#preprocessing\Get tpcf
%2#parameter estimation\Get initial parameter

%initial
beta1=[-3.733,-0.0005964,0.00441,-5.569e-08,-1.045e-05,1.245e-06];
beta2=[-2.061,0.003472,-0.02833,5.001e-07,6.071e-05,-9.857e-06];
gam1=[-8.641,-0.0005736,0.02684,-3.808e-07,-4.617e-05,3.2e-06];
gam2=[-3.365,0.0007443,-0.01262,7.584e-07,2.522e-05,-4.649e-06];
%3#GibbsMH
[beta,sigmaKSI,sigmaDelta,sigmaE,sigmaKSI,ksi,gam,sigmaDelta]=GibbsMH(beta1,beta2,gam1,gam2)


%GibbsHMC
%GibbsHMC()


%4#showResult
showResult(beta,sigmaKSI,sigmaDelta,sigmaE,sigmaKSI,ksi,gam,sigmaDelta)

%5#reconstruction




