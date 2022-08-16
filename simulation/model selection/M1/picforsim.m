
%pic for simulation 
prop=1/2;
startpoint=N*(1-prop);
endpoint=N;
% beta 5
% sigmaE
% gam
% sigmaDelta


meanofbeta=mean(beta(startpoint:endpoint,:));
meanofsigmaE=mean(sigmaE(startpoint:endpoint));
meanofgam=mean(gam(startpoint:endpoint,:));
meanofsigmaDelta=mean(sigmaDelta(startpoint:endpoint,:));



gapofvar
var-gapofvar

beta
sigmaE
gam
sigmaDelta