m=54;
l=16;
ksia;
ksib;
ksi_initial=zeros(m,l,2);
ksi_initial(:,:,1)=ksia;
ksi_initial(:,:,2)=ksib;

sigmaKSIinitial=zeros(m,2);
sigmaKSIinitial(:,1)=sigmaa;
sigmaKSIinitial(:,2)=sigmab;

alpha_initial=zeros(m,m,2);

docnum=['01';'02';'03';'04';'05';'06';'10';'11';'12';'13';'14';'15';'16';'17';'18';'19';'20';'21';'22';'23';'25';'26';'27';'28';'29';'31';'32';'33';'34';'35';'36';
'37';'38';'39';'40';'41';'42';'43';'44';'45';'46';'47';'48';'49';'50';'51';'52';'53';'54';'55';'56';'57';'59';'60'];
for i=1:m    
    filename=['para_' num2str(docnum(i,:)) '.csv'];
    para_temp=csvread(filename);
    alpha_initial(i,:,:)=para_temp;    
end
