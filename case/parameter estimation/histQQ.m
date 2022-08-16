mean_a_QQ=mean(erroramean);
std_a_QQ=std(erroramean);


%sample_a_QQ=normrnd(mean_a_QQ,std_a_QQ,[4,48]);
sample_a_QQ=normrnd(0,0.9,[4,48]);
csvwrite('sample_a_QQ.csv',sample_a_QQ)
mean_a_QQ,std_a_QQ

histogram(sample_a_QQ(1,:))
