function X=Iphi(X)
X(1,:)=exp(X(1,:));
X(2,:)=exp(X(2,:))./(1+exp(X(2,:)));
end