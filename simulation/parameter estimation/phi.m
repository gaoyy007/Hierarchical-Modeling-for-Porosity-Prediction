function X=phi(X)
X(1,:)=log(X(1,:));
X(2,:)=log(X(2,:)./(1-X(2,:)));
end