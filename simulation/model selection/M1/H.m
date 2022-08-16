function mat=H(X,theta,m,l,n,p,q)
    mat=zeros(m*l*n,p*q);
    for i=1:m
        for j=1:l
            mat((((i-1)*l+j-1)*n+1):(((i-1)*l+j)*n),:)=H_1(X,p,n)*H_2(theta(i,:),p,q);
        end
    end
end