function mat=H_1(X,n,alpha)
    mat=zeros(n,1);
    for i=1:n
        mat(i,:)=h1(X(i),alpha)';
    end
end