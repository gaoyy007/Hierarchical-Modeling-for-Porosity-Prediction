function mat=XI(p,n,ksi,m,l,H1_mat)
%ksi m*l*p
    mat=zeros(m*l*n,1);
    for i=1:m
        for j=1:l
            mat((((i-1)*l+j-1)*n+1):(((i-1)*l+j)*n),:)=H1_mat*reshape(ksi(i,j,:),p,1);
        end
    end
end