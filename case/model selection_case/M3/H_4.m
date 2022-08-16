function mat=H_4(theta,m,r)
    mat=zeros(m,r);
    for i=1:m
        mat(i,:)=h3(theta(i,:),r)';
    end
end