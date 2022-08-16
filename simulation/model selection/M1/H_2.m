function mat=H_2(theta,p,q)
    mat=kron(eye(p),h2(theta,q)');
end