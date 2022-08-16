function mat=H_3(theta,p,r)
    mat=kron(eye(p),h3(theta,r)');
end