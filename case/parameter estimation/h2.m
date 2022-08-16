function vec=h2(theta, q)
    vec=zeros(q,1);
    if(q==4)
        vec=[1;theta(1);theta(2);theta(1)*theta(2)];        
    elseif(q==5)
        vec=[theta(1);theta(2);theta(1)^2;theta(2)^2;1];
    elseif(q==6)
        vec=[1;theta(1);theta(2);theta(1)^2;theta(2)^2;theta(1)*theta(2)];
    end
end