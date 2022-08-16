function vec=h3(theta, r)
    vec=zeros(r,1);
    if(r==4)
        vec=[1;theta(1);theta(2);theta(1)*theta(2)];  
    elseif(r==5)
        vec=[theta(1);theta(2);theta(1)^2;theta(2)^2;1];
    elseif(r==6)
        vec=[1;theta(1);theta(2);theta(1)^2;theta(2)^2;theta(1)*theta(2)];
    end
end