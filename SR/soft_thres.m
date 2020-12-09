function x=soft_thres(u,lambda)
    sign_u=u./abs(u+eps);   
    x=sign_u.*max(0,abs(u)-lambda);
end