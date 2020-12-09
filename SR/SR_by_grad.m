function [HR_depth_by_grad,HR_depth_by_grad_rescaled,calc_f]=SR_by_grad(sparse_img,gamma,l1_flag,lambda,kmax,tol,verbose)
% Y=reshape(mat2gray(sparse_img),[],1);
Y=reshape(mat2gray(sparse_img,[0 max(sparse_img(:))]),[],1); % modified on Feb 18 2020 (for case of negative-elevation points)
Ny=size(sparse_img,1);
Nx=size(sparse_img,2);

Idx=find(Y==0);
% kmax=500;
% tol=1e-5;

k=1; 
t_old=1;
y_old=Y;
x_old=Y;%zeros(size(Y));
x_new=zeros(size(Y));
y_new=Y;
% gamma=0.01;%0.02; %fixed step size

if tol==0
    tol=1e-5*norm(Y)
end
while k<=kmax
%     [delta_f,~,~]=gradient_Dxy_fast(y_old,Ny);
    if ~l1_flag
        [delta_f,~,~]=gradient_Dxy(y_old,Ny);
        x_new(Idx)=y_old(Idx) - gamma.*delta_f(Idx);
    else
        [delta_f,~,~]=gradient_Dxy(y_old,Ny);
        x_new(Idx)=soft_thres(y_old(Idx) - gamma.*delta_f(Idx), gamma*lambda);
    end
    t_new=0.5*(1+sqrt(1+4*t_old^2));
    
    y_new(Idx)=x_new(Idx)+(t_old-1)/t_new*(x_new(Idx)-x_old(Idx));
    
%     Y_res=Y;
%     Y_res(Idx)=y_new(Idx);  
%     demo(k,:)={Y_res};
    
    if verbose
        if ~mod(k,50)
            display(['iter=' num2str(k) ', |y_new-y_old|=' num2str(norm(y_new-y_old))])
        end

        if any(y_old==0) && all(y_new~=0)
            display(['Last iteration that has 0 value: ' num2str(k-1)])
        end
    end
    
    % Stopping criterion
    if norm(y_new-y_old)<tol
        break
    end
    
    % Update
    t_old=t_new;
    y_old(Idx)=y_new(Idx);
    x_old(Idx)=x_new(Idx);
    k=k+1;
end
Y_res=Y;
Y_res(Idx)=y_new(Idx);

calc_f=cost_func_Dxy(Y_res,Ny);

HR_depth_by_grad=reshape(Y_res,Ny,Nx);

Y_res_rescaled=(max(sparse_img(:))-min(sparse_img(:))).*Y_res + min(sparse_img(:));
HR_depth_by_grad_rescaled=reshape(Y_res_rescaled,Ny,Nx);


if verbose
    figure
    imagesc(HR_depth_by_grad)
    axis equal
    axis off
end
% save('/Users/thanhhuynguyen/Storage/demo_SR_20201126_int.mat','demo','-v7.3')
end