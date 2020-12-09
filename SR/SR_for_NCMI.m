%% 
% [M, N, ~] = size(imgRGB)
% P: a transformation model
function [int_img,z_img]=SR_for_NCMI(P,pts,avg,M,N,rx,ry,kmax)
    
    % Project all points using a given transformation model
    pts_cent=pts(:,1:3)'-avg(1:3)'*ones(1,size(pts,1));
    pts_proj_cent=P*[pts_cent; ones(1,size(pts,1))]; 
    pts_proj=zeros(size(pts,1),4);
    pts_proj(:,1:2)=(pts_proj_cent(1:2,:)+avg(1:2)'*ones(1,size(pts_proj_cent,2)))';
    pts_proj(:,3)=pts(:,3);
    pts_proj(:,4)=pts(:,4);
    
    % Create sparse images
    int_sparse=zeros(M,N); 
    z_sparse=zeros(M,N); 
    x=zeros(1,size(pts_proj,1));
    y=zeros(1,size(pts_proj,1));
    for i=1:size(pts_proj,1)
        if pts_proj(i,1)<max(rx) && pts_proj(i,1)>min(rx) &&...
           pts_proj(i,2)<max(ry) && pts_proj(i,2)>min(ry)
            x(i)=min(find(rx>=pts_proj(i,1)));
            y(i)=min(find(ry>=pts_proj(i,2)));
            int_sparse(y(i),x(i))=pts_proj(i,4); % intensity data
            z_sparse(y(i),x(i))=pts_proj(i,3); % altitude data
        end
    end
    
    %% Modify the sparse image by removing occluded points e.g. points on vertical plane (added on May 14)
    pts_idx=find(int_sparse~=0);
    [row,col]=ind2sub(size(int_sparse),pts_idx);
    
    wind_sz=5;
    diff_val=20; % value of difference (for i-image: 10 values)
    int_sparse2=zeros(size(int_sparse));
    for i=1:length(pts_idx)
        patch=int_sparse(max(1,row(i)-wind_sz):min(size(int_sparse,1),row(i)+wind_sz),...
            max(1,col(i)-wind_sz):min(size(int_sparse,2),col(i)+wind_sz));
        val=int_sparse(row(i),col(i));
        
        idx=find(patch~=0);
        if length(idx)>1
            count_near=sum(patch(idx)>val-diff_val & patch(idx)<val+diff_val); % count of pixels with z near to the center value
            count_far=sum(patch(idx)>val+diff_val | patch(idx)<val-diff_val); % count of pixels with z bigger than the center value+3 (or 20)
            if count_far>count_near
                int_sparse2(row(i),col(i))=0;
            else
                int_sparse2(row(i),col(i))=int_sparse(row(i),col(i));
            end
        else
            int_sparse2(row(i),col(i))=int_sparse(row(i),col(i));
        end
    end
    
    %
    pts_idx=find(z_sparse~=0);
    [row,col]=ind2sub(size(z_sparse),pts_idx);
    
    wind_sz=5;
    diff_val=3; % value of difference (for z-image: 3 meters)
    z_sparse2=zeros(size(z_sparse));
    for i=1:length(pts_idx)
        patch=z_sparse(max(1,row(i)-wind_sz):min(size(z_sparse,1),row(i)+wind_sz),...
            max(1,col(i)-wind_sz):min(size(z_sparse,2),col(i)+wind_sz));
        val=z_sparse(row(i),col(i));
        
        idx=find(patch~=0);
        if length(idx)>1
            count_near=sum(patch(idx)>val-diff_val & patch(idx)<val+diff_val); % count of pixels with z near to the center value
            count_far=sum(patch(idx)>val+diff_val | patch(idx)<val-diff_val); % count of pixels with z bigger than the center value+3 (or 20)
            if count_far>count_near
                z_sparse2(row(i),col(i))=0;
            else
                z_sparse2(row(i),col(i))=z_sparse(row(i),col(i));
            end
        else
            z_sparse2(row(i),col(i))=z_sparse(row(i),col(i));
        end
    end
    
    %% Run SRs
    
    % SR parameters
    % kmax=250; %max iteration number
    tol=0; %preset tolerance
    gamma=0.01;%0.02; %fixed step size
    lambda=5e-4; % soft threshold multiplication param
    l1_flag=1; % l1-norm regularization term flag

    % run SR
    [HR_int_by_grad,~]=SR_by_grad(mat2gray(int_sparse2,[0 256]),gamma,l1_flag,lambda,kmax,tol,0);
    [HR_z_by_grad,~]=SR_by_grad(z_sparse2,gamma,1,lambda,kmax,tol,0);
    
    int_img=mat2gray(flipud(HR_int_by_grad));
    z_img=mat2gray(flipud(HR_z_by_grad)*max(pts_proj(:,3)),[min(pts_proj(:,3)) max(pts_proj(:,3))]);
end