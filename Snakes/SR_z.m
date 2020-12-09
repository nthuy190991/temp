%% 
% [M, N, ~] = size(imgRGB)
function [HR_z_by_grad,HR1,HR0]=SR_z(pts,M,N,rx,ry,kmax,l1_flag)
    
%     % Project all points using a given transformation model
%     pts_cent=pts(:,1:3)'-avg(1:3)'*ones(1,size(pts,1));
%     pts_proj_cent=P*[pts_cent; ones(1,size(pts,1))]; 
%     pts_proj=zeros(size(pts,1),4);
%     pts_proj(:,1:2)=(pts_proj_cent(1:2,:)+avg(1:2)'*ones(1,size(pts_proj_cent,2)))';
%     pts_proj(:,3)=pts(:,3);
%     pts_proj(:,4)=pts(:,4);

    pts_proj=pts(:,1:3);

    % Create sparse image
    HR0=zeros(M,N); 
    x=zeros(1,size(pts_proj,1));
    y=zeros(1,size(pts_proj,1));
    for i=1:size(pts_proj,1)
        if pts_proj(i,1)<max(rx) && pts_proj(i,1)>min(rx) &&...
           pts_proj(i,2)<max(ry) && pts_proj(i,2)>min(ry)
            x(i)=min(find(rx>=pts_proj(i,1)));
            y(i)=min(find(ry>=pts_proj(i,2)));
            HR0(y(i),x(i))=pts_proj(i,3); % altitude data
        end
    end
    
    % Modify the sparse image by removing occluded points e.g. points on vertical plane (added on May 14)
    pts_idx=find(HR0~=0);
    [row,col]=ind2sub(size(HR0),pts_idx);
    
    wind_sz=5;
    diff_val=3; % value of difference (for z-image: 3 meters)
    HR1=zeros(size(HR0));
    for i=1:length(pts_idx)
        patch=HR0(max(1,row(i)-wind_sz):min(size(HR0,1),row(i)+wind_sz),...
            max(1,col(i)-wind_sz):min(size(HR0,2),col(i)+wind_sz));
        val=HR0(row(i),col(i));
        
        idx=find(patch~=0);
        if length(idx)>1
            count_near=sum(patch(idx)>val-diff_val & patch(idx)<val+diff_val); % count of pixels with z near to the center value
            count_far=sum(patch(idx)>val+diff_val | patch(idx)<val-diff_val); % count of pixels with z bigger than the center value+3 (or 20)
            if count_far>count_near
                HR1(row(i),col(i))=0;
            else
                HR1(row(i),col(i))=HR0(row(i),col(i));
            end
        else
            HR1(row(i),col(i))=HR0(row(i),col(i));
        end
    end
    
    % SR parameters
    % kmax=250; %max iteration number
    tol=0; %preset tolerance
    gamma=0.05;%0.02; %fixed step size
    lambda=5e-4; % soft threshold multiplication param
%     l1_flag=1; % l1-norm regularization term flag

%     HR_z_by_grad=[];
    % run SR
    d=pwd;
%     cd('/Users/thanhhuynguyen/Dropbox/Thesis/dev/SR')
    cd('../SR')
    [~,HR_z_by_grad]=SR_by_grad(HR1,gamma,l1_flag,lambda,kmax,tol,1);
    %% 
    cd(d)
end