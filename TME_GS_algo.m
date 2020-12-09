function [P_est,pts_proj_res,avg]=TME_GS_algo(pts,point2d,point3d,verbose)

% begin the algo
point2d_all=[point2d];
point3d_all=[point3d];

% data centering
avg=zeros(1,3);
for i=1:3
    avg(i)=mean(point3d_all(i,:));
end
point2d_cent=point2d_all-avg(1:2)'*ones(1,size(point2d_all,2));
point3d_cent=point3d_all-avg'*ones(1,size(point3d_all,2));


% Gold Standard algo
[P_est,mse] = gold_standard_algo(point2d_cent,point3d_cent,verbose);
if (verbose)
	display(mse)
end

% Project all points using the estimated transformation model
pts_cent=pts(:,1:3)'-avg(1:3)'*ones(1,size(pts,1));
pts_proj_cent=P_est*[pts_cent; ones(1,size(pts,1))]; 
pts_proj=zeros(3,size(pts,1));
pts_proj(1:2,:)=pts_proj_cent(1:2,:)+avg(1:2)'*ones(1,size(pts_proj_cent,2));
pts_proj(3,:)=pts(:,3)';

% new (19/03): adding classfication and/or intensity
N=size(pts,2);
pts_proj=zeros(N,size(pts,1));
pts_proj(3:N,:)=pts(:,3:N)';
pts_proj(1:2,:)=pts_proj_cent(1:2,:)+avg(1:2)'*ones(1,size(pts_proj_cent,2));

pts_proj_res=pts_proj';

