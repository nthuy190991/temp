function res=calc_mi_with_theta_6params(x0)
load('temp_6params_9patches.mat','opt_img','pts_small','K','row_col')

R=rotz(x0(1,3))*roty(x0(1,2))*rotx(x0(1,1)); %first line of x0: yaw, pitch, roll
C=x0(2,:)'; %second line of x0: X, Y, Z from translation

P=K*R*[eye(3) -C(1:3)];

% disp(P)

XYZ=pts_small(:,1:3)';
n=size(XYZ,2);
xy_est=P*[XYZ; ones(1,n)];

for i = 1:n
    xy_est(:,i) = xy_est(:,i)/xy_est(end,i);
end
pts_proj=xy_est';
pts_proj(:,4)=pts_small(:,4);

[M,N]=size(opt_img);

% Create sparse image
HR0=zeros(M,N);
x=zeros(1,size(pts_proj,1));
y=zeros(1,size(pts_proj,1));
for i=1:size(pts_proj,1)
    if pts_proj(i,1)-row_col(2)<N && pts_proj(i,1)-row_col(2)>0 &&...
            pts_proj(i,2)-row_col(1)<M && pts_proj(i,2)-row_col(1)>0
        x(i)=ceil(pts_proj(i,1)-row_col(2));
        y(i)=ceil(pts_proj(i,2)-row_col(1));
        HR0(y(i),x(i))=pts_proj(i,4); % intensity data
    end
end

%     figure
%     imagesc(HR0,[0 256])

% if ~any(HR0==0)
% SR parameters
kmax=500; %max iteration number
tol=0; %preset tolerance
gamma=0.01;%0.02; %fixed step size
lambda=5e-4; % soft threshold multiplication param
l1_flag=1; % l1-norm regularization term flag

% run SR
[HR_int_by_grad,~]=SR_by_grad(mat2gray(HR0,[0 256]),gamma,l1_flag,lambda,kmax,tol,0);
int_img=mat2gray((HR_int_by_grad),[0 0.5]);

res=-mi(opt_img,int_img);
%     figure
%     subplot(121), imshow(opt_img)
%     subplot(122), imshow(int_img)

% else
%     res=0;
% end
display(['mi=' num2str(-res)])
end
