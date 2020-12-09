% Depth map super resolution
% Input: 
% - a low-resolution 3D point cloud
% - a high-resolution image, giving size and coordinates, and for comparison
% - a transformation model from 3D point cloud -> image frame of ref
%
% Output:
%   an image whose size=size(HR image), with pixels representing altitude,
%   interpreted/interpolated from 3D poinr cloud


%% Load datasets
% Load the optical image
[A, R] = geotiffread('/Users/thanhhuynguyen/Storage/DATA/Orthos_2016/2016_CMQ_243500-5181500.tif');
figure
mapshow(A,R)
hold on

% Load LiDAR point cloud
fileID = fopen('/Users/thanhhuynguyen/Storage/DATA/LiDAR/nineth2p_txyzic.txt','r');
formatSpec = '%f,%f,%f,%f,%d,%d';
sizeA=[6 Inf];

tmp=fscanf(fileID, formatSpec, sizeA);
fclose(fileID);

ptCloudDBMat = tmp(2:5,:)'; %[XYZI]
ptCloudDBClass = tmp(6,:)'; %[C]

figure
pcshow(ptCloudDBMat(:,1:3),ptCloudDBClass,'MarkerSize',20)
view(2)
title('1: Unclassified, 2: Ground, 3: Low vege., 4: Med. vege., 5: High vege.')

% figure
% pcshow(ptCloudDBMat(ptCloudDBClass==3|ptCloudDBClass==4|ptCloudDBClass==5,1:3),...
%     'MarkerSize',20)
% view(2)

points=ptCloudDBMat(:,1:3);

%% Get datasets on the considered region, i.e. reg2

% 1st row: residential, 2nd: avg-sized buildings, 3rd: big and complex buildings
row=[1 2000; 1 1550; 3300 5700];
col=[1 2000; 2350 4000; 1500 3500];
iReg=2;

imgRGB = A(row(iReg,1):row(iReg,2), col(iReg,1):col(iReg,2), :);
rx=R.XWorldLimits(1)+(col(iReg,1):col(iReg,2))*R.CellExtentInWorldX;
ry=R.YWorldLimits(1)+(size(A,1)-fliplr(row(iReg,1):row(iReg,2)))*R.CellExtentInWorldY; %note: directions of row and y are opposite

%small raster
ri=R;
ri.XWorldLimits = [rx(1), rx(end)];
ri.YWorldLimits = [ry(1), ry(end)];
ri.RasterSize = size(imgRGB);

figure
    mapshow(rgb2gray(A),R)
    hold on
    mapshow(imgRGB,ri)

[pts, ~]=findPoints(points,[rx(1), rx(end)],[ry(1), ry(end)],0);
figure
    mapshow(imgRGB,ri)
    hold on
    pcshow(pts(:,1:3),'MarkerSize',20)
    view(2)

%% Get the transformation model
load corr4TME_reg2_byGTM
[P_est,pts_proj]=TME_GS_algo_IGARSS(imgRGB,ri,pts,point2d',point3d');
close all

figure
mapshow(imgRGB, ri)
hold on
scatter(pts_proj(:,1),pts_proj(:,2),10,pts_proj(:,3),'filled')




%% Create the HR depth map (getting the pixel value of nearby pixels)
HR_depth=zeros(size(imgRGB,1),size(imgRGB,2));
mask=zeros(size(imgRGB,1),size(imgRGB,2));

x=zeros(1,size(pts_proj,1));
y=zeros(1,size(pts_proj,1));
for i=1:size(pts_proj,1)
    if pts_proj(i,1)<max(rx) && pts_proj(i,1)>min(rx) && pts_proj(i,2)<max(ry) && pts_proj(i,2)>min(ry)
        x(i)=min(find(rx>=pts_proj(i,1)));
        y(i)=min(find(ry>=pts_proj(i,2)));

        HR_depth(y(i),x(i))=pts_proj(i,3);
        mask(y(i),x(i))=i;
    end
end

figure
imagesc(flipud(HR_depth),[min(pts_proj(:,3)) max(pts_proj(:,3))]), axis equal, axis off

% [X,Y]=meshgrid(rx,ry);
% Vq = interp2(X,Y,V,Xq,Yq,'cubic');

for i=1:size(imgRGB,1)
    i
    for j=1:size(imgRGB,2)
        if mask(i,j)==0
            window=5;
            i0=max(1,i-window);
            j0=max(1,j-window);
            imax=min(size(imgRGB,1),i+window);
            jmax=min(size(imgRGB,2),j+window);
            idx=find(mask(i0:imax,j0:jmax)~=0);
            if ~isempty(idx)
                % get avalable nearest non-zero pixels
                dist=zeros(1,length(idx));
                for k=1:length(idx)
                    [a,b]=ind2sub(size(mask(i0:imax,j0:jmax)),idx(k));
                    dist(k)=norm([a+i0-1 b+j0-1]-[i j]);
                end 
                [min_dists,idx2]=mink(dist,3); % get the 3 nearest pixel
                
                if ~isempty(idx2)
                    patch=HR_depth(i0:imax,j0:jmax);
                    [a,b]=ind2sub(size(patch),idx(idx2));
                    w=min_dists./sum(min_dists);
%                     HR_depth(i,j)=mean(diag(patch(a,b))); % mean of altitude values
                    HR_depth(i,j)=w*(diag(patch(a,b))); % weighted sum of altitude values
                end
            end
        end
    end
end

HR_depth_8bit=flipud(mat2gray(HR_depth,[min(pts_proj(:,3)) max(pts_proj(:,3))]));


figure
subplot(121)
imagesc(flipud(HR_depth),[min(pts_proj(:,3)) max(pts_proj(:,3))]), axis equal, axis off
subplot(122)
imshow(imgRGB)

figure
imshowpair(HR_depth_8bit,rgb2gray(imgRGB),'checkerboard')

figure
imshowpair(HR_depth_8bit,rgb2gray(imgRGB),'montage')

%%%%%%%%%%%%%%%%%%%
%% March 04-08 2019
%%%%%%%%%%%%%%%%%%%
%% Create HR depth map, by propagate pixel values dictated by gradient
HR_depth_by_grad=zeros(size(imgRGB,1),size(imgRGB,2));
HR_depth0=zeros(size(imgRGB,1),size(imgRGB,2)); % Sparse depth map

x=zeros(1,size(pts_proj,1));
y=zeros(1,size(pts_proj,1));
for i=1:size(pts_proj,1)
    if pts_proj(i,1)<max(rx) && pts_proj(i,1)>min(rx) && pts_proj(i,2)<max(ry) && pts_proj(i,2)>min(ry)
        x(i)=min(find(rx>=pts_proj(i,1)));
        y(i)=min(find(ry>=pts_proj(i,2)));

        HR_depth0(y(i),x(i))=pts_proj(i,3);
        mask(y(i),x(i))=i;
    end
end

norm(HR_depth0)
% ans =
%    1.4372e+04

norm(HR_depth0)*1e-5
% ans =
%     0.1437


%% Remove points that are projected, but should be occluded, e.g. points on vertical plane

pts_idx=find(HR_depth0~=0);
[row,col]=ind2sub(size(HR_depth0),pts_idx);

point_spacing=0.5;
wind_sz=5;
HR_depth1=zeros(size(HR_depth0));
for i=1:length(pts_idx)
    patch=HR_depth0(max(1,row(i)-wind_sz):min(size(HR_depth0,1),row(i)+wind_sz),...
        max(1,col(i)-wind_sz):min(size(HR_depth0,2),col(i)+wind_sz));
    val=HR_depth0(row(i),col(i));
    
    idx=find(patch~=0);
    if length(idx)>=1
        count_near=sum(patch(idx)>val-3 & patch(idx)<val+3); % count of pixels with z near to the center value
        count_far=sum(patch(idx)>val+3); % count of pixels with z bigger than the center value+3
        if count_far>(count_near)
            HR_depth1(row(i),col(i))=0;
        else
            HR_depth1(row(i),col(i))=HR_depth0(row(i),col(i));
        end
    else
        HR_depth1(row(i),col(i))=HR_depth0(row(i),col(i));
    end
end

figure
imagesc(flipud(HR_depth1),[min(pts_proj(:,3)) max(pts_proj(:,3))])
axis equal
axis off

save('HR_depth_for_grad_des.mat','HR_depth1','HR_depth0')

%%
HR_depth1_for_disp=HR_depth1;
[x,y]=find(HR_depth1~=0);
for i=1:length(x)
    HR_depth1_for_disp(x(i)-1:x(i)+1,y(i))=HR_depth1(x(i),y(i));
end

figure
imagesc(flipud(HR_depth1_for_disp),[min(pts_proj(:,3)) max(pts_proj(:,3))])
axis equal
axis off
myColorMap = parula(256);
myColorMap(1,:) = 1;
colormap(myColorMap);
colorbar

%% Run descent gradient, minimization of horizontal and vertical gradients and l1-norm regularization
load HR_depth_for_grad_des

l1_flag=0; % l1-norm regularization term flag
lambda=5e-4; % soft threshold multiplication param

LR_img=HR_depth1;%(1:1000,1:1000);
Y=reshape(mat2gray(LR_img),[],1);
Ny=size(LR_img,1);
Nx=size(LR_img,2);

% test=[0 0 1 0; 0 0 0 0; 0 1 0 0; 0 0 0 0];
% Y=reshape(test,[],1);
% Nx=size(test,1); 
% Ny=size(test,2);

Idx=find(Y==0);
kmax=500;
tol=1e-5;

k=1; 
t_old=1;
y_old=Y;
x_old=Y;%zeros(size(Y));
x_new=zeros(size(Y));
y_new=Y;
gamma=0.01;%0.02; %fixed step size
while k<=kmax
    if ~l1_flag
        [delta_f,~,~]=gradient_Dxy(y_old,Ny);
        x_new(Idx)=y_old(Idx) - gamma.*delta_f(Idx);
    else
        [delta_f,~,~]=gradient_Dxy(y_old,Ny);
        x_new(Idx)=soft_thres(y_old(Idx) - gamma.*delta_f(Idx), gamma*lambda);
    end
    t_new=0.5*(1+sqrt(1+4*t_old^2));
    
    y_new(Idx)=x_new(Idx)+(t_old-1)/t_new*(x_new(Idx)-x_old(Idx));
    
    if ~mod(k,10)
        display(['iter=' num2str(k) ', |y_new-y_old|=' num2str(norm(y_new-y_old))])
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

figure
imagesc(HR_depth_by_grad)
axis equal
axis off

save('HR_depth_for_grad_des_no_l1.mat','HR_depth1','HR_depth_by_grad')

figure
imagesc(LR_img)

Y_res_rescaled=(max(LR_img(:))-min(LR_img(:))).*Y_res + min(LR_img(:));
HR_depth_by_grad_rescaled=reshape(Y_res_rescaled,Ny,Nx);

figure
imagesc(flipud(HR_depth_by_grad_rescaled),[min(pts_proj(:,3)) max(pts_proj(:,3))])
axis equal
axis off

figure
imshowpair(rgb2gray(imgRGB),flipud(HR_depth_by_grad_rescaled),'checkerboard')

figure
imshowpair((imgRGB),flipud(HR_depth_by_grad_rescaled),'checkerboard')

Z=flipud(HR_depth_by_grad_rescaled);
x=1:500;
y=1:500;
figure
surf(Z(x,y),imgRGB(x,y,:))
% zlim([100 125])
axis equal


%% colorize the HR depth map
I1= flipud(HR_depth_by_grad_rescaled);
I2 = mat2gray(I1,[min(pts_proj(:,3)) max(pts_proj(:,3))]);
figure
imshow(I2)

HRD_rgb(:,:,1)=I2;
HRD_rgb(:,:,2)=I2;
HRD_rgb(:,:,3)=I2;

figure
imshow(HRD_rgb)

figure
imshowpair(imgRGB,I2,'checkerboard')

figure
imagesc(HRD_rgb)

%% Measure SSIM
load HR_depth_for_grad_des_no_l1

Y_res_rescaled=(max(HR_depth1(:))-min(HR_depth1(:))).*HR_depth_by_grad(:) + min(HR_depth1(:));
HR_depth_by_grad_rescaled=reshape(Y_res_rescaled,size(HR_depth_by_grad,1),size(HR_depth_by_grad,2));

figure
imagesc(flipud(HR_depth_by_grad_rescaled),[min(pts_proj(:,3)) max(pts_proj(:,3))])
axis equal
axis off



Z=mat2gray(flipud(HR_depth_by_grad_rescaled),[min(pts_proj(:,3)) max(pts_proj(:,3))]);
I=mat2gray(rgb2gray(imgRGB));


[ssimval, ssimmap] = ssim(Z,I);
  
fprintf('The SSIM value is %0.4f.\n',ssimval);

figure
subplot(121), imshow(Z)
subplot(122), imshow(I)

figure, imshow(ssimmap,[]);
title(sprintf('ssim Index Map - Mean ssim Value is %0.4f',ssimval));

[peaksnr,snr] = psnr(Z,I)

mse=immse(Z,I)

[peaksnr,snr] = psnr(A,B)
