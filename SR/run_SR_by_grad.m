% Run super-resolution by minimizing SSDG to create z-image, i-image, c-image

%% Load datasets
iReg=2;
[pts,imgRGB,~,rx,ry]=load_datasets(iReg,1,0);


%% Get the transformation model
load corr4TME_reg2_byGTM
[~,pts_proj]=TME_GS_algo(pts,point2d',point3d',0);
close all

%% Run SR_by_grad with altitude data

cd SR

load HR_depth_for_grad_des

LR_img=HR_depth1;

kmax=1000;
tol=1e-5;
gamma=0.01;%0.02; %fixed step size
l1_flag=1; % l1-norm regularization term flag
lambda=5e-4; % soft threshold multiplication param

[HR_depth_by_grad,Z,~]=SR_by_grad(LR_img,gamma,l1_flag,lambda,kmax,0,1);


% save('HR_depth_by_grad_250iter.mat','HR_depth_by_grad_250iter')
figure
imagesc(HR_depth_by_grad,[0.8 1])
axis equal
axis off

% Control point evaluation code --> 'discrepancy_eval.m'

%% Run SR_by_grad with intensity data

HR_depth0=zeros(size(imgRGB,1),size(imgRGB,2)); % Sparse depth map

x=zeros(1,size(pts_proj,1));
y=zeros(1,size(pts_proj,1));
for i=1:size(pts_proj,1)
    if pts_proj(i,1)<max(rx) && pts_proj(i,1)>min(rx) && pts_proj(i,2)<max(ry) && pts_proj(i,2)>min(ry)
        x(i)=min(find(rx>=pts_proj(i,1)));
        y(i)=min(find(ry>=pts_proj(i,2)));

        HR_depth0(y(i),x(i))=pts_proj(i,4); %% intensity data
    end
end
figure
imagesc(mat2gray(HR_depth0,[0 256]))


cd SR
l1_flag=0; % l1-norm regularization term flag
[HR_int_by_grad,~]=SR_by_grad(mat2gray(HR_depth0,[0 256]),gamma,l1_flag,lambda,kmax,0,1);

[HR_int_by_grad_at,~]=SR_by_grad_aniso_tensor(mat2gray(HR_depth0,[0 256]),flipud(imgRGB),gamma,l1_flag,lambda,1000,0,1);


figure
imagesc(flipud(HR_int_by_grad_at))
axis equal
axis off
colormap(gray)

int_img=mat2gray(flipud(HR_int_by_grad),[0 0.5]);

[intensity_LR_image, R_int] = geotiffread('/Users/thanhhuynguyen/Storage/DATA/LiDAR/11_2435181F07.tif');

[col,row]=convertGeoCoorMatCoor(R_int.XWorldLimits,R_int.YWorldLimits,[rx(1), rx(end)],[ry(end), ry(1)],R_int.RasterSize);
display(['Size of LR intensity image: ' num2str(floor(row(end))-floor(row(1))) ', ' num2str(floor(col(end))-floor(col(1)))])

figure
imshow(intensity_LR_image(floor(row(1)):floor(row(end)),floor(col(1)):floor(col(end))))

figure
subplot(121), title('LR intensity image')
mapshow(intensity_LR_image, R_int)
xlim([rx(1), rx(end)])
ylim([ry(1), ry(end)])
axis off


subplot(122), title('Reconstructed intensity image')
imagesc(flipud(HR_int_by_grad),[0 0.5])
axis equal
axis off
colormap(gray)

[col, row]=convertGeoCoorMatCoor(R_int.XWorldLimits, R_int.YWorldLimits,[rx(1), rx(end)],[ry(1), ry(end)],R_int.RasterSize);

int_lr_image=intensity_LR_image(floor(row(end)):floor(row(1)),floor(col(1)):floor(col(end)),:);


figure
subplot(121),
imshow(int_lr_image)
title('LR intensity image')

subplot(122), 
imshow(int_img)
title('Reconstructed intensity image')



%% Run SR_by_grad with classification data 

HR_depth0=zeros(size(imgRGB,1),size(imgRGB,2)); % Sparse depth map

x=zeros(1,size(pts_proj,1));
y=zeros(1,size(pts_proj,1));
for i=1:size(pts_proj,1)
    if pts_proj(i,1)<max(rx) && pts_proj(i,1)>min(rx) && pts_proj(i,2)<max(ry) && pts_proj(i,2)>min(ry)
        x(i)=min(find(rx>=pts_proj(i,1)));
        y(i)=min(find(ry>=pts_proj(i,2)));

        HR_depth0(y(i),x(i))=pts_proj(i,5);
    end
end
figure
imagesc(mat2gray(HR_depth0,[0 256]))
l1_flag=0;
[HR_class_by_grad,~]=SR_by_grad(mat2gray(HR_depth0,[0 256]),gamma,l1_flag,lambda,200,0,1);

figure
imagesc(flipud(HR_class_by_grad)*5)
axis equal
axis off
title('1: Unclassified, 2: Ground, 3: Low vege., 4: Med. vege., 5: High vege.')


class=flipud(HR_class_by_grad)*5;
[x,y]=find(class>4.5);
figure
imshow(imgRGB)
hold on
plot(y,x,'b.')

%% Display sparse and dense representation

% Altitude

HR_depth1_for_disp=HR_depth1;
[x,y]=find(HR_depth1~=0);
for i=1:length(x)
HR_depth1_for_disp(x(i):x(i)+1,y(i):y(i))=HR_depth1(x(i),y(i));
end

figure
imagesc(flipud(HR_depth1_for_disp),[min(pts_proj(:,3)) max(pts_proj(:,3))])
axis equal
axis off
myColorMap = parula(256);
myColorMap(1,:) = 1;
colormap(myColorMap);
colorbar


figure
imagesc(flipud(HR_depth_by_grad.*max(pts_proj(:,3))),[min(pts_proj(:,3)) max(pts_proj(:,3))])
axis equal
axis off
myColorMap = parula(256);
myColorMap(1,:) = 1;
colormap(myColorMap);
colorbar

% Intensity
sparse=mat2gray(HR_depth0,[0 256]);

HR_depth1=mat2gray(sparse,[0 0.7])*256;
HR_depth1_for_disp=HR_depth1;
[x,y]=find(HR_depth1~=0);
for i=1:length(x)
HR_depth1_for_disp(x(i):x(i)+1,y(i):y(i))=HR_depth1(x(i),y(i));
end

figure
imagesc(flipud(HR_depth1_for_disp))
axis equal
axis off
myColorMap = jet(256);
myColorMap(1,:) = 1;
colormap(myColorMap);
colorbar

int_img=mat2gray(flipud(HR_int_by_grad),[0 0.7])*256;

figure
imagesc(int_img)
axis equal
axis off
colormap(jet)
colorbar