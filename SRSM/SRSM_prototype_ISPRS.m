clear
clc
% close all

format compact

verbose=0;
version='new';

% points > pts_non_gnd > pts_non_vege

%% Read LiDAR data
dir_data='/Users/thanhhuynguyen/Storage/DATA/ISPRS_Benchmark/';

iArea=1;
[~,imgRGB,ri,rx,ry]=load_datasets_isprs(iArea,1,0);

%% Parameters (listed in order of using)

ndvi_thres=0.17;
% evi2_thres=0.08;

tol=0.25; % tolerance used to determine *vegetation* lidar points

Trf=4; % elevation threshold parameter

seg_min_size=20; % building segments minimum size allowed (pixel size: 14cm)

%% NDVI-based processing: to remove vegetation points from lidar
imgD=im2double(imgRGB);

NDVI=(imgD(:,:,1)-imgD(:,:,2))./(imgD(:,:,1)+imgD(:,:,2));
[x,y]=find(NDVI>ndvi_thres);
if verbose
    figure
    imshow(imgD)
    hold on
    plot(y,x,'g.')
end

ndvi_bw=zeros(size(imgRGB,1),size(imgRGB,2));
for i=1:length(x)
    ndvi_bw(x(i),y(i))=1;
end

% figure
% imshow(ndvi_bw)

se = strel('diamond',3);
ndvi_bw2 = imclose(ndvi_bw,se);

% figure
% imshow(ndvi_bw2)

[L,n]=bwlabel(ndvi_bw2,8);

L2=zeros(size(L));
ndvi_bw3=zeros(size(L));
k=1;

for i=1:n
    [r, c] = find(L==i);
    if (numel(r)>5)
        for j=1:numel(r)
            L2(r(j),c(j))=k;
            ndvi_bw3(r(j),c(j))=1;
        end
        k=k+1;
    end
end
n2=k-1;

% figure
% imagesc(ndvi_bw3)

%% Ground filtering
% Non-ground and ground points are separated using the method called Cloth
% Simulation Filtering (CSF). It is provided as a plug-in in Cloud Compare.
%
% Reference:
% Zhang, W., Qi, J., Wan, P., Wang, H., Xie, D., Wang, X., & Yan, G. (2016).
% An Easy-to-Use Airborne LiDAR Data Filtering Method Based on Cloth Simulation.
% Remote Sensing, 8(6), 501. https://doi.org/10.3390/rs8060501

fileID = fopen(strcat(dir_data,'ALS/area1_non_gnd.txt'),'r');
formatSpec = '%f,%f,%f,%f,%f';
size_pts=[5 Inf];

tmp=fscanf(fileID, formatSpec, size_pts);
fclose(fileID);

pts_non_gnd_CSF=tmp(1:4,:)'; %[XYZI]
clear tmp

fileID = fopen(strcat(dir_data,'ALS/area1_gnd.txt'),'r');
formatSpec = '%f,%f,%f,%f,%f';
size_pts=[5 Inf];

tmp=fscanf(fileID, formatSpec, size_pts);
fclose(fileID);

pts_gnd=tmp(1:4,:)'; %[XYZI]
clear tmp


%% Elevation Thresholding
cents=[mean(pts_gnd(:,1)), mean(pts_gnd(:,2))];

x=pts_gnd(:,1)-cents(1);
y=pts_gnd(:,2)-cents(2);
z=pts_gnd(:,3);
f1 = fit([x,y],z,'poly22');

if verbose
    [X,Y] = meshgrid(linspace(min(x),max(x),100),linspace(min(y),max(y),100));
    Z=f1(X,Y);
    X=X+cents(1);
    Y=Y+cents(2);
    figure
    pcshow(pts_non_gnd_CSF(:,1:3))
    hold on
    surf(X,Y,Z)
end


x_non_gnd=pts_non_gnd_CSF(:,1)-cents(1);
y_non_gnd=pts_non_gnd_CSF(:,2)-cents(2);
T=f1(x_non_gnd,y_non_gnd)+Trf;

pts_non_gnd=pts_non_gnd_CSF(pts_non_gnd_CSF(:,3)>T(:),:);

figure
pcshow(pts_non_gnd(:,1:3))

%% Determine and remove vegetation points from lidar
li_class_pts_non_gnd=zeros(length(pts_non_gnd),1); % vege: class=1, others: class=0

for i=1:length(pts_non_gnd)
    loc=pts_non_gnd(i,1:2);
    
    idx=find(rx>loc(1)-tol & rx<loc(1)+tol);
    idy=find(ry>loc(2)-tol & ry<loc(2)+tol);
    count=sum(sum(ndvi_bw3(idy,idx)));
    if (count>3)
        li_class_pts_non_gnd(i)=1;
    end
end

figure
pcshow(pts_non_gnd(:,1:3))
hold on
pcshow(pts_non_gnd(li_class_pts_non_gnd==1,1:3),'r')

pts_vege=[]; % to stock the vegetation points
pts_non_vege=[]; % to stock the remaining non-vegetation points

% if a point is vegetation, find points in its proximity tol and remove them
for i=1:length(pts_non_gnd)
    if li_class_pts_non_gnd(i)==0
        p=pts_non_gnd(i,1:3);
        
        [~,idx] = findPoints(pts_non_gnd(li_class_pts_non_gnd==1,1:3),...
            [p(1)-tol, p(1)+tol],[p(2)-tol, p(2)+tol],0);
        if (isempty(idx))
            pts_non_vege = [pts_non_vege; pts_non_gnd(i,:)];
        else
            pts_vege = [pts_vege; pts_non_gnd(i,:)];
        end
    else
        pts_vege = [pts_vege; pts_non_gnd(i,:)];
    end
end

figure
mapshow(imgRGB,ri)
hold on
pcshow(pts_non_vege(:,1:3),'MarkerSize',30)
view(2)


%% Point cloud denoising - Use Matlab built-in pointcloud denoising function
ptCloudIn = pointCloud(pts_non_vege(:,1:3));
[~, inl] = pcdenoise(ptCloudIn);

pts_non_gnd_denoised=pts_non_vege(inl,:);

figure
pcshow(pts_non_gnd_denoised(:,1:3))
view(2)

%% Create binary grid (building pixel=1, other=0)

% without camera pose parameters (i.e. using vertical projection)

xmax=max(pts_non_gnd_denoised(:,1));
xmin=min(pts_non_gnd_denoised(:,1));
ymax=max(pts_non_gnd_denoised(:,2));
ymin=min(pts_non_gnd_denoised(:,2));

resolution = [1 1];
Ncol=round((xmax-xmin)/resolution(1));
Nrow=round((ymax-ymin)/resolution(2));
step=[(xmax-xmin)/Ncol, (ymax-ymin)/Nrow];

mask = zeros(Nrow,Ncol);
for k=1:size(pts_non_gnd_denoised,1)
    P=pts_non_gnd_denoised(k,:);
    
    col=ceil((P(1)-xmin)/step(1))+1;
    row=ceil((P(2)-ymin)/step(2))+1;
    mask(row,col)=1;
end

figure
subplot(121)
imagesc(mask)
colormap(gray)
axis equal
axis xy

xWorldLimits = [xmin xmax];
yWorldLimits = [ymin ymax];
mask_ref=maprefcells(xWorldLimits, yWorldLimits, [Nrow,Ncol]);

%% applying a opening morphological operation
se = strel('diamond',1);
mask2 = imopen(mask,se);

subplot(122)
imagesc(mask2)
colormap(gray)
axis equal
axis xy

se = strel('diamond',1);
mask3=imclose(mask2,se);

%% Labeling the building binary mask
[L,n]=bwlabel(mask3,8);
% [L,n]=bwlabel(mask2,8);

figure
imagesc(L)

L2=zeros(size(L));
n2=zeros(size(n));
k=1;

for i=1:n
    [r, c] = find(L==i);
    if (numel(r)>seg_min_size) %only take regions that have more than 20 pixels
        for j=1:numel(r)
            L2(r(j),c(j))=k;
        end
        k=k+1;
    end
end
n2=k-1;

figure
imagesc(L2)
axis xy

%% Preliminary building boundary extraction
% close all
% building_mask = mask3;
% building_label = L2;
% nb_labels = n2;
% save(strcat('building_mask_area1',version,'.mat'),'building_mask','mask_ref',...
%      'pts_non_gnd_denoised','building_label','nb_labels')
%
% [bounds,points] = buildingBoundaryExtraction(strcat('building_mask_area1',version,'.mat'),0);

points={};
for i=1:n2
    [r,c] = find(L2==i);
    x=c*resolution(1)+xmin;
    y=r*resolution(2)+ymin;
    points(i,:)={[x,y]};
end

%% Calculate z-image gradient
[pts,imgRGB,ri,rx,ry]=load_datasets_isprs(iArea,1,0);

[M,N,~]=size(imgRGB);
[HR_z_by_grad,HR1,HR0]=SR_z(pts,M,N,rx,fliplr(ry),600,1);

HR_z_by_grad_rescale=mat2gray(flipud(HR_z_by_grad),[min(pts(:,3)), max(pts(:,3))]);

% [Gx,Gy] = imgradientxy(HR_z_by_grad_rescale);

figure
imagesc(HR_z_by_grad_rescale)

% figure
% imagesc(sqrt(Gx.^2+Gy.^2),[0 0.05])
% axis equal
% axis off

% [x,y]=find(HR0~=0);
% HR0_=zeros(size(HR0));
% for i=1:length(x)
%     HR0_(x(i)-1:x(i),y(i)-1:y(i))=HR0(x(i),y(i));
% end
%
% figure
% imagesc(HR0_)
% axis equal
% myColorMap = parula(256);
% myColorMap(1,:) = 1;
% colormap(myColorMap);
% colorbar
% xlim([1 1000])
% ylim([1 1000])
% caxis([270 297])
% axis xy
% axis off
%
% figure
% imagesc(HR_z_by_grad,[270 297])
% colorbar
% xlim([1 1500])
% ylim([1 1500])
% caxis([270 297])
% axis xy
% axis off

%%
tb_idx_building=1:length(bounds);

xArea=[min(rx) max(rx)];
yArea=[min(ry) max(ry)];

% tb_idx_building_in_reg=[];
% for idx=1:length(tb_idx_building)
%     idx_building=tb_idx_building(idx);
%
%     b=bounds{idx_building};
%     if all(b(:,1)>min(xArea) & b(:,1)<max(xArea) & b(:,2)>min(yArea) & b(:,2)<max(yArea))
%         tb_idx_building_in_reg=[tb_idx_building_in_reg idx_building];
%     end
% end
tb_idx_building_in_reg=tb_idx_building;


[outline,r_outline]=geotiffread(strcat(dir_data,'/Reference_object_detection/Reference_Area/','borderline_area_1.tif'));
[r,c]=find(outline==0);
b=boundary(c,r);
[x,y]=convertMatCoorGeoCoor(r_outline.XWorldLimits, r_outline.YWorldLimits, c(b), r(b), r_outline.RasterSize);

figure
mapshow(rgb2gray(imgRGB),ri)
hold on
for idx=tb_idx_building_in_reg
    b=points{idx};
    b=b(boundary(b(:,1:2),0.75),:);
    plot(b(:,1),b(:,2),'-')
    text(mean(b(:,1)),mean(b(:,2)),num2str(idx),'Color','r')
end
plot(x,y,'y','LineWidth',2)
% xlim([min(xArea) max(xArea)])
% ylim([min(yArea) max(yArea)])

tb_idx_outside=[];
for idx=tb_idx_building_in_reg
    b=bounds{idx};
    in=inpolygon(b(:,1),b(:,2),x,y);
    if size(b,1)~=size(b(in),1)
        tb_idx_outside=[tb_idx_outside idx];
    end
end

% remove buildings outside of the borderline of the area
tb_idx_building_in_reg=setdiff(tb_idx_building_in_reg,tb_idx_outside);


figure
mapshow(rgb2gray(imgRGB),ri)
hold on
for idx=tb_idx_building_in_reg
    b=points{idx};
    b=b(boundary(b(:,1:2),0.75),:);
    plot(b(:,1),b(:,2),'-')
    text(mean(b(:,1)),mean(b(:,2)),num2str(idx),'Color','r')
end
plot(x,y,'y','LineWidth',2)
xlim([min(xArea) max(xArea)])
ylim([min(yArea) max(yArea)])

%% Snake option

SnakeOptions=struct;
SnakeOptions.Verbose=true;

SnakeOptions.Gamma=1;        % time step, default 1
SnakeOptions.Iterations=150; % number of iterations, default 100

SnakeOptions.Wline=0.04;     % attraction to line, default 0.04
SnakeOptions.Wedge=2;        % attraction to edge, default 2
SnakeOptions.Wterm=0.01;     % attraction to terminations of lines and corners, default 0.01
SnakeOptions.Sigma1=5; % sigma used to calculate image derivative, default 10
SnakeOptions.Sigma2=5; % sigma used to calculate the gradient of the edge energy image (which gives the image force), default 20
SnakeOptions.Kappa=2;   % weight of external IMAGE force, default 2

% GVF Params
SnakeOptions.Mu=0.2;     % trade between real edge vectors, noise vectors, default 0.2
SnakeOptions.Sigma3=1.0; % sigma used to calculate laplacian in GVF, default 1.0
SnakeOptions.GIterations=200; % Number of GVF iterations, default 0

SnakeOptions.Alpha=0.01; % membrane energy (1st order), default 0.2 -> amount of stretch in the snake (i.e. tension)
% SnakeOptions.Beta=0.01;  % thin-plate energy (2nd order), default 0.2 -> amount of smoothness in the snake (i.e. rigidity)
SnakeOptions.Beta0=0.01; % thin-plate energy (2nd order), default 0.2 -> amount of smoothness in the snake (i.e. rigidity)
SnakeOptions.Delta=0.1; % Balloon force, default 0.1

% flag of using beta/delta as a matrix (pixel wise beta/delta)
% SnakeOptions.betaMatFlag=0;
% SnakeOptions.deltaMatFlag=0;

% New terms
% Proposed shape sim term for snake
SnakeOptions.ShapeSimFlag=0;
SnakeOptions.ShapeSimCoeffs=[0.1 50];

% Regional Variance Energy (Kabolizade's model)
SnakeOptions.RVEFlag=0;
SnakeOptions.RVECoeff=0.5;

%%
tol=50;

%% Run snake and polygonization
% close all
SnakeOptions
snake_mask=zeros(size(imgRGB,1),size(imgRGB,2));

patch_corner=zeros(length(tb_idx_building_in_reg),4);
for idx=1:length(tb_idx_building_in_reg)
    %     idx_building=1
    idx
    idx_building=tb_idx_building_in_reg(idx);
    
    p=points{idx_building};
    p_proj=p(boundary(p(:,1:2),0.5),:);
    
    x=[min(p_proj(:,1)) max(p_proj(:,1))];
    y=[min(p_proj(:,2)) max(p_proj(:,2))];
    
    [col,row]=convertGeoCoorMatCoor(ri.XWorldLimits,ri.YWorldLimits,x,y,ri.RasterSize);
    
    patch_corner(idx,:)=[max(1,round(min(col))-tol),max(1,round(min(row))-tol),...
        min(size(imgRGB,2),round(max(col))+tol),min(size(imgRGB,1),round(max(row))+tol)];
    
    Z = HR_z_by_grad_rescale(...
        patch_corner(idx,2):patch_corner(idx,4),...
        patch_corner(idx,1):patch_corner(idx,3),:);
    
    % Get initial seeds
    [x,y]=convertGeoCoorMatCoor(ri.XWorldLimits,ri.YWorldLimits,p_proj(:,1),p_proj(:,2),ri.RasterSize);
    P=[y(:)-max(1,min(row)-tol), x(:)-max(1,min(col)-tol)];
    
    if (verbose)
        figure
        imshow(I)
        hold on
        plot(P(:,2), P(:,1), 'r.', 'MarkerSize', 15)
        % for i=1:size(P,1)
        %     text(P(i,2), P(i,1), num2str(i), 'Color','b')
        % end
    end
    SnakeOptions.nPoints=length(p_proj)*2;
    
    phi_z=HR_z_by_grad_rescale(...
        patch_corner(idx,2):patch_corner(idx,4),...
        patch_corner(idx,1):patch_corner(idx,3));
    
    bw=imbinarize(phi_z,0.78);
    SnakeOptions.deltaMat=bw*2-1;
    
    SnakeOptions.ShapeSimRef=P;
    
    %% Run the snake model
    % SnakeOptions.Verbose=1;
    [O,J]=Snake2D_pixel_wise(Z,P,SnakeOptions);
    
    snake_seeds(idx,:)={P};
    snake_results(idx,:)={O};
    
    ind=find(snake_mask(patch_corner(idx,2):patch_corner(idx,4),patch_corner(idx,1):patch_corner(idx,3))==1);
    if ~isempty(ind)
        J(ind)=1;
    end
    snake_mask(patch_corner(idx,2):patch_corner(idx,4),patch_corner(idx,1):patch_corner(idx,3))=J;
    
    %     cd ..
    %     [vertices,~,~,~]= polygonization_dutty(O,P,Vspec,Mspec,0,verbose);
    %     % lidar_mbr(idx,:)={minBoundingBox([P(:,1)'; P(:,2)'])};
    %     cd Snakes
    % close all
    
    %     building_vertices(idx,:)={vertices};
    
    if (verbose)
        figure(idx)
        imshow(I);
        hold on;
        plot([P(:,2);P(1,2)],[P(:,1);P(1,1)],'c.','LineWidth',2,'MarkerSize',10);
        plot([O(:,2);O(1,2)],[O(:,1);O(1,1)],'r.','LineWidth',2,'MarkerSize',10);
    end
end

%%

% iou_snake=sum(sum(snake_mask&ref_mask))/sum(sum(snake_mask|ref_mask))
% com_snake=sum(sum(snake_mask&ref_mask))/sum(sum(ref_mask))
% cor_snake=sum(sum(snake_mask&ref_mask))/sum(sum(snake_mask))

%% Metric result
dir2='/Users/thanhhuynguyen/Storage/DATA/';
filename=strcat(dir2,'ISPRS_Benchmark/Reference_3d_reconstruction/Reference_Buildings/');
[S,A]=shaperead(strcat(filename,'building_outline_area_1.shp'));

I=rgb2gray(imgRGB);

% Ground truth building mask
pt=zeros(numel(I),2);
for i=1:size(I,1)
    for j=1:size(I,2)
        pt((i-1)*size(I,2)+j,:)=[i,j];
    end
end

ref_mask=zeros(size(I,1),size(I,2));
for i=1:length(S)
    if (mean(S(i).X(1:end-1))<=ri.XWorldLimits(2) && mean(S(i).X(1:end-1))>=ri.XWorldLimits(1) &&...
            mean(S(i).Y(1:end-1))<=ri.YWorldLimits(2) && mean(S(i).Y(1:end-1))>=ri.YWorldLimits(1))
        
        [x,y]=convertGeoCoorMatCoor(ri.XWorldLimits,ri.YWorldLimits,...
            S(i).X(1:end-1),S(i).Y(1:end-1),ri.RasterSize);
        inBox = reshape(inpolygon(pt(:,1),pt(:,2),y,x),size(I,2),size(I,1));
        ref_mask = ref_mask | inBox';
    end
end

figure
imshow(ref_mask)

iou_snake=sum(sum(snake_mask&ref_mask))/sum(sum(snake_mask|ref_mask))
com_snake=sum(sum(snake_mask&ref_mask))/sum(sum(ref_mask))
cor_snake=sum(sum(snake_mask&ref_mask))/sum(sum(snake_mask))

acc_snake=[iou_snake com_snake cor_snake];
%

figure
mapshow((imgRGB),ri)
hold on
for i=35%1:length(S)
    plot(S(i).X,S(i).Y,'LineWidth',5,'Color','y')
    %     text(mean(S(i).X,'omitnan'),mean(S(i).Y,'omitnan'),num2str(i),'Color','g')
    %     i
    %     mean(S(i).X,'omitnan')
end

figure
imshow(snake_mask)
hold on
for i=1:length(S)
    [x,y]=convertGeoCoorMatCoor(ri.XWorldLimits,ri.YWorldLimits,...
        S(i).X(1:end-1),S(i).Y(1:end-1),ri.RasterSize);
    
    plot(x,y,'LineWidth',3)
    text(mean(x,'omitnan'),mean(y,'omitnan'),num2str(i),'Color','g')
    
end

figure
imshow(snake_mask)
hold on
for i=1:length(S)
    [x,y]=convertGeoCoorMatCoor(ri.XWorldLimits,ri.YWorldLimits,...
        S(i).X(1:end-1),S(i).Y(1:end-1),ri.RasterSize);
    
    if (polyarea(S(i).X(1:end-1),S(i).Y(1:end-1))>50)
        plot(x,y,'LineWidth',3)
        text(mean(x,'omitnan'),mean(y,'omitnan'),num2str(i),'Color','g')
    end
end

save('ISPRS_area1_result.mat','snake_mask','snake_results','acc_snake','patch_corner')



%% LiDAR boundary result
lidar_mask=zeros(size(I,1),size(I,2));
for i=1:length(snake_seeds)
    s=snake_seeds{i};
    inBox = reshape(inpolygon(pt(:,1),pt(:,2),s(:,1)+patch_corner(i,2),s(:,2)+patch_corner(i,1)),size(I,2),size(I,1));
    lidar_mask = lidar_mask | inBox';
end

figure
imshow(lidar_mask)

iou_lidar=sum(sum(lidar_mask&ref_mask))/sum(sum(lidar_mask|ref_mask))
com_lidar=sum(sum(lidar_mask&ref_mask))/sum(sum(ref_mask))
cor_lidar=sum(sum(lidar_mask&ref_mask))/sum(sum(lidar_mask))

ref_mask_color=zeros(size(I,1),size(I,2));
ref_mask_color(ref_mask==1&lidar_mask==1)=1;
ref_mask_color(ref_mask==0&lidar_mask==1)=2;
ref_mask_color(ref_mask==1&lidar_mask==0)=3;
figure
imagesc(ref_mask_color)


%%
figure
imshow(imgRGB)

for i=1:length(snake_seeds)
    s=snake_seeds{i};
    r=snake_results{i};
    %     p=building_vertices{i};
    hold on
    plot(s(:,2)+patch_corner(i,1),s(:,1)+patch_corner(i,2),'b.')
    plot(r(:,2)+patch_corner(i,1),r(:,1)+patch_corner(i,2),'g.')
    %     plot([p(2,:) p(2,1)]+patch_corner(i,1),[p(1,:) p(1,1)]+patch_corner(i,2),'r-')
end

figure
imshow(snake_mask)


%% Create shape file from snake results
create_shp_from_snake_result(version,1); % the text is added to discriminate different version

%% Measure accuracy
call_eval_snake;
