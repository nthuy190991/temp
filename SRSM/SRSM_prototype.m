clear
clc
% close all

format compact

verbose=0;

% Notes:
% Even though the LiDAR 2017 Ville de Québec has point classification,
% including buildings and ground, the following script performs
% the proposed SRSM without such information.
%
% In the scope of running this algorithm in many tiles, I recommend the use
% of point classification to reduce the computing time.
%
% Since the LiDAR point cloud has very good vegetation accuracy, the optical
% image is not needed. Vegetation points are excluded in the Preliminary process
% which gives the snake initial points, using the point classification.

flag_use_classif=0;

%% Read LiDAR data
filename = '17_4381_dc.txt';
txt_filename=strcat(filename(1:end-3),'txt');

if exist(strcat(dir_data,txt_filename),'file')~=2
    las2txt_matlab(filename, '', ' --parse txyzic');
end

fileID = fopen(strcat(dir_data,txt_filename),'r'); % Field order: [t, x, y, z, i, c]
formatSpec = '%f,%f,%f,%f,%d,%d'; % comma-separated
size_pts=[6 Inf];

tmp=fscanf(fileID, formatSpec, size_pts);
fclose(fileID);

points=tmp(2:6,:)'; %[XYZIC]
clear tmp

if verbose
    figure
    pcshow(points(points(:,5)==3|points(:,5)==4|points(:,5)==5,1:3))
end

%% Parameters (listed in order of using)
% tol=0.1; % tolerance used to determine *vegetation* lidar points
% Trf = 1.5; % height threshold parameter

seg_min_size=18; % building segments minimum size allowed (equivalent 10m2)

%% Remove points that are close to a tree

% Non-ground and ground points are separated using the method called Cloth
% Simulation Filtering. It is provided as a plug-in in Cloud Compare.
%
% Ref:
% Zhang, W., Qi, J., Wan, P., Wang, H., Xie, D., Wang, X., & Yan, G. (2016).
% An Easy-to-Use Airborne LiDAR Data Filtering Method Based on Cloth Simulation.
% Remote Sensing, 8(6), 501. https://doi.org/10.3390/rs8060501

fileID = fopen('4381_non_gnd.txt','r'); % Field order: [x, y, z, t, i, c]
formatSpec = '%f %f %f %f %d %d'; % space-separated
size_pts=[6 Inf];

tmp=fscanf(fileID, formatSpec, size_pts);
fclose(fileID);

pts_non_gnd=tmp([1:3, 5, 6],:)'; %[XYZIC]
clear tmp

points_bd=pts_non_vege(idx_bd,:);


%%

if flag_use_classif
    li_class=points(:,5)==3|points(:,5)==4|points(:,5)==5;
    
    pts_non_vege=points(~li_class,1:3);
    
    pts_non_gnd_denoised=points(points(:,5)==6,1:3);
else
    li_class=points(:,5)==3|points(:,5)==4|points(:,5)==5;
    pts_non_gnd_denoised=[];
end

if verbose
    figure
    pcshow(points(~li_class,1:3))
end

%% Create building mask (binary image)
str=strsplit(filename,'17_');
tile_no=str{2}(1:4);

strx=tile_no(1:2);
stry=tile_no(3:4);

xmin=2e5+str2double(strx)*1e3;
xmax=2e5+str2double(strx)*1e3+1e3;
if str2double(stry)>74
    ymin=51e5+str2double(stry)*1e3;
    ymax=51e5+str2double(stry)*1e3+1e3;
elseif str2double(stry)<6
    ymin=52e5+str2double(stry)*1e3;
    ymax=52e5+str2double(stry)*1e3+1e3;
end

[pts_non_gnd_denoised,~]=findPoints(pts_non_gnd_denoised,[xmin,xmax],[ymin,ymax],0);

resolution = [0.75 0.75];
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

if verbose
    figure
    subplot(121)
    imagesc(mask)
    colormap(gray)
    axis equal
    axis xy
end

xWorldLimits = [xmin xmax];
yWorldLimits = [ymin ymax];
mask_ref=maprefcells(xWorldLimits, yWorldLimits, [Nrow,Ncol]);

%% applying a opening morphological operation
se = strel('diamond',1);
% se = strel('square',2);
mask2 = imopen(mask,se);

if verbose
    subplot(122)
    imagesc(mask2)
    colormap(gray)
    axis equal
    axis xy
end

se = strel('diamond',1);
mask3=imclose(mask2,se);

%% Labeling the building binary mask
[L,n]=bwlabel(mask3,8);
% [L,n]=bwlabel(mask2,8);

if verbose
    figure
    imagesc(L)
end

L2=zeros(size(L));
n2=zeros(size(n));
k=1;

for i=1:n
    [r, c] = find(L==i);
    if (numel(r)>=seg_min_size)
        %only take regions that have more than 20 pixels
        for j=1:numel(r)
            L2(r(j),c(j))=k;
        end
        k=k+1;
    end
end
n2=k-1;

if verbose
    figure
    imagesc(L2)
    axis xy
end

%% building boundary extraction
bounds={};
for i=1:n2
    [r,c] = find(L2==i);
    b=boundary(c,r,1);
    
    [x,y]=convertMatCoorGeoCoor(mask_ref.XWorldLimits,mask_ref.YWorldLimits,c(b),Nrow-r(b),mask_ref.RasterSize);
    bounds(i,:)={[x,y]};
end

bounds_no_nan={};
for i=1:length(bounds)
    b=bounds{i,:};
    if all(all(~isnan(b)))
        bounds_no_nan(i,:)={b};
    else
        [x,y]=find(isnan(b));
        b_=b;
        b_(x,:)=[];
        bounds_no_nan(i,:)={b_};
    end
end

idx_all_nan=[];
for i=1:length(bounds_no_nan)
    b=bounds_no_nan{i,:};
    if size(b,1)==0
        idx_all_nan=[idx_all_nan i];
    end
end
bounds_no_nan(idx_all_nan)=[];

if verbose
    figure
    %     mapshow(imgRGB,ri)
    hold on
    for i=1:length(bounds)
        b=bounds{i,:};
        if all(all(~isnan(b)))
            plot(b(:,1),b(:,2))
        else
            [x,y]=find(isnan(b));
            b_=b;
            b_(x,:)=[];
            plot(b_(:,1),b_(:,2),'r')
        end
    end
end

%% Load LiDAR point cloud
clear li_class pts_non_gnd_denoised mask mask2 mask3

if verbose
    figure
    pcshow(points(points(:,5)==6,1:3))
    
    figure
    pcshow(points(1:100:end,1:3))
end

% pts=points(:,1:3);

%% z-image metadata (size and resolution)
N=floor((xmax-xmin)/0.15); % resolution: 15cm
M=floor((ymax-ymin)/0.15);

rx=linspace(xmin,xmax,N);
ry=linspace(ymin,ymax,M);
ri=maprefcells([xmin,xmax],[ymin,ymax],[M,N]);
%% Calculate z-image gradient

[HR_z_by_grad,HR1,HR0]=SR_z(points,M,N,rx,ry,300,1);

HR_z_by_grad_rescale=mat2gray(flipud(HR_z_by_grad),[min(points(:,3)), max(points(:,3))]);

if verbose
    figure
    imagesc(HR_z_by_grad_rescale)
end

%%
tb_idx_building_in_reg=1:length(bounds_no_nan);

if verbose
    figure
    for i=tb_idx_building_in_reg
        b=bounds_no_nan{i,:};
        plot(b(:,1),b(:,2))
        text(mean(b(:,1)),mean(b(:,2)),num2str(i),'Color','g','FontSize',12)
    end
end

%% Snake option
SnakeOptions=struct;
SnakeOptions.Verbose=verbose;

SnakeOptions.Gamma=1;        % time step, default 1
SnakeOptions.Iterations=150; % number of iterations, default 100

SnakeOptions.Wline=0.04;    % attraction to line, default 0.04
SnakeOptions.Wedge=2;       % attraction to edge, default 2
SnakeOptions.Wterm=0.01;    % attraction to terminations of lines and corners, default 0.01
SnakeOptions.Sigma1=5;      % sigma used to calculate image derivative, default 10
SnakeOptions.Sigma2=5;      % sigma used to calculate the gradient of the edge energy image (which gives the image force), default 20
SnakeOptions.Kappa=4;       % weight of external IMAGE force, default 2

% GVF Params
SnakeOptions.Mu=0.2;        % trade between real edge vectors, noise vectors, default 0.2
SnakeOptions.Sigma3=1.0;    % sigma used to calculate laplacian in GVF, default 1.0
SnakeOptions.GIterations=200; % Number of GVF iterations, default 0

SnakeOptions.Alpha=0.005;   % membrane energy (1st order), default 0.2 -> amount of stretch in the snake (i.e. tension)
% SnakeOptions.Beta=0.01;   % thin-plate energy (2nd order), default 0.2 -> amount of smoothness in the snake (i.e. rigidity)
SnakeOptions.Beta=0.005;   % thin-plate energy (2nd order), default 0.2 -> amount of smoothness in the snake (i.e. rigidity)
SnakeOptions.Delta=0.01;    % Balloon force, default 0.1
% (Note: alpha and beta will be multiplied with kappa inside the algo)

% % flag of using beta/delta as a matrix (pixel wise beta/delta)
% SnakeOptions.betaMatFlag=0;
% SnakeOptions.deltaMatFlag=0;
%
% % New terms
% % Proposed shape sim term for snake
% SnakeOptions.ShapeSimFlag=1;
% SnakeOptions.ShapeSimCoeffs=[0.1 50];
%
% % Regional Variance Energy (Kabolizade's model)
% SnakeOptions.RVEFlag=0;
% SnakeOptions.RVECoeff=0.5;

%%
tol=50;

%% Run snake and polygonization
% close all
% SnakeOptions
% snake_mask=zeros(size(imgRGB,1),size(imgRGB,2));
snake_mask=zeros(M,N);
snake_seeds={};
snake_results={};
patch_corner=zeros(length(tb_idx_building_in_reg),4);

%%
for idx=1:length(tb_idx_building_in_reg)
    idx_building=tb_idx_building_in_reg(idx);
    
    p_proj=bounds_no_nan{idx_building};
    
    x=[min(p_proj(:,1)) max(p_proj(:,1))];
    y=[min(p_proj(:,2)) max(p_proj(:,2))];
    
    [col,row]=convertGeoCoorMatCoor(ri.XWorldLimits,ri.YWorldLimits,x,y,ri.RasterSize);
    
    patch_corner(idx,:)=[max(1,round(min(col))-tol),max(1,round(min(row))-tol),...
        min(N,round(max(col))+tol),min(M,round(max(row))+tol)];
    
    Z = HR_z_by_grad_rescale(...
        patch_corner(idx,2):patch_corner(idx,4),...
        patch_corner(idx,1):patch_corner(idx,3),:);
    
    % Get initial seeds
    [x,y]=convertGeoCoorMatCoor(ri.XWorldLimits,ri.YWorldLimits,p_proj(:,1),p_proj(:,2),ri.RasterSize);
    P=[y(:)-max(1,min(row)-tol), x(:)-max(1,min(col)-tol)];
    
    if (verbose)
        figure
        imshow(Z)
        hold on
        plot(P(:,2), P(:,1), 'r.', 'MarkerSize', 15)
        % for i=1:size(P,1)
        %     text(P(i,2), P(i,1), num2str(i), 'Color','b')
        % end
    end
    
    if size(p_proj,1)<4 % SRSM requires more than 4 pts in the initial pts
        snake_seeds(idx,:)={P};
        snake_results(idx,:)={P};
        disp('Error: there must be more than 4 pts in the initial pts.')
    else
        SnakeOptions.nPoints=length(p_proj)*2;
        
        % Run the snake model
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
            imshow(Z);
            hold on;
            plot([P(:,2);P(1,2)],[P(:,1);P(1,1)],'c.','LineWidth',2,'MarkerSize',10);
            plot([O(:,2);O(1,2)],[O(:,1);O(1,1)],'r.','LineWidth',2,'MarkerSize',10);
        end
    end
end

%% END
toc


%%
figure
% mapshow(imgRGB,ri)
imshow(HR_z_by_grad_rescale)
hold on
for i=1:length(snake_results)
    s=snake_results{i,:};
    plot(s(:,2)+patch_corner(i,1),s(:,1)+patch_corner(i,2),'r.')
end
save(strcat(dir_data,'res_BE_LiDAR2017_',tile_no,'_New.mat'),'patch_corner','snake_results','snake_seeds','snake_mask')




%% Create shape file from snake results
create_shp_from_snake_result('New',1); % the text is added to discriminate different version

%% Measure accuracy
call_eval_snake;
