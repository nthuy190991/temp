clear
close all
clc

dir='/Users/thanhhuynguyen/Dropbox/Thesis/dev/';
dir2='/Users/thanhhuynguyen/Storage/DATA/';

% dir='C:\Users\thnguyen\Dropbox\Thesis\dev\';
% dir2='D:\Data\';
% cd(strcat(dir,'Snakes'))

verbose=0;
saveFlag=0;

snake_seeds={}; % seeds for snake model (i.e. lidar projected points)
snake_results={}; % result of snake model
building_vertices={}; % result of polygonization
% lidar_mbr={}; % lidar points' MBR used for the polygonization

% [imgRGB, R] = geotiffread(strcat(dir2,'/Orthos_2016/2016_CMQ_243500-5181500.tif'));

row=[9800 11950; 14300 16450; 15820 18320];
col=[9200 10570; 6790 8640; 8800 10450];



load(strcat(dir,'building_boundary_area1'))
% load(strcat(dir,'building_boundary_nineth2p'))
tb_idx_building=1:length(bounds);%[151 195 202 203 207];

iReg=1;
[pts,imgRGB,ri,rx,ry]=load_datasets_isprs(iReg,1,0);
%%

xRegSel=[min(rx) max(rx)];
yRegSel=[min(ry) max(ry)];


tb_idx_building_in_reg=[];
for idx=1:length(tb_idx_building)
    idx_building=tb_idx_building(idx);
    
    b=bounds{idx_building};
    if all(b(:,1)>min(xRegSel) & b(:,1)<max(xRegSel) & b(:,2)>min(yRegSel) & b(:,2)<max(yRegSel))
        tb_idx_building_in_reg=[tb_idx_building_in_reg idx_building];
    end
end

% remove buildings outside of area 1
tb_idx_building_in_reg=setdiff(tb_idx_building_in_reg,[13 19:21 24 26]);

figure
mapshow(rgb2gray(imgRGB),ri)
hold on
for idx=tb_idx_building_in_reg
    b=points{idx};
    b=b(boundary(b(:,1:2),0.75),:);
    plot(b(:,1),b(:,2),'-')
    text(mean(b(:,1)),mean(b(:,2)),num2str(idx),'Color','r')
end
xlim([min(xRegSel) max(xRegSel)])
ylim([min(yRegSel) max(yRegSel)])

%% Snake option
SnakeOptions=struct;
SnakeOptions.Verbose=false;

SnakeOptions.Gamma=1;        % time step, default 1
SnakeOptions.Iterations=150; % number of iterations, default 100

Options.Wline=0.04;    % attraction to line, default 0.04
Options.Wedge=2;        % attraction to edge, default 2
Options.Wterm=0.01;      % attraction to terminations of lines and corners, default 0.01
SnakeOptions.Sigma1=10;       % sigma used to calculate image derivative, default 10
SnakeOptions.Sigma2=20;       % sigma used to calculate the gradient of the edge energy image (which gives the image force), default 20
SnakeOptions.Kappa=0.25;       % weight of external IMAGE force, default 2

% GVF Params
SnakeOptions.Mu=0.2;         % trade between real edge vectors, noise vectors, default 0.2
SnakeOptions.Sigma3=1.0;     % sigma used to calculate laplacian in GVF, default 1.0
SnakeOptions.GIterations=200; % Number of GVF iterations, default 0

SnakeOptions.Alpha=0.005;     % membrane energy (1st order), default 0.2 -> amount of stretch in the snake (i.e. tension)
SnakeOptions.Beta=0.005;      % thin-plate energy (2nd order), default 0.2 -> amount of smoothness in the snake (i.e. rigidity)
SnakeOptions.Delta=0.025;       % Balloon force, default 0.1

% Regional Variance Energy (kabolizade's model) // Unused atm
SnakeOptions.RVEflag=0;
SnakeOptions.RVEcoeff=0.025;

tol=50;

%%
cd Snakes

%% Polygonization parameters
Vspec=10;
Mspec=50;

%% Run snake and polygonization

snake_mask=zeros(size(imgRGB,1),size(imgRGB,2));

patch_corner=zeros(length(tb_idx_building_in_reg),4);
for idx=1:length(tb_idx_building_in_reg)
    %     idx_building=1
    idx_building=tb_idx_building_in_reg(idx);
    
    
    %     p_proj=bounds{idx_building};
    p=points{idx_building};
    p_proj=p(boundary(p(:,1:2),1),:);
    
    % Altitude thresholding (remove points on vertical surface)
    %     h=histogram(p(:,3));
    %     [h,edges]=histcounts(p(:,3),100);
    %     hist_cum=cumsum(h);
    %     idx_hist_80percent=min(find(hist_cum./hist_cum(end)>0.15));
    %     thres=edges(idx_hist_80percent-1);
    
    %     p_proj=p(p(:,3)>=thres,:);
    
    x=[min(p_proj(:,1)) max(p_proj(:,1))];
    y=[min(p_proj(:,2)) max(p_proj(:,2))];
    
    [col,row]=convertGeoCoorMatCoor(ri.XWorldLimits,ri.YWorldLimits,x,y,ri.RasterSize);
    
    patch_corner(idx,:)=[max(0,min(col)-tol),max(0,min(row)-tol),...
        min(size(imgRGB,2),max(col)+tol),...
        min(size(imgRGB,1),max(row)+tol)];
    
    I = im2double(rgb2gray(imgRGB(patch_corner(idx,2):patch_corner(idx,4),...
        patch_corner(idx,1):patch_corner(idx,3),:)));
    
    % Get initial seeds
    %     b2=p_proj(boundary(p_proj(:,1),p_proj(:,2),0.5),:);
    %     [x,y]=convertGeoCoorMatCoor(ri.XWorldLimits,ri.YWorldLimits,b2(:,1),b2(:,2),ri.RasterSize);
    [x,y]=convertGeoCoorMatCoor(ri.XWorldLimits,ri.YWorldLimits,p_proj(:,1),p_proj(:,2),ri.RasterSize);
    P=[y(:)-min(row)+tol, x(:)-min(col)+tol];
    
    if (verbose)
        figure
        imshow(I)
        hold on
        plot(P(:,2), P(:,1), 'r.', 'MarkerSize', 15)
        % for i=1:size(P,1)
        %     text(P(i,2), P(i,1), num2str(i), 'Color','b')
        % end
    end
    
    if (size(I,1)*size(I,2)>1e5)
        SnakeOptions.nPoints=200;    % number of contour points, default 100
    else
        SnakeOptions.nPoints=100;    % number of contour points, default 100
        %         disp(size(I))
    end
    
    % Proposed new term for snake
    SnakeOptions.DFflag=1;
    SnakeOptions.DFseeds=P;
    SnakeOptions.DFcoeffs=[0.01 100];
    
    % Run the snake model
    % SnakeOptions.Verbose=1;
    [O,J]=Snake2D(I,P,SnakeOptions);
    
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

% save('BE_by_snake_PIA19_isprs_area1.mat','tb_idx_building_in_reg','bounds','snake_seeds','snake_results','snake_mask','patch_corner','building_vertices')

iou_snake=sum(sum(snake_mask&ref_mask))/sum(sum(snake_mask|ref_mask));
com_snake=sum(sum(snake_mask&ref_mask))/sum(sum(ref_mask));
cor_snake=sum(sum(snake_mask&ref_mask))/sum(sum(snake_mask));

display([num2str(iou_snake) ', ' num2str(com_snake) ', ' num2str(cor_snake)])


%%

figure
imshow(imgRGB)

for i=1:length(snake_seeds)
    s=snake_seeds{i};
    r=snake_results{i};
    p=building_vertices{i};
    hold on
    plot(s(:,2)+patch_corner(i,1),s(:,1)+patch_corner(i,2),'b.')
    plot(r(:,2)+patch_corner(i,1),r(:,1)+patch_corner(i,2),'g.')
    plot([p(2,:) p(2,1)]+patch_corner(i,1),[p(1,:) p(1,1)]+patch_corner(i,2),'r-')
end

figure
imshow(snake_mask)

figure
mapshow(rgb2gray(imgRGB),ri)
hold on
for i=1:length(snake_seeds)
    s=snake_seeds{i};
    r=snake_results{i};
    
    [x,y]=convertMatCoorGeoCoor(ri.XWorldLimits,ri.YWorldLimits,...
        s(:,2)+patch_corner(i,1),...
        s(:,1)+patch_corner(i,2),...
        ri.RasterSize);
    
    plot(x,y,'b.')
    
    [x,y]=convertMatCoorGeoCoor(ri.XWorldLimits,ri.YWorldLimits,...
        r(:,2)+patch_corner(i,1),...
        r(:,1)+patch_corner(i,2),...
        ri.RasterSize);
    %     plot(x,y,'g.')
end

%% Evaluation based on ground truth

I=rgb2gray(imgRGB);

filename=strcat(dir2,'ISPRS_Benchmark/Reference_3d_reconstruction/Reference_Buildings/');
[S,A]=shaperead(strcat(filename,'building_outline_area_1.shp'));

figure(100)
mapshow(I,ri)
hold on
grid on
for i=1:length(S)
    plot(S(i).X,S(i).Y,'r','LineWidth',3)
end
% for i=1:length(S)
%     plot([S(i).BoundingBox(1,1) S(i).BoundingBox(1,1)],...
%          [S(i).BoundingBox(1,2) S(i).BoundingBox(2,2)],'g')
%     plot([S(i).BoundingBox(2,1) S(i).BoundingBox(2,1)],...
%          [S(i).BoundingBox(1,2) S(i).BoundingBox(2,2)],'g')
%     plot([S(i).BoundingBox(1,1) S(i).BoundingBox(2,1)],...
%          [S(i).BoundingBox(1,2) S(i).BoundingBox(1,2)],'g')
%     plot([S(i).BoundingBox(1,1) S(i).BoundingBox(2,1)],...
%          [S(i).BoundingBox(2,2) S(i).BoundingBox(2,2)],'g')
% end
for i=1:length(snake_seeds)
    s=snake_seeds{i};
    r=snake_results{i};
    
    [x,y]=convertMatCoorGeoCoor(ri.XWorldLimits,ri.YWorldLimits,...
        s(:,2)+patch_corner(i,1),...
        s(:,1)+patch_corner(i,2),...
        ri.RasterSize);
    
    plot(x,y,'b.')
    
    [x,y]=convertMatCoorGeoCoor(ri.XWorldLimits,ri.YWorldLimits,...
        r(:,2)+patch_corner(i,1),...
        r(:,1)+patch_corner(i,2),...
        ri.RasterSize);
    plot(x,y,'g.')
end
axis equal
xlim([min(rx) max(rx)])
ylim([min(ry) max(ry)])


figure
imshow(I)
hold on
for i=1:length(S)
    [x,y]=convertGeoCoorMatCoor(ri.XWorldLimits,ri.YWorldLimits,...
        S(i).X(1:end-1),S(i).Y(1:end-1),ri.RasterSize);
    if ~any(isnan(x))
        plot(x,y,'r')
    end
end

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

load BE_by_snake_PIA19_isprs_area1

% Snake result (w/o polygonization)
iou_snake=sum(sum(snake_mask&ref_mask))/sum(sum(snake_mask|ref_mask))
com_snake=sum(sum(snake_mask&ref_mask))/sum(sum(ref_mask))
cor_snake=sum(sum(snake_mask&ref_mask))/sum(sum(snake_mask))

ref_mask_color=zeros(size(I,1),size(I,2));
ref_mask_color(ref_mask==1&snake_mask==1)=1;
ref_mask_color(ref_mask==0&snake_mask==1)=2;
ref_mask_color(ref_mask==1&snake_mask==0)=3;

figure
imagesc(ref_mask_color)
hold on
[x,y]=find(ref_mask_color==0);
plot(y,x,'w.')
[x,y]=find(ref_mask_color==1);
p1=plot(y,x,'y.');
[x,y]=find(ref_mask_color==2);
p2=plot(y,x,'r.');
[x,y]=find(ref_mask_color==3);
p3=plot(y,x,'b.');
axis equal
axis off
legend([p1,p2,p3],{'True Positive pixels','False Positive pixels','False Negative pixels'})

% Snake result (with polygonization)
snake_mask_poly=zeros(size(I,1),size(I,2));
for i=1:length(building_vertices)
    s=building_vertices{i};
    inBox = reshape(inpolygon(pt(:,1),pt(:,2),...
        [s(1,:) s(1,1)]+patch_corner(i,2),...
        [s(2,:) s(2,1)]+patch_corner(i,1)),size(I,2),size(I,1));
    snake_mask_poly = snake_mask_poly | inBox';
end

iou_snake_poly=sum(sum(snake_mask_poly&ref_mask))/sum(sum(snake_mask_poly|ref_mask))
com_snake_poly=sum(sum(snake_mask_poly&ref_mask))/sum(sum(ref_mask))
cor_snake_poly=sum(sum(snake_mask_poly&ref_mask))/sum(sum(snake_mask_poly))

% Snake result (with QGIS polygonization) 03 July 2019

[SS,AA]=shaperead('../test_simp_tol05.shp');
snake_mask_poly2=zeros(size(I,1),size(I,2));
for i=1:length(SS)
    if (mean(SS(i).X(1:end-1))<=ri.XWorldLimits(2) && mean(SS(i).X(1:end-1))>=ri.XWorldLimits(1) &&...
            mean(SS(i).Y(1:end-1))<=ri.YWorldLimits(2) && mean(SS(i).Y(1:end-1))>=ri.YWorldLimits(1))
        
        [x,y]=convertGeoCoorMatCoor(ri.XWorldLimits,ri.YWorldLimits,...
            SS(i).X(1:end-1),SS(i).Y(1:end-1),ri.RasterSize);
        inBox = reshape(inpolygon(pt(:,1),pt(:,2),y,x),size(I,2),size(I,1));
        snake_mask_poly2 = snake_mask_poly2 | inBox';
    end
end

iou_snake_poly=sum(sum(snake_mask_poly2&ref_mask))/sum(sum(snake_mask_poly2|ref_mask))
com_snake_poly=sum(sum(snake_mask_poly2&ref_mask))/sum(sum(ref_mask))
cor_snake_poly=sum(sum(snake_mask_poly2&ref_mask))/sum(sum(snake_mask_poly2))

mask_color=zeros(size(I,1),size(I,2));
mask_color(ref_mask==1&snake_mask_poly2==1)=1;
mask_color(ref_mask==0&snake_mask_poly2==1)=2;
mask_color(ref_mask==1&snake_mask_poly2==0)=3;
figure
imagesc(mask_color)

% LiDAR boundary result
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