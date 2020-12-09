clear
close all
clc

%% Load segment indices and their MBR 
load segs4compare.mat 
% containing ('lidar_segs_in_reg2','mbr_lidar','img_segs_in_reg2_60highest','mbr_geo')

tmp=squeeze([mbr_lidar{lidar_segs_in_reg2}]);
% tmp=squeeze([mbr_geo{img_segs_in_reg2_60highest}]);

lims=[min(tmp(2,:)) max(tmp(2,:));
      min(tmp(1,:)) max(tmp(1,:))];
centers=[mean(tmp(2,:)) mean(tmp(1,:))];

%% Load lidar segments, image segments and the image
load building_boundary_nineth2p

[A, R] = geotiffread('/Users/thanhhuynguyen/Storage/DATA/Orthos_2016/2016_CMQ_243500-5181500.tif');

%% Display inputs
figure
mapshow(rgb2gray(A),R)
hold on

for i=img_segs_in_reg2_60highest
    b2=mbr_geo{i};    
    plot(b2(2,:), b2(1,:), 'g-','LineWidth',1)
end

for i=lidar_segs_in_reg2
    b=mbr_lidar{i};    
    plot(b(2,:), b(1,:), 'r-','LineWidth',1)
end
xlim(lims(1,:))
ylim(lims(2,:))

%% Test 1: compare segment area
area_img_seg=zeros(1,length(img_segs_in_reg2_60highest));
for i=1:length(img_segs_in_reg2_60highest)
    b2=mbr_geo{img_segs_in_reg2_60highest(i)};    
    area_img_seg(i)=polyarea(b2(2,:), b2(1,:));
end

area_pcl_seg=zeros(1,length(lidar_segs_in_reg2));
for i=1:length(lidar_segs_in_reg2)
    b=mbr_lidar{lidar_segs_in_reg2(i)};    
    area_pcl_seg(i)=polyarea(b(2,:), b(1,:));
end

figure
mapshow(rgb2gray(A),R)
hold on

for i=1:length(img_segs_in_reg2_60highest)
    b2=mbr_geo{img_segs_in_reg2_60highest(i)};
    plot(b2(2,:), b2(1,:), 'g-','LineWidth',1)
    text(mean(b2(2,1:4)),mean(b2(1,1:4)),num2str(area_img_seg(i),'%.2f'),...
        'Color','g','FontWeight','bold')
end

for i=1:length(lidar_segs_in_reg2)
    b=mbr_lidar{lidar_segs_in_reg2(i)};
    plot(b(2,:), b(1,:), 'r-','LineWidth',1)
    text(mean(b(2,1:4)),mean(b(1,1:4)),num2str(area_pcl_seg(i),'%.2f'),...
        'Color','r','FontWeight','bold')
end

%% Test 2: segment center

figure
mapshow(rgb2gray(A),R)
hold on

for i=1:length(img_segs_in_reg2_60highest)
    b2=mbr_geo{img_segs_in_reg2_60highest(i)};
    plot(b2(2,:), b2(1,:), 'g-','LineWidth',1)
    plot(mean(b2(2,1:4)),mean(b2(1,1:4)), 'g+')
end

for i=1:length(lidar_segs_in_reg2)
    b=mbr_lidar{lidar_segs_in_reg2(i)};
    plot(b(2,:), b(1,:), 'r-','LineWidth',1)
    plot(mean(b(2,1:4)),mean(b(1,1:4)),'rx')
end
xlim(lims(1,:))
ylim(lims(2,:))

%% GTM algo
P=[];
indices10=[];
k=1;
for i=1:length(img_segs_in_reg2_60highest)
    b2=mbr_geo{img_segs_in_reg2_60highest(i)};
    
    if ~isnan(mean(b2(1,1:4)))
        p=[mean(b2(2,1:4)),mean(b2(1,1:4))]-centers;
        P=[P; p];
        indices10=[indices10, i];
        k=k+1;
    end
end

% first two rows are for centered X and Y coordinates, third row for Z
load('building_mask_nineth2p.mat','building_label','mask_ref','nb_labels','pts_non_gnd_denoised')
ptCloud=pointCloud(pts_non_gnd_denoised(:,1:3));
P2_init=zeros(length(lidar_segs_in_reg2),3);
for i=1:length(lidar_segs_in_reg2)
    b=mbr_lidar{lidar_segs_in_reg2(i)};
    
    P2_init(i,1:2)=[mean(b(2,1:4)),mean(b(1,1:4))]-centers;
    
    pt=[mean(b(2,1:4)),mean(b(1,1:4)), mean(pts_non_gnd_denoised(:,3))];
    [k,~] = findNearestNeighbors(ptCloud,pt,3);
    P2_init(i,3)=mean(ptCloud.Location(k,3));
end

% find initial P from image that are closest to initial P2
P_init=zeros(size(P2_init,1),2);

D=zeros(length(P2_init),length(P));
for i=1:length(P2_init)
    for j=1:length(P)
        D(i,j)=norm(P2_init(i,1:2)-P(j,:));
    end
end

ind=[];
for i=1:length(P2_init)
    [~,id]=min(D(i,:));
    P_init(i,:)=P(id,:);
    ind=[ind, id];
end

%% Display inital matching
line_wd=2.5;
figure
mapshow(rgb2gray(A),R)
hold on
for i=1:length(img_segs_in_reg2_60highest)
    b2=mbr_geo{img_segs_in_reg2_60highest(i)};
    plot(b2(2,:), b2(1,:), 'g-','LineWidth',line_wd)
    plot(mean(b2(2,1:4)),mean(b2(1,1:4)), 'g.')
end

for i=1:length(lidar_segs_in_reg2)
    b=mbr_lidar{lidar_segs_in_reg2(i)};
    plot(b(2,:), b(1,:), 'r-','LineWidth',line_wd)
    plot(mean(b(2,1:4)),mean(b(1,1:4)),'r.')
end
plot(P_init(:,1)+centers(1),P_init(:,2)+centers(2),'r.','LineWidth',line_wd)
plot(P2_init(:,1)+centers(1),P2_init(:,2)+centers(2),'g.','LineWidth',line_wd)
for i=1:length(P_init)
    plot([P_init(i,1) P2_init(i,1)]+centers(1),[P_init(i,2) P2_init(i,2)]+centers(2),'y','LineWidth',line_wd)
end
xlim(lims(1,:))
ylim(lims(2,:))
axis off


%% Run GTM algo
K=5;
[P_final, P2_final, MKNNG_0, MKNNG2_0, MKNNG, MKNNG2] = GTM_algo_original(P_init,P2_init(:,1:2), K);

indices1=zeros(1,length(P_final));
for i=1:length(P_final)
    id=find(P(:,1)==P_final(i,1) & P(:,2)==P_final(i,2));
    indices1(i)=indices10(id);
end

indices2=zeros(1,length(P2_final));
for i=1:length(P2_final)
    id=find(P2_init(:,1)==P2_final(i,1) & P2_init(:,2)==P2_final(i,2));
    indices2(i)=id;
end

P2_final=[P2_final P2_init(indices2,3)];

% Display GTM result
figure(1)
hold on
plot(P_final(:,1),P_final(:,2),'mx')
plot(P2_final(:,1),P2_final(:,2),'b+')
grid on

%% Display GTM matching
line_wd=2.5;
figure
mapshow(rgb2gray(A),R)
hold on
for i=indices1%1:length(img_segs_in_reg2_60highest)
    b2=mbr_geo{img_segs_in_reg2_60highest(i)};
    plot(b2(2,:), b2(1,:), 'g-','LineWidth',line_wd)
    plot(mean(b2(2,1:4)),mean(b2(1,1:4)), 'g.')
end

for i=indices2%1:length(lidar_segs_in_reg2)
    b=mbr_lidar{lidar_segs_in_reg2(i)};
    plot(b(2,:), b(1,:), 'r-','LineWidth',line_wd)
    plot(mean(b(2,1:4)),mean(b(1,1:4)),'r.')
end
% plot(P_init(:,1)+centers(1),P_init(:,2)+centers(2),'r.','LineWidth',line_wd)
% plot(P2_init(:,1)+centers(1),P2_init(:,2)+centers(2),'g.','LineWidth',line_wd)
for i=1:length(P_final)
    plot([P_final(i,1) P2_final(i,1)]+centers(1),[P_final(i,2) P2_final(i,2)]+centers(2),'y','LineWidth',line_wd)
end
xlim(lims(1,:))
ylim(lims(2,:))
axis off


%% Try RANSAC
[fRANSAC,inliersIndex,status]= estimateFundamentalMatrix(P_init,P2_init(:,1:2),'Method','RANSAC',...
    'DistanceThreshold',1e-4);

P_final_RANSAC=P_init(inliersIndex,:);
P2_final_RANSAC=P2_init(inliersIndex,:);

indices1R=zeros(1,length(P_final_RANSAC));
for i=1:length(P_final_RANSAC)
    id=find(P(:,1)==P_final_RANSAC(i,1) & P(:,2)==P_final_RANSAC(i,2));
    indices1R(i)=indices10(id);
end

indices2R=zeros(1,length(P2_final_RANSAC));
for i=1:length(P2_final_RANSAC)
    id=find(P2_init(:,1)==P2_final_RANSAC(i,1) & P2_init(:,2)==P2_final_RANSAC(i,2));
    indices2R(i)=id;
end

%% Display RANSAC matching
line_wd=2.5;
figure
mapshow(rgb2gray(A),R)
hold on
for i=indices1R%1:length(img_segs_in_reg2_60highest)
    b2=mbr_geo{img_segs_in_reg2_60highest(i)};
    plot(b2(2,:), b2(1,:), 'g-','LineWidth',line_wd)
    plot(mean(b2(2,1:4)),mean(b2(1,1:4)), 'g.')
end

for i=indices2R%1:length(lidar_segs_in_reg2)
    b=mbr_lidar{lidar_segs_in_reg2(i)};
    plot(b(2,:), b(1,:), 'r-','LineWidth',line_wd)
    plot(mean(b(2,1:4)),mean(b(1,1:4)),'r.')
end
% plot(P_init(:,1)+centers(1),P_init(:,2)+centers(2),'r.','LineWidth',line_wd)
% plot(P2_init(:,1)+centers(1),P2_init(:,2)+centers(2),'g.','LineWidth',line_wd)
for i=1:length(P_final_RANSAC)
    plot([P_final_RANSAC(i,1) P2_final_RANSAC(i,1)]+centers(1),[P_final_RANSAC(i,2) P2_final_RANSAC(i,2)]+centers(2),'y','LineWidth',line_wd)
end
xlim(lims(1,:))
ylim(lims(2,:))
axis off

%%

% Display matched segments

% before matching
figure(98)
mapshow(rgb2gray(A),R)
hold on
for i=indices10
    b2=mbr_geo{img_segs_in_reg2_60highest(i)};
    plot(b2(2,:), b2(1,:), 'g-','LineWidth',2)  
end
for i=1:length(lidar_segs_in_reg2)
    b=mbr_lidar{lidar_segs_in_reg2(i)};
    plot(b(2,:), b(1,:), 'r-','LineWidth',2)
end
plot(P_init(:,1)+centers(1),P_init(:,2)+centers(2), 'g+', 'MarkerSize',2)
plot(P2_init(:,1)+centers(1),P2_init(:,2)+centers(2), 'rx', 'MarkerSize',2)
for i=1:length(P2_init)
    plot([P_init(i,1)+centers(1) P2_init(i,1)+centers(1)],...
         [P_init(i,2)+centers(2) P2_init(i,2)+centers(2)], 'y-','LineWidth',2)
end
xlim(lims(1,:))
ylim(lims(2,:))
% title('Before matching')

% after matching
figure(99)
mapshow(rgb2gray(A),R)
hold on
for i=indices1
    b2=mbr_geo{img_segs_in_reg2_60highest(i)};
    plot(b2(2,:), b2(1,:), 'g-','LineWidth',2)  
end
for i=indices2
    b=mbr_lidar{lidar_segs_in_reg2(i)};
    plot(b(2,:), b(1,:), 'r-','LineWidth',2)
end
plot(P_final(:,1)+centers(1),P_final(:,2)+centers(2), 'g+', 'MarkerSize',2)
plot(P2_final(:,1)+centers(1),P2_final(:,2)+centers(2), 'rx', 'MarkerSize',2)
for i=1:length(P2_final)
    plot([P_final(i,1)+centers(1) P2_final(i,1)+centers(1)],...
         [P_final(i,2)+centers(2) P2_final(i,2)+centers(2)], 'y-','LineWidth',2)
end
xlim(lims(1,:))
ylim(lims(2,:))
% title('Result of segment center matching using GTM')

%% calculate areas

area_img_seg=zeros(1,length(img_segs_in_reg2_60highest));
for i=1:length(img_segs_in_reg2_60highest)
    b2=mbr_geo{img_segs_in_reg2_60highest(i)};    
    area_img_seg(i)=polyarea(b2(2,:), b2(1,:));
end
area_pcl_seg=zeros(1,length(lidar_segs_in_reg2));
for i=1:length(lidar_segs_in_reg2)
    b=mbr_lidar{lidar_segs_in_reg2(i)};    
    area_pcl_seg(i)=polyarea(b(2,:), b(1,:));
end


%% saving matching result for the Transformation Model Estimation
% Consider pairwise segment area and orientation
% idx_area=[];
idx_ori=[];
idx=[];
ori_img=zeros(1,length(indices1));
ori_pcl=zeros(1,length(indices2));
for i=1:length(indices1)
    b_img=mbr_geo{img_segs_in_reg2_60highest(indices1(i))};
    b_pcl=mbr_lidar{lidar_segs_in_reg2(indices2(i))};
    
    % segment area not differ more than 0.3
%     if abs(area_img_seg(indices1(i))-area_pcl_seg(indices2(i)))/area_pcl_seg(indices2(i))<0.3
%         idx_area=[idx_area, i];
%     end
    
    
    if norm(b_img(:,1)-b_img(:,2))>norm(b_img(:,3)-b_img(:,2))
        ori_img(i)=atand((b_img(1,2)-b_img(1,1))/(b_img(2,2)-b_img(2,1)));
    else
        ori_img(i)=atand((b_img(1,3)-b_img(1,2))/(b_img(2,3)-b_img(2,2)));
    end
    if norm(b_pcl(:,1)-b_pcl(:,2))>norm(b_pcl(:,3)-b_pcl(:,2))
        ori_pcl(i)=atand((b_pcl(1,2)-b_pcl(1,1))/(b_pcl(2,2)-b_pcl(2,1)));
    else
        ori_pcl(i)=atand((b_pcl(1,3)-b_pcl(1,2))/(b_pcl(2,3)-b_pcl(2,2)));
    end
    
%     if abs(ori_img(i)-ori_pcl(i))/ori_pcl(i)<0.3
%         idx_ori=[idx_ori, i];
%     end

    % segment orientation not differ more than 0.3 && segment area not differ more than 0.3
    if abs(abs(ori_img(i)-ori_pcl(i))/ori_pcl(i))<0.4 && ...
            abs(area_img_seg(indices1(i))-area_pcl_seg(indices2(i)))/area_pcl_seg(indices2(i))<0.4
        idx=[idx, i];
    end
end

figure(10)
mapshow(rgb2gray(A),R)
hold on
for i=indices1(idx)
    b2=mbr_geo{img_segs_in_reg2_60highest(i)};
    p1=plot(b2(2,:), b2(1,:), 'g-','LineWidth',2);
end
for i=indices2(idx)
    b=mbr_lidar{lidar_segs_in_reg2(i)};
    p2=plot(b(2,:), b(1,:), 'r-','LineWidth',2);
end
p3=plot(P_final(idx,1)+centers(1),P_final(idx,2)+centers(2), 'g+', 'MarkerSize',10);
p4=plot(P2_final(idx,1)+centers(1),P2_final(idx,2)+centers(2), 'rx', 'MarkerSize',10);
for i=idx
    plot([P_final(i,1)+centers(1) P2_final(i,1)+centers(1)],...
         [P_final(i,2)+centers(2) P2_final(i,2)+centers(2)], 'y-','LineWidth',2)
%     text(P_final(i,1)+centers(1),P_final(i,2)+centers(2),num2str(i),'Color','y','FontWeight','bold','FontSize',14)
end
xlim(lims(1,:))
ylim(lims(2,:))
legend([p1,p2,p3,p4],'Image building segments','LiDAR building segments','2D correspondences','3D correspondences')

point2d=[P_final(idx,1)+centers(1),P_final(idx,2)+centers(2)];
point3d=[P2_final(idx,1)+centers(1),P2_final(idx,2)+centers(2),P2_final(idx,3)];
% save('corr4TME_reg2_byGTM.mat','point2d','point3d')


%% Final result

line_wd=2.5;
figure
mapshow(rgb2gray(A),R)
hold on
for i=indices1(idx)%1:length(img_segs_in_reg2_60highest)
    b2=mbr_geo{img_segs_in_reg2_60highest(i)};
    plot(b2(2,:), b2(1,:), 'g-','LineWidth',line_wd)
    plot(mean(b2(2,1:4)),mean(b2(1,1:4)), 'g.')
end

for i=indices2(idx)%1:length(lidar_segs_in_reg2)
    b=mbr_lidar{lidar_segs_in_reg2(i)};
    plot(b(2,:), b(1,:), 'r-','LineWidth',line_wd)
    plot(mean(b(2,1:4)),mean(b(1,1:4)),'r.')
end
% plot(P_init(:,1)+centers(1),P_init(:,2)+centers(2),'r.','LineWidth',line_wd)
% plot(P2_init(:,1)+centers(1),P2_init(:,2)+centers(2),'g.','LineWidth',line_wd)
for i=idx
    plot([P_final(i,1)+centers(1) P2_final(i,1)+centers(1)],...
         [P_final(i,2)+centers(2) P2_final(i,2)+centers(2)], 'y-','LineWidth',line_wd)
%     text(P_final(i,1)+centers(1),P_final(i,2)+centers(2),num2str(i),'Color','y','FontWeight','bold','FontSize',14)
end
xlim(lims(1,:))
ylim(lims(2,:))
axis off


%% Display matched segments by RANSAC
figure(100)
mapshow(rgb2gray(A),R)
hold on
for i=indices1R
    b2=mbr_geo{img_segs_in_reg2_60highest(i)};
    plot(b2(2,:), b2(1,:), 'g-','LineWidth',2)
end
for i=indices2R
    b=mbr_lidar{lidar_segs_in_reg2(i)};
    plot(b(2,:), b(1,:), 'r-','LineWidth',2) 
end
plot(P_final_RANSAC(:,1)+centers(1),P_final_RANSAC(:,2)+centers(2), 'g+', 'MarkerSize',2)
plot(P2_final_RANSAC(:,1)+centers(1),P2_final_RANSAC(:,2)+centers(2), 'rx', 'MarkerSize',2)
for i=1:length(P2_final_RANSAC)
    plot([P_final_RANSAC(i,1)+centers(1) P2_final_RANSAC(i,1)+centers(1)],...
         [P_final_RANSAC(i,2)+centers(2) P2_final_RANSAC(i,2)+centers(2)], 'y-','LineWidth',2)
end
xlim(lims(1,:))
ylim(lims(2,:))
title('Result of segment center matching using RANSAC')