clear
close all

filename = '/Users/thanhhuynguyen/Storage/DATA/Orthos_2016/2016_CMQ_243500-5181500.tif';
[A, R] = geotiffread(filename);

% 1st row: residential, 2nd: avg-sized buildings, 3rd: big and complex buildings
row=[1 2000; 1 1550; 3300 5700];
col=[1 2000; 2350 4000; 1500 3500];

% convert into geospatial coordinates
[x,y]=convertMatCoorGeoCoor(R.XWorldLimits,R.YWorldLimits,col(:),row(:),R.RasterSize);
x=reshape(x,[3 2]);
y=reshape(y,[3 2]);

%% Load meanshift result
matfile = 'reg1_bw2.5_0.mat';
load(matfile)

seg_npts_all = seg_area_all;
clear seg_area_all
nb_segs = size(bound_all,1);

s=strsplit(matfile,'_');
if strcmp(s{1},'reg1')
    iReg=1;
elseif strcmp(s{1},'reg2')
    iReg=2;
else
    iReg=3;
end

I = A(row(iReg,1):row(iReg,2), col(iReg,1):col(iReg,2), :);

load reg1_res.mat

%% Get MBR and calculate their fill percentage - FROM IMAGE
segment_thres=[1000 2e4];

k=1;
tb_segs_inrange=[];
mbr_area=[];
mbr_orie=[];
seg_area=[]; 
seg_area_from_mask=[];
mbr_fill_percent=[];
mbr_fill_percent2=[];

for i=1:nb_segs
    if seg_npts_all(i)>segment_thres(1) && seg_npts_all(i)<segment_thres(2)
        
        tb_segs_inrange=[tb_segs_inrange i];
       
        % get min bounding box
        bb2=mbr{k,:};

        % calculate bounding box area
        mbr_area = [mbr_area, polyarea(bb2(1,:), bb2(2,:))];
        
        % calculate segment area
        b=bound_all{i};
        seg_area=[seg_area, polyarea(b(:,2),b(:,1))];
        seg_area_from_mask=[seg_area_from_mask size(mask_all{i},1)];
         
%         % calculate the % of filling of each segment inside its bounding box
        mbr_fill_percent=[mbr_fill_percent, seg_area(k)/mbr_area(k)];
        mbr_fill_percent2=[mbr_fill_percent2, seg_area_from_mask(k)/mbr_area(k)];
        
        % calculate orientation of the box
        axe1=norm(bb2(:,2)-bb2(:,1));
        axe2=norm(bb2(:,2)-bb2(:,3));
        if axe2>axe1
            orie=atan((bb2(2,3)-bb2(2,2))/(bb2(1,3)-bb2(1,2)));
        else
            orie=atan((bb2(2,2)-bb2(2,1))/(bb2(1,2)-bb2(1,1)));
        end
        mbr_orie=[mbr_orie orie];
        
        k=k+1;
    end
end

nb_segs_inrange=size(tb_segs_inrange,2);


n=100;
[~,idx]=maxk(mbr_fill_percent,n);
figure(99)
subplot(121)
imshow(rgb2gray(I))
hold on
for i=idx
%     p=mask_all{tb_segs_inrange(i),:};
%     plot(p(:,2),p(:,1),'.','Color',rand(3,1));

    p=bound_all{tb_segs_inrange(i),:};
    plot(p(:,2),p(:,1),'g.');
    
    bb=mbr{i};
    plot(bb(1,:), bb(2,:), 'r-','LineWidth',1)
end
title([num2str(n) ' segments that have highest MBR filling percentage'])


%% Get MBR and calculate their fill percentage - FROM POINT CLOUD
load building_boundary_nineth2
load('building_mask_nineth2.mat','building_label','mask_ref','nb_labels','pts_non_gnd_denoised')

mbr_lidar={};
mbr_lidar_area=zeros(1,nb_labels);
seg_lidar_area=zeros(1,nb_labels);

for i=1:nb_labels
    pt=points{i,:};
    b=bounds{i,:};
    
    bb3 = minBoundingBox([pt(:,2)'; pt(:,1)']);
    bb4=[bb3 bb3(:,1)];
    mbr_lidar(i,:)={bb4};
    
    mbr_lidar_area(i)=polyarea(bb4(1,:), bb4(2,:));
    seg_lidar_area(i)=polyarea(b(:,2),b(:,1)); 
end
mbr_lidar_fill_percent=seg_lidar_area./mbr_lidar_area;

lidar_segs_in_reg1=[]; % lidar segments in the region, to be used afterwards

n=100;
[~,idx]=maxk(mbr_lidar_fill_percent,n);

figure(99)
subplot(122)
mapshow(rgb2gray(I),R)
hold on
for i=idx
    b=bounds{i,:};
    plot(b(:,1),b(:,2),'.','Color',rand(3,1));

    bb4=mbr_lidar{i};
    plot(bb4(2,:), bb4(1,:), 'r-','LineWidth',1)
    
    center=[mean(b(:,1)), mean(b(:,2))];
    if (center(1)>x(iReg,1) && center(1)<x(iReg,2)) &&...
            (center(2)>y(iReg,2) && center(2)<y(iReg,1))
        lidar_segs_in_reg1=[lidar_segs_in_reg1 i];
    end
end
title([num2str(n) ' segments that have highest MBR filling percentage'])

xlim([x(iReg,1) x(iReg,2)])
ylim([y(iReg,2) y(iReg,1)])


%% Calculate area of segment (in distance units) from image

n=75;
[~,idx]=maxk(mbr_fill_percent,n);
img_segs_in_reg1_75highest=idx;


figure(100)
mapshow(rgb2gray(A),R)
hold on
for i=idx
    bb2=mbr{i,:};
    [bbx,bby]=convertMatCoorGeoCoor(R.XWorldLimits,R.YWorldLimits,...
        bb2(1,:)+col(iReg,1),bb2(2,:)+row(iReg,1),R.RasterSize);
    p1=plot(bbx, bby, 'y-','LineWidth',1.5);
    
    center=[mean(bbx(1:4)), mean(bby(1:4))];
%     text(center(1),center(2),num2str(polyarea(bbx, bby),'%.2f'),'Color','g')
    text(center(1),center(2),num2str(mbr_orie(i),'%.2f'),'Color','g')
end
xlim([x(1,1) x(1,2)])
ylim([y(1,2) y(1,1)])

%% Calculate area of region from lidar pointcloud

figure(101)
% pcshow(pts_non_gnd_denoised(:,1:3))
mapshow(rgb2gray(A),R)
hold on
for i=1:nb_labels
    b=points{i,:};
    bb4=mbr_lidar{i};
    
    plot(bb4(2,:), bb4(1,:), 'r-','LineWidth',1)

    center=[mean(b(:,1)), mean(b(:,2))];
%     text(center(1),center(2),num2str(polyarea(b(:,1), b(:,2)),'%.2f'),'Color','g')
    text(center(1),center(2),num2str(polyarea(bb4(1,:), bb4(2,:)),'%.2f'),'Color','c')
end
grid on
view(2)
xlim([x(iReg,1) x(iReg,2)])
ylim([y(iReg,2) y(iReg,1)])

%% Display both segmentation results from lidar and image

figure(102)
mapshow(rgb2gray(A),R)
hold on

% from image, 100 highest MBR filling percentage segments
n=150;
[~,idx]=maxk(mbr_fill_percent,n);

for i=idx
    bb2=mbr{i,:};
    [bbx,bby]=convertMatCoorGeoCoor(R.XWorldLimits,R.YWorldLimits,...
        bb2(1,:)+col(iReg,1),bb2(2,:)+row(iReg,1),R.RasterSize);
    p1=plot(bbx, bby, 'g-','LineWidth',1.5);
    
%     center=[mean(bbx(1:4)), mean(bby(1:4))];
%     text(center(1),center(2),num2str(polyarea(bbx, bby),'%.2f'),'Color','g')
end

% from lidar point cloud
% [~,idx]=maxk(mbr_lidar_fill_percent,n);

for i=1:nb_labels    
    bb4=mbr_lidar{i};
    p2=plot(bb4(2,:), bb4(1,:), 'r-','LineWidth',1.5);

%     center=[mean(b(:,1)), mean(b(:,2))];
%     text(center(1),center(2),num2str(polyarea(bb4(1,:), bb4(2,:)),'%.2f'),'Color','c')
end
view(2)
xlim([x(iReg,1) x(iReg,2)])
ylim([y(iReg,2) y(iReg,1)])
legend([p1,p2],'from optical image','from LiDAR pointcloud')

%% Save segments id for later registration dev
mbr_geo={};
for i=1:nb_segs_inrange
    bb2=mbr{i,:};
    
    [bbx,bby]=convertMatCoorGeoCoor(R.XWorldLimits,R.YWorldLimits,...
        bb2(1,:)+col(iReg,1),bb2(2,:)+row(iReg,1),R.RasterSize);
    mbr_geo(i,:)={[bby; bbx]};
end

save('segs4compare_reg1.mat','lidar_segs_in_reg1','mbr_lidar','img_segs_in_reg1_75highest','mbr_geo')
